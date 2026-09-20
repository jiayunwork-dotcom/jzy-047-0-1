"""FastAPI 入口：装配校验、Lamé 计算与 HTTP 接口。

服务无任何持久化与共享可变状态：每个请求独立构造自己的
LameCylinder 并纯函数式地完成计算，并发请求之间互不影响。
"""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.encoders import jsonable_encoder
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from . import __version__
from .lame import LameCylinder
from .models import (
    ProfileRequest,
    ProfileResponse,
    StressPoint,
    StressRequest,
    StressResponse,
)
from .thinwall import thin_wall_hoop_stress
from .validation import (
    CalculationError,
    validate_geometry,
    validate_poisson_ratio,
    validate_pressures,
    validate_query_radius,
    validate_yield_strength,
)

app = FastAPI(
    title="厚壁圆筒 Lamé 应力分析服务",
    version=__version__,
    description=(
        "给定内/外半径、内/外压力与轴向约束方式，按 Lamé 经典弹性解 "
        "σr = A − B/r²、σθ = A + B/r² 计算厚壁圆筒应力场。"
        "应力符号约定：拉为正；压力以压为正。"
    ),
)


@app.exception_handler(CalculationError)
async def calculation_error_handler(
    _request: Request, exc: CalculationError
) -> JSONResponse:
    """领域校验失败 → 400 + 结构化错误原因。"""
    return JSONResponse(
        status_code=400,
        content={"error": {"code": exc.code, "message": exc.message}},
    )


@app.exception_handler(RequestValidationError)
async def request_validation_error_handler(
    _request: Request, exc: RequestValidationError
) -> JSONResponse:
    """schema 校验失败（如压力未给全、类型错误）→ 422 + 结构化错误。"""
    return JSONResponse(
        status_code=422,
        content={
            "error": {
                "code": "INVALID_REQUEST",
                "message": "请求参数缺失或类型不合法（内/外压力均须按约定显式给出）",
                "details": jsonable_encoder(exc.errors()),
            }
        },
    )


def _build_cylinder(
    params: StressRequest | ProfileRequest,
) -> LameCylinder:
    """统一入口：先校验，再构造 Lamé 解。校验不过直接抛 CalculationError。"""
    validate_geometry(params.inner_radius, params.outer_radius)
    validate_pressures(params.internal_pressure, params.external_pressure)
    validate_poisson_ratio(params.poisson_ratio)
    return LameCylinder.from_pressures(
        params.inner_radius,
        params.outer_radius,
        params.internal_pressure,
        params.external_pressure,
    )


def _to_point(state) -> StressPoint:
    return StressPoint(
        r=state.r,
        sigma_r=state.sigma_r,
        sigma_theta=state.sigma_theta,
        sigma_z=state.sigma_z,
        von_mises=state.von_mises,
    )


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/stress", response_model=StressResponse)
def compute_stress(request: StressRequest) -> StressResponse:
    """单根厚壁圆筒的应力计算：内外壁环向应力 + 指定半径处应力 + 屈服判定。"""
    cylinder = _build_cylinder(request)
    validate_query_radius(
        request.query_radius, request.inner_radius, request.outer_radius
    )
    if request.yield_strength is not None:
        validate_yield_strength(request.yield_strength)

    query_state = cylinder.stress_at(
        request.query_radius, request.end_condition, request.poisson_ratio
    )
    max_vm = cylinder.max_von_mises(
        request.end_condition, request.poisson_ratio
    )
    yielded = (
        max_vm > request.yield_strength
        if request.yield_strength is not None
        else None
    )

    return StressResponse(
        lame_A=cylinder.A,
        lame_B=cylinder.B,
        inner_hoop_stress=cylinder.hoop_stress(cylinder.inner_radius),
        outer_hoop_stress=cylinder.hoop_stress(cylinder.outer_radius),
        axial_stress=cylinder.axial_stress(
            request.end_condition, request.poisson_ratio
        ),
        query=_to_point(query_state),
        max_von_mises=max_vm,
        yield_strength=request.yield_strength,
        yielded=yielded,
    )


@app.post("/profile", response_model=ProfileResponse)
def compute_profile(request: ProfileRequest) -> ProfileResponse:
    """沿壁厚等间距采样，返回 σr(r)、σθ(r) 点列（逐点来自 Lamé 解）。"""
    cylinder = _build_cylinder(request)

    a, b = request.inner_radius, request.outer_radius
    n = request.n_points
    step = (b - a) / (n - 1)
    points = [
        _to_point(
            cylinder.stress_at(
                a + i * step, request.end_condition, request.poisson_ratio
            )
        )
        for i in range(n)
    ]

    return ProfileResponse(
        lame_A=cylinder.A,
        lame_B=cylinder.B,
        points=points,
        thin_wall_hoop_stress=thin_wall_hoop_stress(
            request.internal_pressure,
            request.external_pressure,
            cylinder.mean_radius,
            cylinder.thickness,
        ),
    )
