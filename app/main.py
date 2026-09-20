"""FastAPI 入口：厚壁圆筒 Lamé 弹性应力计算服务。

路由：
    GET  /health                       健康检查
    POST /api/v1/stress                单根筒完整计算（含屈服判定）
    POST /api/v1/stress/profile        沿壁厚等间距采样点列
    GET  /api/v1/examples/internal-pressure  仅内压示范算例

服务完全无状态：每次请求独立计算，可安全并发。
"""

from __future__ import annotations

from typing import Any

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from . import service
from .schemas import ProfileRequest, StressRequest
from .validation import ValidationError

app = FastAPI(
    title="厚壁圆筒 Lamé 应力计算服务",
    version="1.0.0",
    description="给定几何、内外压与轴向约束，计算厚壁圆筒弹性应力场（拉为正、压为负）。",
)


@app.exception_handler(ValidationError)
async def _domain_validation_handler(request: Request, exc: ValidationError) -> JSONResponse:
    return JSONResponse(status_code=422, content=_error_body([exc.to_detail()]))


@app.exception_handler(RequestValidationError)
async def _request_validation_handler(
    request: Request, exc: RequestValidationError
) -> JSONResponse:
    details = [
        {
            "field": ".".join(str(x) for x in err["loc"] if x not in ("body",)),
            "reason": err["msg"],
        }
        for err in exc.errors()
    ]
    return JSONResponse(status_code=422, content=_error_body(details))


def _error_body(details: list[dict[str, str | None]]) -> dict[str, Any]:
    return {"error": "validation_failed", "details": details}


@app.get("/health")
async def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/api/v1/stress")
async def calculate_stress(req: StressRequest) -> dict[str, Any]:
    # 纯 CPU 计算且无共享状态；service 内部还会做一次领域级校验。
    return service.compute_stress(
        a=req.a,
        b=req.b,
        p_i=req.p_i,
        p_o=req.p_o,
        end_condition=req.end_condition,
        poisson=req.poisson,
        r=req.r,
        yield_strength=req.yield_strength,
    )


@app.post("/api/v1/stress/profile")
async def calculate_profile(req: ProfileRequest) -> dict[str, Any]:
    return service.compute_profile(
        a=req.a,
        b=req.b,
        p_i=req.p_i,
        p_o=req.p_o,
        end_condition=req.end_condition,
        poisson=req.poisson,
        points=req.points,
    )


# ─────────────────────────────────────────────────────────────────────────
# 示范算例（仅内压）：a=100, b=200, p_i=100, p_o=0，闭口
#   内壁环向 σθ(a) = p_i (b²+a²)/(b²−a²) = 100·50000/30000 = 166.67（拉，最大）
#   外壁环向 σθ(b) = p_i·2a²/(b²−a²)     = 100·20000/30000 =  66.67（拉）
#   内壁径向 σr(a) = −100；外壁径向 σr(b) = 0
#   闭口轴向 σz    = p_i a²/(b²−a²)       = 100·10000/30000 =  33.33（拉）
# 所有数字均可用 Lamé 公式手工核对，且内壁环向应力为全场最大拉应力。
# ─────────────────────────────────────────────────────────────────────────
EXAMPLE_PARAMS = {
    "a": 100.0,
    "b": 200.0,
    "p_i": 100.0,
    "p_o": 0.0,
    "end_condition": "closed",
    "poisson": None,
}


@app.get("/api/v1/examples/internal-pressure")
async def internal_pressure_example() -> dict[str, Any]:
    result = service.compute_stress(
        r=EXAMPLE_PARAMS["a"],  # 同时示范指定半径查询（取内壁）
        yield_strength=250.0,
        **EXAMPLE_PARAMS,
    )
    result["example_notes"] = {
        "description": "仅受内压的闭口厚壁筒（a=100, b=200, p_i=100, p_o=0）",
        "hand_check": {
            "inner_hoop": "p_i(b²+a²)/(b²−a²) = 100·50000/30000 = 166.67（全场最大拉应力）",
            "outer_hoop": "2·p_i·a²/(b²−a²) = 66.67",
            "inner_radial": "−p_i = −100（受压）",
            "outer_radial": "0（无外压，自由外表面）",
            "axial_closed": "p_i·a²/(b²−a²) = 33.33",
        },
        "max_tensile_location": "inner_wall hoop",
    }
    return result
