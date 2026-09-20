"""输入合法性校验：在计算入口处挡回一切不合法数据。

所有校验失败都抛出 CalculationError，携带结构化错误码与明确原因，
由 Web 层统一转换为带原因的 JSON 错误响应，绝不带病进入计算。
"""

from __future__ import annotations

import math


class CalculationError(Exception):
    """计算入口校验失败。code 为机器可读错误码，message 说明原因。"""

    def __init__(self, code: str, message: str) -> None:
        self.code = code
        self.message = message
        super().__init__(message)


def validate_geometry(inner_radius: float, outer_radius: float) -> None:
    """几何约束：内外半径必须为正的有限数，且外半径严格大于内半径。"""
    for name, value in (
        ("inner_radius", inner_radius),
        ("outer_radius", outer_radius),
    ):
        if not math.isfinite(value) or value <= 0.0:
            raise CalculationError(
                "INVALID_GEOMETRY",
                f"{name} 必须为正的有限数值，收到 {value!r}",
            )
    if outer_radius <= inner_radius:
        raise CalculationError(
            "INVALID_GEOMETRY",
            f"外半径必须严格大于内半径：outer_radius={outer_radius} "
            f"不大于 inner_radius={inner_radius}",
        )


def validate_pressures(
    internal_pressure: float, external_pressure: float
) -> None:
    """压力约束：两个压力都必须按约定给全且为有限数（压为正）。"""
    for name, value in (
        ("internal_pressure", internal_pressure),
        ("external_pressure", external_pressure),
    ):
        if not math.isfinite(value):
            raise CalculationError(
                "INVALID_PRESSURE",
                f"{name} 必须为有限数值（约定：压力以压为正），收到 {value!r}",
            )


def validate_poisson_ratio(poisson_ratio: float) -> None:
    """泊松比约束：各向同性弹性体的物理范围为 −1 < ν ≤ 0.5。

    平面应变条件下 σz = 2νA 直接依赖 ν，超出物理范围的取值
    （如 ν > 0.5 或 ν ≤ −1）必须在此挡回。
    """
    if not math.isfinite(poisson_ratio) or not (-1.0 < poisson_ratio <= 0.5):
        raise CalculationError(
            "INVALID_POISSON_RATIO",
            f"泊松比必须满足 −1 < ν ≤ 0.5（各向同性弹性体的物理范围），"
            f"收到 {poisson_ratio!r}",
        )


def validate_query_radius(
    query_radius: float, inner_radius: float, outer_radius: float
) -> None:
    """查询半径必须落在筒壁范围内 [a, b]。"""
    if not math.isfinite(query_radius) or not (
        inner_radius <= query_radius <= outer_radius
    ):
        raise CalculationError(
            "QUERY_RADIUS_OUT_OF_RANGE",
            f"查询半径必须位于壁厚区间 [{inner_radius}, {outer_radius}] 内，"
            f"收到 {query_radius!r}",
        )


def validate_yield_strength(yield_strength: float) -> None:
    """屈服强度必须为正的有限数。"""
    if not math.isfinite(yield_strength) or yield_strength <= 0.0:
        raise CalculationError(
            "INVALID_YIELD_STRENGTH",
            f"屈服强度必须为正的有限数值，收到 {yield_strength!r}",
        )
