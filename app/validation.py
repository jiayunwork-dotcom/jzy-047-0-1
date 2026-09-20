"""输入校验：在任何计算发生之前把不合法的输入挡回。

校验失败时抛出 :class:`ValidationError`，由 HTTP 层转成结构化错误响应，
绝不让坏数据流入 Lamé 计算得到“看似正常”的无意义结果。
"""

from __future__ import annotations

import math

from .domain import Cylinder, EndCondition

# 各向同性弹性材料泊松比的物理范围（开区间）：-1 < ν < 0.5；
# 常规工程材料还满足 0 <= ν < 0.5，这里取理论上界并略作收紧以避开奇异极限。
POISSON_MIN = -1.0
POISSON_MAX = 0.5

# 沿壁厚采样点数的合理上下限。
PROFILE_MIN_POINTS = 2
PROFILE_MAX_POINTS = 10_000


class ValidationError(ValueError):
    """带字段名与原因的结构化输入错误。"""

    def __init__(self, message: str, field: str | None = None) -> None:
        super().__init__(message)
        self.message = message
        self.field = field

    def to_detail(self) -> dict[str, str | None]:
        return {"field": self.field, "reason": self.message}


def _is_finite_number(x: object) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(float(x))


def _check_pressure(value: object, name: str) -> float:
    if not _is_finite_number(value):
        raise ValidationError(f"{name} 必须为有限数值", name)
    value = float(value)
    # 约定：压力以受压为正。统一要求非负，避免符号语义歧义。
    if value < 0.0:
        raise ValidationError(
            f"{name} 必须以受压为正、取非负值，收到 {value}", name
        )
    return value


def _check_poisson(value: object, *, required: bool) -> float | None:
    if value is None:
        if required:
            raise ValidationError(
                "平面应变状态必须提供泊松比 poisson", "poisson"
            )
        return None
    if not _is_finite_number(value):
        raise ValidationError("泊松比必须为有限数值", "poisson")
    value = float(value)
    if not (POISSON_MIN < value < POISSON_MAX):
        raise ValidationError(
            f"泊松比必须落在物理区间 ({POISSON_MIN}, {POISSON_MAX}) 内，收到 {value}",
            "poisson",
        )
    return value


def validate_cylinder(
    a: object,
    b: object,
    p_i: object,
    p_o: object,
    end_condition: object,
    poisson: object = None,
) -> Cylinder:
    """把原始入参构造成经过完整校验的 :class:`Cylinder`。"""

    if not _is_finite_number(a):
        raise ValidationError("内半径 a 必须为正的有限数值", "a")
    if not _is_finite_number(b):
        raise ValidationError("外半径 b 必须为正的有限数值", "b")
    a, b = float(a), float(b)
    if a <= 0.0:
        raise ValidationError(f"内半径 a 必须为正，收到 {a}", "a")
    if b <= 0.0:
        raise ValidationError(f"外半径 b 必须为正，收到 {b}", "b")
    if b <= a:
        raise ValidationError(
            f"外半径 b 必须严格大于内半径 a（收到 a={a}, b={b}）", "b"
        )

    p_i = _check_pressure(p_i, "p_i")
    p_o = _check_pressure(p_o, "p_o")

    try:
        end = EndCondition(end_condition)
    except (ValueError, TypeError):
        allowed = ", ".join(c.value for c in EndCondition)
        raise ValidationError(
            f"轴向约束方式必须是 {allowed} 之一，收到 {end_condition!r}",
            "end_condition",
        ) from None

    poisson = _check_poisson(
        poisson, required=(end is EndCondition.PLANE_STRAIN)
    )

    return Cylinder(
        a=a, b=b, p_i=p_i, p_o=p_o, end_condition=end, poisson=poisson
    )


def validate_radius(r: object, cyl: Cylinder) -> float:
    """待查询的半径必须落在筒壁 [a, b] 内。"""

    if not _is_finite_number(r):
        raise ValidationError("查询半径 r 必须为有限数值", "r")
    r = float(r)
    if not (cyl.a <= r <= cyl.b):
        raise ValidationError(
            f"查询半径 r={r} 超出筒壁范围 [{cyl.a}, {cyl.b}]", "r"
        )
    return r


def validate_profile_points(n: object) -> int:
    if not _is_finite_number(n) or float(n) != int(float(n)):
        raise ValidationError("采样点数必须为整数", "points")
    n = int(n)
    if not (PROFILE_MIN_POINTS <= n <= PROFILE_MAX_POINTS):
        raise ValidationError(
            f"采样点数必须在 [{PROFILE_MIN_POINTS}, {PROFILE_MAX_POINTS}] 内，收到 {n}",
            "points",
        )
    return n


def validate_yield_strength(sigma_y: object) -> float:
    if not _is_finite_number(sigma_y):
        raise ValidationError("屈服强度必须为正的有限数值", "yield_strength")
    sigma_y = float(sigma_y)
    if sigma_y <= 0.0:
        raise ValidationError(f"屈服强度必须为正，收到 {sigma_y}", "yield_strength")
    return sigma_y
