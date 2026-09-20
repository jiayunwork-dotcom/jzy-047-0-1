"""计算编排层：把校验、Lamé 解、薄壁对照组装成无状态的纯计算。

本模块不持有任何模块级可变状态：每次调用都新建 Cylinder 与 Lamé 系数，
所有中间量都是函数内的局部量，因此并发请求之间天然互不干扰。
"""

from __future__ import annotations

from typing import Any

from . import lame, thin_wall
from .domain import Cylinder, EndCondition
from .validation import (
    validate_cylinder,
    validate_profile_points,
    validate_radius,
    validate_yield_strength,
)


def build_cylinder(
    a: float,
    b: float,
    p_i: float,
    p_o: float,
    end_condition: str | EndCondition,
    poisson: float | None = None,
) -> Cylinder:
    """入口校验（在任何计算之前完成）。"""

    return validate_cylinder(a, b, p_i, p_o, end_condition, poisson)


def _point_dict(p: lame.StressPoint) -> dict[str, float]:
    return {
        "r": p.r,
        "sigma_r": p.sigma_r,
        "sigma_theta": p.sigma_theta,
        "sigma_z": p.sigma_z,
        "von_mises": lame.von_mises(p),
    }


def compute_stress(
    a: float,
    b: float,
    p_i: float,
    p_o: float,
    end_condition: str | EndCondition,
    poisson: float | None = None,
    r: float | None = None,
    yield_strength: float | None = None,
) -> dict[str, Any]:
    """单根筒的完整应力计算。

    返回内外壁应力（含两处环向应力）、Lamé 系数 A/B、
    可选的指定半径处应力、薄壁对照，以及可选的屈服判定。
    """

    cyl = build_cylinder(a, b, p_i, p_o, end_condition, poisson)
    coeffs = lame.solve_lame_coefficients(cyl)
    walls = lame.wall_stresses(cyl, coeffs)

    result: dict[str, Any] = {
        "inputs": {
            "a": cyl.a,
            "b": cyl.b,
            "p_i": cyl.p_i,
            "p_o": cyl.p_o,
            "end_condition": cyl.end_condition.value,
            "poisson": cyl.poisson,
        },
        "lame_coefficients": {"A": coeffs.A, "B": coeffs.B},
        "inner_wall": _point_dict(walls.inner),
        "outer_wall": _point_dict(walls.outer),
        # 单独再报一次内外壁环向应力，方便上游直接取用。
        "hoop_stress": {
            "inner": walls.inner.sigma_theta,
            "outer": walls.outer.sigma_theta,
        },
        "thin_wall_reference": thin_wall.thin_wall_comparison(cyl),
    }

    if r is not None:
        r_checked = validate_radius(r, cyl)
        result["at_radius"] = _point_dict(lame.stress_at(cyl, r_checked, coeffs))

    if yield_strength is not None:
        sigma_y = validate_yield_strength(yield_strength)
        # 内外壁两处取 von Mises 大者作为判定值（内压工况通常在内壁）。
        inner_vm = lame.von_mises(walls.inner)
        outer_vm = lame.von_mises(walls.outer)
        if inner_vm >= outer_vm:
            max_vm, max_location = inner_vm, "inner_wall"
        else:
            max_vm, max_location = outer_vm, "outer_wall"
        result["yield_check"] = {
            "yield_strength": sigma_y,
            "max_von_mises": max_vm,
            "max_location": max_location,
            "yielded": max_vm > sigma_y,
            "criterion": "von_mises",
        }

    return result


def compute_profile(
    a: float,
    b: float,
    p_i: float,
    p_o: float,
    end_condition: str | EndCondition,
    poisson: float | None = None,
    points: int = 21,
) -> dict[str, Any]:
    """沿壁厚等间距采样的 σr(r)、σθ(r)、σz(r) 点列。"""

    cyl = build_cylinder(a, b, p_i, p_o, end_condition, poisson)
    n = validate_profile_points(points)
    coeffs = lame.solve_lame_coefficients(cyl)
    pts = lame.stress_profile(cyl, n, coeffs)
    return {
        "inputs": {
            "a": cyl.a,
            "b": cyl.b,
            "p_i": cyl.p_i,
            "p_o": cyl.p_o,
            "end_condition": cyl.end_condition.value,
            "poisson": cyl.poisson,
        },
        "points": n,
        "radii": [p.r for p in pts],
        "sigma_r": [p.sigma_r for p in pts],
        "sigma_theta": [p.sigma_theta for p in pts],
        "sigma_z": [p.sigma_z for p in pts],
    }
