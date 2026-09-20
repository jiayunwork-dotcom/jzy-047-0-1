"""薄壁极限对照：用薄壁薄膜公式给 Lamé 精确解提供量级参照。

薄壁圆管理论把环向应力沿壁厚视为均匀：

    σθ ≈ p · r_m / t,   r_m = (a + b)/2（平均半径），t = b − a（壁厚）。

内压与外压同时作用时，净压差 Δp = p_i − p_o 是产生环向薄膜力的来源，
故对照式写为

    σθ_thin = (p_i − p_o) · r_m / t.

当壁厚比 K = b/a → 1 时，Lamé 内壁环向应力
    σθ(a) = [p_i(a²+b²) − 2 p_o b²] / (b²−a²)
应收敛到该薄壁值（对净压差而言，相对偏差约为 t/(2r_m) 量级）。
厚壁比 K 明显大于 1 时，内壁环向应力显著高于薄壁估计——
该对照模块只负责如实给出两者，供调用方/测试判定，绝不用近似值冒充精确解。
"""

from __future__ import annotations

from .domain import Cylinder


def wall_ratio(cyl: Cylinder) -> float:
    """厚径比 K = b/a（>1；越接近 1 越趋近薄壁）。"""

    return cyl.b / cyl.a


def thin_wall_hoop(cyl: Cylinder) -> float:
    """薄壁薄膜公式下的环向应力（净压差 × 平均半径 / 壁厚）。"""

    r_m = 0.5 * (cyl.a + cyl.b)
    t = cyl.b - cyl.a
    return (cyl.p_i - cyl.p_o) * r_m / t


def thin_wall_comparison(cyl: Cylinder) -> dict[str, float]:
    """返回薄壁对照量：厚径比、薄壁环向应力、与 Lamé 内壁值的相对偏差。"""

    from .lame import wall_stresses  # 局部导入避免循环依赖

    sigma_thin = thin_wall_hoop(cyl)
    sigma_inner_lame = wall_stresses(cyl).inner.sigma_theta
    if sigma_thin != 0.0:
        relative_error = (sigma_inner_lame - sigma_thin) / abs(sigma_thin)
    else:
        # 均匀静水压等情形下薄壁值为 0，无“相对”偏差可言，报绝对差。
        relative_error = float("inf") if sigma_inner_lame != 0.0 else 0.0
    return {
        "wall_ratio": wall_ratio(cyl),
        "thin_wall_hoop": sigma_thin,
        "lame_inner_hoop": sigma_inner_lame,
        "relative_error": relative_error,
    }
