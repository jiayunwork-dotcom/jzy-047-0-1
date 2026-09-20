"""薄壁极限对照。

薄壁圆筒环向应力的经典近似：

    σθ ≈ (p_i − p_o) · r_m / t

其中 r_m 为平均半径、t 为壁厚、(p_i − p_o) 为净内压。
该式是 Lamé 解在 t/r → 0 时的极限，仅用于对照校验，
服务对外给出的应力一律来自 Lamé 精确解，不使用本近似。
"""

from __future__ import annotations


def thin_wall_hoop_stress(
    internal_pressure: float,
    external_pressure: float,
    mean_radius: float,
    thickness: float,
) -> float:
    """薄壁公式给出的环向应力估计值。"""
    net_pressure = internal_pressure - external_pressure
    return net_pressure * mean_radius / thickness
