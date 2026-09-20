"""Lamé 厚壁圆筒弹性应力场。

应力分量（拉为正）：

    σr(r) = A − B/r²     径向
    σθ(r) = A + B/r²     环向
    σz                  轴向，由端部约束方式决定

轴向约束三种情形分别处理：
- CLOSED（闭口）：端盖承压，筒身承受均匀轴向应力 σz = A
  （轴向力平衡：σz·π(b²−a²) = p_i·πa² − p_o·πb²）；
- OPEN（开口）：轴向自由，σz = 0；
- PLANE_STRAIN（平面应变）：εz = 0，由广义胡克定律
  εz = (σz − ν(σr + σθ))/E = 0 得 σz = ν(σr + σθ) = 2νA，
  泊松效应仅在此情形下进入应力结果。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from .boundary import solve_lame_constants


class EndCondition(str, Enum):
    """筒体两端的轴向约束方式。"""

    CLOSED = "closed"
    OPEN = "open"
    PLANE_STRAIN = "plane_strain"


@dataclass(frozen=True)
class StressState:
    """某一半径处的完整应力状态（量纲与输入压力一致）。"""

    r: float
    sigma_r: float
    sigma_theta: float
    sigma_z: float
    von_mises: float


@dataclass(frozen=True)
class LameCylinder:
    """一根厚壁圆筒的 Lamé 弹性解。

    不可变对象：构造后所有量均为纯函数计算，无任何共享可变状态，
    可被并发请求安全地各自创建、各自使用。
    """

    inner_radius: float
    outer_radius: float
    internal_pressure: float
    external_pressure: float
    A: float
    B: float

    @classmethod
    def from_pressures(
        cls,
        inner_radius: float,
        outer_radius: float,
        internal_pressure: float,
        external_pressure: float,
    ) -> "LameCylinder":
        A, B = solve_lame_constants(
            inner_radius, outer_radius, internal_pressure, external_pressure
        )
        return cls(
            inner_radius=inner_radius,
            outer_radius=outer_radius,
            internal_pressure=internal_pressure,
            external_pressure=external_pressure,
            A=A,
            B=B,
        )

    @property
    def thickness(self) -> float:
        return self.outer_radius - self.inner_radius

    @property
    def mean_radius(self) -> float:
        return 0.5 * (self.inner_radius + self.outer_radius)

    def radial_stress(self, r: float) -> float:
        """σr(r) = A − B/r²（B 项取负号）。"""
        return self.A - self.B / (r * r)

    def hoop_stress(self, r: float) -> float:
        """σθ(r) = A + B/r²（B 项取正号，与径向相反）。"""
        return self.A + self.B / (r * r)

    def axial_stress(
        self,
        end_condition: EndCondition,
        poisson_ratio: float = 0.3,
    ) -> float:
        """按轴向约束方式给出 σz（沿壁厚为常数）。"""
        if end_condition is EndCondition.CLOSED:
            return self.A
        if end_condition is EndCondition.OPEN:
            return 0.0
        if end_condition is EndCondition.PLANE_STRAIN:
            return 2.0 * poisson_ratio * self.A
        raise ValueError(f"未知的轴向约束方式: {end_condition!r}")

    def von_mises(
        self,
        r: float,
        end_condition: EndCondition,
        poisson_ratio: float = 0.3,
    ) -> float:
        """von Mises 等效应力（三向主应力即 σr、σθ、σz）。"""
        sr = self.radial_stress(r)
        st = self.hoop_stress(r)
        sz = self.axial_stress(end_condition, poisson_ratio)
        return math.sqrt(
            0.5 * ((sr - st) ** 2 + (st - sz) ** 2 + (sz - sr) ** 2)
        )

    def stress_at(
        self,
        r: float,
        end_condition: EndCondition,
        poisson_ratio: float = 0.3,
    ) -> StressState:
        """指定半径处的完整应力状态。"""
        return StressState(
            r=r,
            sigma_r=self.radial_stress(r),
            sigma_theta=self.hoop_stress(r),
            sigma_z=self.axial_stress(end_condition, poisson_ratio),
            von_mises=self.von_mises(r, end_condition, poisson_ratio),
        )

    def max_von_mises(
        self,
        end_condition: EndCondition,
        poisson_ratio: float = 0.3,
        samples: int = 201,
    ) -> float:
        """沿壁厚扫描求 von Mises 等效应力的最大值。

        对等间距采样点逐点计算 Lamé 解后取最大，用于屈服判定。
        """
        a, b = self.inner_radius, self.outer_radius
        step = (b - a) / (samples - 1)
        return max(
            self.von_mises(a + i * step, end_condition, poisson_ratio)
            for i in range(samples)
        )
