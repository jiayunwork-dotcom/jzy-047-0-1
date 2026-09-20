"""HTTP 接口的请求 / 响应模型（Pydantic）。

internal_pressure 与 external_pressure 均不设默认值：
压力符号约定（压为正）要求调用方把两个压力显式给全，
缺漏会在 schema 层以 422 结构化错误挡回。
"""

from __future__ import annotations

from pydantic import BaseModel, Field

from .lame import EndCondition


class CylinderParameters(BaseModel):
    """一根厚壁圆筒的完整参数。"""

    inner_radius: float = Field(description="内半径 a，> 0")
    outer_radius: float = Field(description="外半径 b，> a")
    internal_pressure: float = Field(description="内壁压力 p_i，压为正")
    external_pressure: float = Field(description="外壁压力 p_o，压为正")
    end_condition: EndCondition = Field(
        default=EndCondition.CLOSED,
        description="轴向约束：closed / open / plane_strain",
    )
    poisson_ratio: float = Field(
        default=0.3, description="泊松比 ν，仅平面应变条件参与应力计算"
    )


class StressRequest(CylinderParameters):
    """单点应力计算请求。"""

    query_radius: float = Field(description="查询半径 r，须位于 [a, b]")
    yield_strength: float | None = Field(
        default=None,
        description="屈服强度（可选）。给出时按 von Mises 等效应力判断是否屈服",
    )


class StressPoint(BaseModel):
    """某一半径处的应力分量。"""

    r: float
    sigma_r: float
    sigma_theta: float
    sigma_z: float
    von_mises: float


class StressResponse(BaseModel):
    """单点应力计算响应。"""

    lame_A: float = Field(description="Lamé 系数 A（量纲：应力）")
    lame_B: float = Field(description="Lamé 系数 B（量纲：应力·长度²）")
    inner_hoop_stress: float = Field(description="内壁环向应力 σθ(a)")
    outer_hoop_stress: float = Field(description="外壁环向应力 σθ(b)")
    axial_stress: float = Field(description="轴向应力 σz（沿壁厚为常数）")
    query: StressPoint = Field(description="指定半径处的应力状态")
    max_von_mises: float = Field(description="全壁厚 von Mises 等效应力最大值")
    yield_strength: float | None = Field(description="回显的屈服强度，未给为 null")
    yielded: bool | None = Field(
        description="是否屈服（max_von_mises > yield_strength），未给屈服强度为 null"
    )


class ProfileRequest(CylinderParameters):
    """沿壁厚采样的应力分布请求。"""

    n_points: int = Field(
        default=50, ge=2, le=5000, description="等间距采样点数"
    )


class ProfileResponse(BaseModel):
    """沿壁厚的应力分布响应。"""

    lame_A: float
    lame_B: float
    points: list[StressPoint] = Field(
        description="等间距采样点列，逐点由 Lamé 解计算"
    )
    thin_wall_hoop_stress: float = Field(
        description="薄壁公式 (p_i−p_o)·r_m/t 给出的对照值，仅供比较"
    )
