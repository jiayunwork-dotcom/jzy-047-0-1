"""领域模型：贯穿各计算模块的纯数据结构。

约定（全服务统一，务必遵守）：
    * 压力 p_i / p_o 以“受压为正”输入（工程约定）；
    * 输出应力采用弹性力学约定：拉应力为正、压应力为负；
    * 各量只需量纲自洽即可（长度同量纲、应力/压力同量纲），服务不绑定单位制。
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class EndCondition(str, Enum):
    """筒体两端的轴向约束方式。

    CLOSED     —— 闭口筒：端盖承受压力，筒身存在恒定轴向应力；
    OPEN       —— 开口筒：轴向自由，轴向应力恒为零；
    PLANE_STRAIN —— 平面应变：轴向应变被约束为零，轴向应力由泊松效应产生。
    """

    CLOSED = "closed"
    OPEN = "open"
    PLANE_STRAIN = "plane_strain"


@dataclass(frozen=True)
class Cylinder:
    """一根厚壁圆筒的几何与载荷定义。"""

    a: float          # 内半径 r_i（>0）
    b: float          # 外半径 r_o（> a）
    p_i: float        # 内壁压力，受压为正（>=0）
    p_o: float        # 外壁压力，受压为正（>=0）
    end_condition: EndCondition
    poisson: float | None = None   # 泊松比，仅平面应变时必填


@dataclass(frozen=True)
class LameCoeffs:
    """Lamé 解 σr = A − B/r²、σθ = A + B/r² 中的两个常数。

    B 在径向分量取负号、在环向分量取正号，这两个符号相反，
    是 Lamé 解最容易写翻的地方（见 app/lame.py 的推导）。
    """

    A: float
    B: float
    a: float
    b: float
    p_i: float
    p_o: float


@dataclass(frozen=True)
class StressPoint:
    """单个半径处的三向应力（拉为正、压为负）。"""

    r: float
    sigma_r: float
    sigma_theta: float
    sigma_z: float


@dataclass(frozen=True)
class WallStresses:
    """内壁、外壁处的应力汇总。"""

    inner: StressPoint
    outer: StressPoint
