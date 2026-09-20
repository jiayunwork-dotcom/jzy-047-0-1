"""Lamé 边界条件求解：由内外壁压力确定系数 A、B。

符号约定
--------
- 拉应力为正、压应力为负；
- 压力 p_i、p_o 以"作用在壁面上的法向压力大小"为正输入，
  因此壁面处径向应力等于压力的负值：σr(a) = −p_i，σr(b) = −p_o。

推导（此处解释一次，全服务共用）
------------------------------
轴对称厚壁圆筒在无体力条件下的弹性解（Lamé, 1852）：

    σr(r) = A − B/r²        ……径向应力，B 项取 **负** 号
    σθ(r) = A + B/r²        ……环向应力，B 项取 **正** 号

注意 B 项在两个分量上的符号相反：径向为 −B/r²、环向为 +B/r²。
若两处写成同号，径向平衡方程 dσr/dr = (σθ − σr)/r 将不再成立，
整个应力场被翻转，结果完全错误。

将两个压力边界条件代入 σr 的表达式：

    A − B/a² = −p_i          ……(1) 内壁
    A − B/b² = −p_o          ……(2) 外壁

(2) − (1) 消去 A：

    B·(1/a² − 1/b²) = p_i − p_o
    B·(b² − a²)/(a²b²) = p_i − p_o

    ⇒  B = (p_i − p_o)·a²·b² / (b² − a²)

代回 (1)：

    A = B/a² − p_i = (p_i − p_o)·b²/(b² − a²) − p_i
    ⇒  A = (p_i·a² − p_o·b²) / (b² − a²)

两个系数完全由几何 (a, b) 与压力 (p_i, p_o) 决定，与材料常数无关
（材料只影响变形，不影响静定应力场）。
"""

from __future__ import annotations


def solve_lame_constants(
    inner_radius: float,
    outer_radius: float,
    internal_pressure: float,
    external_pressure: float,
) -> tuple[float, float]:
    """由内外壁压力边界条件解出 Lamé 系数 (A, B)。

    参数
    ----
    inner_radius : 内半径 a（> 0）
    outer_radius : 外半径 b（> a）
    internal_pressure : 内壁压力 p_i（压为正）
    external_pressure : 外壁压力 p_o（压为正）

    返回
    ----
    (A, B)：A 的量纲为应力，B 的量纲为应力·长度²。
    调用方须先完成输入合法性校验（见 validation 模块）。
    """
    a = inner_radius
    b = outer_radius
    p_i = internal_pressure
    p_o = external_pressure

    denominator = b * b - a * a  # b > a，保证为正
    A = (p_i * a * a - p_o * b * b) / denominator
    B = (p_i - p_o) * a * a * b * b / denominator
    return A, B
