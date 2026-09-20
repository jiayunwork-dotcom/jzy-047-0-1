"""Lamé 厚壁圆筒弹性解：边界条件求解与三向应力场。

────────────────────────────────────────────────────────────────────────
一、控制方程与 Lamé 形式
────────────────────────────────────────────────────────────────────────
轴对称、轴向不变的弹性问题，径向平衡方程（无体积力）为

    dσr/dr = (σθ − σr) / r.

轴对称相容方程配合广义胡克定律可导出 σr + σθ 沿壁厚为常数，
故应力场的通解写成 Lamé 形式：

    σr(r) = A − B / r²,
    σθ(r) = A + B / r².          （B 项符号与径向相反！）

特别注意：B/r² 在 σr 中取**负号**、在 σθ 中取**正号**，
两边符号相反。若两处同号，整个应力场的拉压分布会翻过来。

────────────────────────────────────────────────────────────────────────
二、由压力边界条件求解 A、B
────────────────────────────────────────────────────────────────────────
边界上的面力与径向应力连续。注意输入压力按“受压为正”的工程约定，
而应力按“拉为正、压为负”输出，因此：

    r = a（内壁）：σr(a) = −p_i = A − B/a²,
    r = b（外壁）：σr(b) = −p_o = A − B/b².

两式相减消去 A：
    −p_i + p_o = −B/a² + B/b² = −B(1/a² − 1/b²)
    ⇒ B = (p_i − p_o) / (1/a² − 1/b²)
         = a² b² (p_i − p_o) / (b² − a²).

回代求 A：
    A = −p_i + B/a²
      = (a² p_i − b² p_o) / (b² − a²).

可直接验证：A − B/b² = −p_o，即外壁边界条件同样满足。

符号自检：仅受内压时 p_i > 0、p_o = 0，
    B > 0，于是 σθ(a) = A + B/a² = p_i (b²+a²)/(b²−a²) > 0（内壁环向受拉，
且为全场最大拉应力），σr(a) = −p_i < 0（径向受压）——物理上正确。

────────────────────────────────────────────────────────────────────────
三、轴向应力（取决于端部约束）
────────────────────────────────────────────────────────────────────────
σr、σθ 与端部条件无关；σz 在截面上恒为常数：

* 闭口筒 CLOSED：端盖承压，由端部合力平衡
      σz · π(b²−a²) = p_i π a² − p_o π b²
      ⇒ σz = (p_i a² − p_o b²) / (b² − a²) = A.
  （闭口条件下 σz 数值恰等于 A，这是巧合也是可验证的恒等式。）

* 开口筒 OPEN：端部自由无轴力，σz ≡ 0。

* 平面应变 PLANE_STRAIN：εz = 0，由胡克定律
      εz = (1/E)[σz − ν(σr + σθ)] = 0
      ⇒ σz = ν(σr + σθ) = 2ν A   （沿壁厚为常数）。
"""

from __future__ import annotations

from .domain import Cylinder, EndCondition, LameCoeffs, StressPoint, WallStresses


def solve_lame_coefficients(cyl: Cylinder) -> LameCoeffs:
    """按内外壁压力边界条件求 Lamé 常数 A、B。

    A = (a²·p_i − b²·p_o) / (b² − a²)
    B = a²·b²·(p_i − p_o) / (b² − a²)
    """

    a, b = cyl.a, cyl.b
    a2, b2 = a * a, b * b
    denom = b2 - a2  # 由校验保证 b > a > 0，分母严格为正

    A = (a2 * cyl.p_i - b2 * cyl.p_o) / denom
    B = a2 * b2 * (cyl.p_i - cyl.p_o) / denom
    return LameCoeffs(A=A, B=B, a=a, b=b, p_i=cyl.p_i, p_o=cyl.p_o)


def _axial_stress(cyl: Cylinder, A: float) -> float:
    """按端部约束计算沿壁厚恒定的轴向应力（拉为正）。"""

    if cyl.end_condition is EndCondition.CLOSED:
        # 端盖合力平衡：σz = (p_i a² − p_o b²)/(b² − a²) = A
        return A
    if cyl.end_condition is EndCondition.OPEN:
        return 0.0
    # 平面应变：εz = 0 ⇒ σz = ν(σr + σθ) = 2νA
    assert cyl.poisson is not None  # 校验阶段已保证
    return 2.0 * cyl.poisson * A


def stress_at(cyl: Cylinder, r: float, coeffs: LameCoeffs | None = None) -> StressPoint:
    """计算半径 r 处的 (σr, σθ, σz)。

    σr = A − B/r²；σθ = A + B/r²。B 项两处符号相反，切勿写同号。
    """

    if coeffs is None:
        coeffs = solve_lame_coefficients(cyl)
    inv_r2 = 1.0 / (r * r)
    sigma_r = coeffs.A - coeffs.B * inv_r2
    sigma_theta = coeffs.A + coeffs.B * inv_r2
    sigma_z = _axial_stress(cyl, coeffs.A)
    return StressPoint(r=r, sigma_r=sigma_r, sigma_theta=sigma_theta, sigma_z=sigma_z)


def wall_stresses(cyl: Cylinder, coeffs: LameCoeffs | None = None) -> WallStresses:
    """内壁 r=a 与外壁 r=b 两处的应力。"""

    if coeffs is None:
        coeffs = solve_lame_coefficients(cyl)
    return WallStresses(
        inner=stress_at(cyl, cyl.a, coeffs),
        outer=stress_at(cyl, cyl.b, coeffs),
    )


def stress_profile(
    cyl: Cylinder, points: int, coeffs: LameCoeffs | None = None
) -> list[StressPoint]:
    """沿壁厚在 [a, b] 上取 points 个等间距半径，逐点计算 Lamé 应力。

    这些点全部由 stress_at 逐点算出（即真来自 Lamé 解），
    不使用任何预设的衰减形状。
    """

    if points < 2:
        raise ValueError("采样点数至少为 2")
    if coeffs is None:
        coeffs = solve_lame_coefficients(cyl)
    step = (cyl.b - cyl.a) / (points - 1)
    return [
        stress_at(cyl, cyl.a + i * step, coeffs)
        for i in range(points)
    ]


def von_mises(point: StressPoint) -> float:
    """三主应力 (σr, σθ, σz) 下的 von Mises 等效应力。

    σ_vm = sqrt(0.5·[(σr−σθ)² + (σθ−σz)² + (σz−σr)²])
    """

    sr, st, sz = point.sigma_r, point.sigma_theta, point.sigma_z
    return (0.5 * ((sr - st) ** 2 + (st - sz) ** 2 + (sz - sr) ** 2)) ** 0.5


def deviatoric_diameter(point: StressPoint) -> float:
    """环向与径向之差 σθ − σr = 2B/r²，即偏应力的一个直接度量。"""

    return point.sigma_theta - point.sigma_r
