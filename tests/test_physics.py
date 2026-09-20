"""物理判据测试（逐条自动化守）：

1. 径向平衡方程 dσr/dr = (σθ−σr)/r 在筒壁内处处成立（数值差分抽查）；
2. 薄壁极限：b/a=1.1 时 Lamé 与薄膜公式很接近；b/a=2 时
   内壁环向应力必须显著高于薄壁估计（禁止全程套薄壁近似冒充精确解）；
3. 仅内压时外壁径向应力恰好为零。
"""

from __future__ import annotations

import math

import pytest

from app.lame import solve_lame_coefficients, stress_at, wall_stresses
from app.thin_wall import thin_wall_comparison, thin_wall_hoop, wall_ratio
from app.validation import validate_cylinder


# ── 判据 1：径向平衡方程 ────────────────────────────────────────────────
def equilibrium_residual(cyl, r, h_rel=1e-5):
    """中心差分计算 dσr/dr 与 (σθ−σr)/r 之差。"""

    c = solve_lame_coefficients(cyl)
    h = max(r * h_rel, 1e-12)
    r_lo, r_hi = max(cyl.a, r - h), min(cyl.b, r + h)
    h = r_hi - r_lo  # 边界附近自适应缩短后重新取对称步长
    r_mid = 0.5 * (r_lo + r_hi)
    sr_lo = stress_at(cyl, r_lo, c).sigma_r
    sr_hi = stress_at(cyl, r_hi, c).sigma_r
    d_sr_dr = (sr_hi - sr_lo) / h
    p_mid = stress_at(cyl, r_mid, c)
    rhs = (p_mid.sigma_theta - p_mid.sigma_r) / r_mid
    return d_sr_dr - rhs


@pytest.mark.parametrize(
    "p_i,p_o",
    [(100.0, 0.0), (0.0, 100.0), (120.0, 40.0), (75.0, 75.0)],
)
def test_radial_equilibrium_at_sample_radii(p_i, p_o):
    cyl = validate_cylinder(100.0, 200.0, p_i, p_o, "closed")
    # 在内壁、外壁附近及壁厚中部多点抽查（含三种端部条件的代表）。
    radii = [
        cyl.a * (1 + 1e-4),
        cyl.a + 0.25 * (cyl.b - cyl.a),
        0.5 * (cyl.a + cyl.b),
        cyl.a + 0.75 * (cyl.b - cyl.a),
        cyl.b * (1 - 1e-4),
    ]
    scale = max(p_i, p_o, 1.0) / cyl.a  # 残差的特征量级
    for r in radii:
        assert abs(equilibrium_residual(cyl, r)) < 1e-6 * scale


def test_radial_equilibrium_plane_strain():
    # σr/σθ 与端部条件无关，但平面应变路径也单独走一遍。
    cyl = validate_cylinder(50.0, 75.0, 90.0, 20.0, "plane_strain", 0.25)
    for r in [55.0, 62.5, 70.0]:
        assert abs(equilibrium_residual(cyl, r)) < 1e-8


# ── 判据 2：薄壁极限 ────────────────────────────────────────────────────
def test_thin_wall_close_at_ratio_1_1():
    # b/a = 1.1：Lamé 内壁环向应力与 p·r_m/t 的相对偏差应在 ~5% 以内。
    cyl = validate_cylinder(100.0, 110.0, 100.0, 0.0, "open")
    comp = thin_wall_comparison(cyl)
    assert wall_ratio(cyl) == pytest.approx(1.1)
    assert abs(comp["relative_error"]) < 0.05
    # 解析上：Lamé/薄壁 = 2(K²+1)/(K+1)²，K=1.1 时偏差 = 4.42/4.41−1 ≈ 0.00227。
    K = 1.1
    expected = 2.0 * (K * K + 1.0) / (K + 1.0) ** 2 - 1.0
    assert comp["relative_error"] == pytest.approx(expected, rel=1e-3)


def test_thin_wall_significantly_different_at_ratio_2():
    # b/a = 2：内壁环向必须显著高于薄壁估计（此处 166.7 vs 150，差 11.1%）。
    cyl = validate_cylinder(100.0, 200.0, 100.0, 0.0, "open")
    comp = thin_wall_comparison(cyl)
    assert comp["lame_inner_hoop"] > 1.10 * comp["thin_wall_hoop"]
    assert comp["lame_inner_hoop"] == pytest.approx(166.6666667, rel=1e-6)
    assert comp["thin_wall_hoop"] == pytest.approx(150.0)


def test_thin_wall_convergence_trend():
    # b/a 从 2 趋近 1，相对偏差应单调趋于 0。
    errs = []
    for ratio in [2.0, 1.5, 1.2, 1.05, 1.01]:
        cyl = validate_cylinder(100.0, 100.0 * ratio, 100.0, 0.0, "open")
        errs.append(abs(thin_wall_comparison(cyl)["relative_error"]))
    assert errs == sorted(errs, reverse=True)
    assert errs[-1] < 0.005


def test_service_reports_lame_not_thin_wall_at_thick_case():
    # 防伪装：厚壁工况下服务给出的内壁环向必须是 Lamé 值而非薄膜值。
    cyl = validate_cylinder(100.0, 200.0, 100.0, 0.0, "closed")
    w = wall_stresses(cyl)
    assert w.inner.sigma_theta != pytest.approx(thin_wall_hoop(cyl))
    assert w.inner.sigma_theta == pytest.approx(5.0 / 3.0 * 100.0)


# ── 判据 3：仅内压时外壁径向应力为零 ───────────────────────────────────
@pytest.mark.parametrize("ratio", [1.01, 1.1, 1.5, 2.0, 5.0])
def test_outer_radial_stress_exactly_zero_internal_only(ratio):
    cyl = validate_cylinder(100.0, 100.0 * ratio, 100.0, 0.0, "closed")
    assert wall_stresses(cyl).outer.sigma_r == pytest.approx(0.0, abs=1e-12)


def test_outer_radial_equals_negative_external_pressure():
    # 推广：有外压时该值应恰好等于 -p_o。
    cyl = validate_cylinder(100.0, 200.0, 100.0, 35.0, "open")
    assert wall_stresses(cyl).outer.sigma_r == pytest.approx(-35.0, abs=1e-12)
