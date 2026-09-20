"""Lamé 解的直接正确性测试：边界条件、系数、符号、手工算例。"""

from __future__ import annotations

import math

import pytest

from app.domain import EndCondition
from app.lame import (
    deviatoric_diameter,
    solve_lame_coefficients,
    stress_at,
    stress_profile,
    von_mises,
    wall_stresses,
)
from app.validation import validate_cylinder


def make(a=100.0, b=200.0, p_i=100.0, p_o=0.0, end="closed", poisson=None):
    return validate_cylinder(a, b, p_i, p_o, end, poisson)


class TestLameCoefficients:
    def test_internal_pressure_hand_calc(self):
        # a=100,b=200,p_i=100,p_o=0：
        # A = a²p_i/(b²-a²) = 1e4·100/3e4 = 33.333...
        # B = a²b²p_i/(b²-a²) = 1e4·4e4·100/3e4 = 4e6/3 ≈ 1.3333e6
        cyl = make()
        c = solve_lame_coefficients(cyl)
        assert c.A == pytest.approx(100.0 / 3.0)
        assert c.B == pytest.approx(4.0e6 / 3.0)

    def test_boundary_conditions_satisfied(self):
        # 任意内外压组合下，σr(a)=-p_i、σr(b)=-p_o 必须精确成立。
        for p_i, p_o in [(100.0, 0.0), (0.0, 80.0), (120.0, 45.0)]:
            cyl = make(p_i=p_i, p_o=p_o)
            w = wall_stresses(cyl)
            assert w.inner.sigma_r == pytest.approx(-p_i)
            assert w.outer.sigma_r == pytest.approx(-p_o)

    def test_B_sign_opposite_between_radial_and_hoop(self):
        # 关键防回归：B 项在 σr 取负、σθ 取正。
        # 仅内压时 B>0，故 σθ−σr = 2B/r² > 0，内壁处相差 2p_i·b²/(b²−a²)。
        cyl = make(p_i=100.0, p_o=0.0)
        c = solve_lame_coefficients(cyl)
        assert c.B > 0
        p = stress_at(cyl, 125.0, c)
        assert p.sigma_theta - p.sigma_r == pytest.approx(2.0 * c.B / 125.0**2)
        # 内壁 σθ>0（拉）、σr<0（压）：同号错误会让这条翻转。
        assert p.sigma_theta > 0 and p.sigma_r < 0

    def test_external_pressure_compressive_hoop(self):
        # 仅外压时环向应力处处为压；内壁处压得最重：
        # σθ(a) = -2 p_o b²/(b²-a²) = -266.67，σθ(b) = -p_o(b²+a²)/(b²-a²) = -166.67
        cyl = make(p_i=0.0, p_o=100.0)
        w = wall_stresses(cyl)
        assert w.inner.sigma_theta < 0
        assert w.outer.sigma_theta < 0
        assert w.outer.sigma_theta == pytest.approx(-500.0 / 3.0)
        assert w.inner.sigma_theta == pytest.approx(-800.0 / 3.0)


class TestExampleCase:
    """示范算例：内壁环向应力 166.67 为全场最大拉应力，可手工核对。"""

    def test_wall_values(self):
        cyl = make()
        w = wall_stresses(cyl)
        assert w.inner.sigma_theta == pytest.approx(166.6666667, rel=1e-6)
        assert w.outer.sigma_theta == pytest.approx(66.6666667, rel=1e-6)
        assert w.inner.sigma_r == pytest.approx(-100.0)
        assert w.outer.sigma_r == pytest.approx(0.0, abs=1e-12)

    def test_inner_hoop_is_max_tensile(self):
        cyl = make()
        pts = stress_profile(cyl, 200)
        max_theta = max(p.sigma_theta for p in pts)
        assert pts[0].sigma_theta == pytest.approx(max_theta)
        # 环向应力沿壁厚单调下降。
        thetas = [p.sigma_theta for p in pts]
        assert all(thetas[i] > thetas[i + 1] for i in range(len(thetas) - 1))

    def test_radial_monotonic_compression(self):
        cyl = make()
        pts = stress_profile(cyl, 200)
        rs = [p.sigma_r for p in pts]
        # 径向应力从 -p_i 单调升到 0，处处为压或零。
        assert all(rs[i] < rs[i + 1] for i in range(len(rs) - 1))
        assert all(v <= 1e-12 for v in rs)


class TestAxialConditions:
    def test_closed_axial_equals_force_balance(self):
        cyl = make(p_i=100.0, p_o=30.0, end="closed")
        c = solve_lame_coefficients(cyl)
        p = stress_at(cyl, 150.0, c)
        expected = (100.0 * 100.0**2 - 30.0 * 200.0**2) / (200.0**2 - 100.0**2)
        assert p.sigma_z == pytest.approx(expected)
        assert p.sigma_z == pytest.approx(c.A)  # 闭口恒等式 σz = A

    def test_open_axial_zero_everywhere(self):
        cyl = make(end="open")
        pts = stress_profile(cyl, 10)
        assert all(p.sigma_z == 0.0 for p in pts)

    def test_plane_strain_axial_is_2nuA(self):
        nu = 0.3
        cyl = make(end="plane_strain", poisson=nu)
        c = solve_lame_coefficients(cyl)
        for p in stress_profile(cyl, 10, c):
            assert p.sigma_z == pytest.approx(2.0 * nu * c.A)
            # 平面应变下轴向应变应为零（由构造保证，数值复核一次）。
            # εz 无法在此直接算 E，但 σz = ν(σr+σθ) 即等价条件：
            assert p.sigma_z == pytest.approx(nu * (p.sigma_r + p.sigma_theta))

    def test_radial_hoop_independent_of_end_condition(self):
        ref = wall_stresses(make(end="closed"))
        for end in ("open", "plane_strain"):
            w = wall_stresses(make(end=end, poisson=0.3))
            assert w.inner.sigma_r == pytest.approx(ref.inner.sigma_r)
            assert w.inner.sigma_theta == pytest.approx(ref.inner.sigma_theta)
            assert w.outer.sigma_theta == pytest.approx(ref.outer.sigma_theta)


class TestVonMises:
    def test_hydrostatic_zero_vm(self):
        # σr=σθ=σz（静水压球应力）时等效应力为零。
        cyl = make(p_i=50.0, p_o=50.0, end="closed")
        p = stress_at(cyl, 150.0)
        # 均匀内外压：A=-p, B=0 ⇒ σr=σθ=-p；闭口 σz=A=-p。
        assert von_mises(p) == pytest.approx(0.0, abs=1e-12)
        assert deviatoric_diameter(p) == pytest.approx(0.0, abs=1e-12)

    def test_vm_positive_under_internal_pressure(self):
        cyl = make()
        assert von_mises(stress_at(cyl, 100.0)) > 0


class TestProfile:
    def test_profile_endpoints_and_spacing(self):
        cyl = make()
        pts = stress_profile(cyl, 11)
        assert len(pts) == 11
        assert pts[0].r == pytest.approx(100.0)
        assert pts[-1].r == pytest.approx(200.0)
        drs = [pts[i + 1].r - pts[i].r for i in range(10)]
        assert all(d == pytest.approx(10.0) for d in drs)

    def test_profile_recomputes_per_geometry(self):
        # 点列必须逐点来自 Lamé 解：换几何后数值随之改变，而非固定形状。
        p1 = stress_profile(make(b=150.0), 5)
        p2 = stress_profile(make(b=300.0), 5)
        assert [p.sigma_theta for p in p1] != pytest.approx(
            [p.sigma_theta for p in p2]
        )
