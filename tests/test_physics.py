"""物理正确性判据与联动关系测试。

每条判据相互独立，逐条把服务的物理行为卡死：
1. 径向平衡方程 dσr/dr = (σθ − σr)/r 在壁内处处成立（数值差分抽查）；
2. 薄壁极限：厚比 K=1.1 时环向应力贴近薄壁公式，K=2 时内壁环向
   应力显著高于薄壁估计（防止全程偷换薄壁近似）；
3. 仅内压作用时外壁径向应力恰为零；
4. 线性：内压单独放大两倍，整个应力场按比例放大两倍；
5. 静水压：叠加与内压相等的外压后趋于均匀静水状态，偏应力下降；
6. 相似性：几何按同一系数缩放后，以 r/a 为横坐标的应力分布不变；
7. 三种轴向约束分别给出正确的 σz。
"""

import math

import pytest

from app.boundary import solve_lame_constants
from app.lame import EndCondition, LameCylinder
from app.thinwall import thin_wall_hoop_stress

# 测试里独立手写 Lamé 公式作为对照，避免与被测实现同源漂移。
def reference_stresses(a, b, p_i, p_o, r):
    A = (p_i * a**2 - p_o * b**2) / (b**2 - a**2)
    B = (p_i - p_o) * a**2 * b**2 / (b**2 - a**2)
    return A - B / r**2, A + B / r**2


def make_cylinder(a=100.0, b=200.0, p_i=10.0, p_o=0.0):
    return LameCylinder.from_pressures(a, b, p_i, p_o)


def wall_radii(cyl, fractions=(0.0, 0.13, 0.37, 0.5, 0.71, 0.9, 1.0)):
    a, b = cyl.inner_radius, cyl.outer_radius
    return [a + (b - a) * f for f in fractions]


class TestRadialEquilibrium:
    """判据一：dσr/dr = (σθ − σr)/r 在筒壁内处处成立。"""

    @pytest.mark.parametrize(
        "a,b,p_i,p_o",
        [
            (100.0, 200.0, 10.0, 0.0),
            (50.0, 75.0, 8.0, 3.0),
            (100.0, 110.0, 25.0, 5.0),
        ],
    )
    def test_equilibrium_by_finite_difference(self, a, b, p_i, p_o):
        cyl = make_cylinder(a, b, p_i, p_o)
        for r in wall_radii(cyl, fractions=(0.13, 0.37, 0.5, 0.71, 0.9)):
            h = 1e-5 * r
            dsr_dr = (cyl.radial_stress(r + h) - cyl.radial_stress(r - h)) / (
                2.0 * h
            )
            rhs = (cyl.hoop_stress(r) - cyl.radial_stress(r)) / r
            assert math.isclose(dsr_dr, rhs, rel_tol=1e-6, abs_tol=1e-9), (
                f"r={r}: dσr/dr={dsr_dr} != (σθ−σr)/r={rhs}"
            )


class TestThinWallLimit:
    """判据二：薄壁极限对照，且厚壁时明显偏离薄壁估计。"""

    def test_near_thin_wall_matches_thin_wall_formula(self):
        # 厚比 K = b/a = 1.1，已接近薄壁
        a, b, p = 100.0, 110.0, 5.0
        cyl = make_cylinder(a, b, p, 0.0)
        thin = thin_wall_hoop_stress(p, 0.0, cyl.mean_radius, cyl.thickness)
        # 薄壁公式 p·r_m/t 与 Lamé 内壁环向应力应非常接近
        assert abs(cyl.hoop_stress(a) / thin - 1.0) < 0.01

    def test_thick_wall_inner_hoop_significantly_exceeds_thin_wall(self):
        # 厚比 K = 2.0，明显厚壁：内壁环向应力必须显著高于薄壁估计
        a, b, p = 100.0, 200.0, 5.0
        cyl = make_cylinder(a, b, p, 0.0)
        thin = thin_wall_hoop_stress(p, 0.0, cyl.mean_radius, cyl.thickness)
        assert cyl.hoop_stress(a) > 1.05 * thin
        # 精确比值 (a²+b²)/(b²−a²) / (r_m/t) = (5/3)/1.5 ≈ 1.111
        assert math.isclose(cyl.hoop_stress(a) / thin, 10.0 / 9.0, rel_tol=1e-12)


class TestBoundaryConditions:
    """判据三：仅内压时外壁径向应力恰为零；内壁径向应力等于负内压。"""

    def test_outer_radial_stress_vanishes_with_internal_pressure_only(self):
        cyl = make_cylinder(a=100.0, b=200.0, p_i=10.0, p_o=0.0)
        assert cyl.radial_stress(cyl.outer_radius) == pytest.approx(0.0, abs=1e-12)

    def test_wall_radial_stresses_match_applied_pressures(self):
        p_i, p_o = 12.0, 4.0
        cyl = make_cylinder(a=80.0, b=160.0, p_i=p_i, p_o=p_o)
        assert cyl.radial_stress(cyl.inner_radius) == pytest.approx(-p_i, abs=1e-10)
        assert cyl.radial_stress(cyl.outer_radius) == pytest.approx(-p_o, abs=1e-10)


class TestAgainstReferenceFormulas:
    """与独立手写的 Lamé 公式逐点对照（含 B 项符号：径向负、环向正）。"""

    def test_stress_field_matches_reference(self):
        a, b, p_i, p_o = 100.0, 200.0, 10.0, 3.0
        cyl = make_cylinder(a, b, p_i, p_o)
        for r in wall_radii(cyl):
            sr_ref, st_ref = reference_stresses(a, b, p_i, p_o, r)
            assert cyl.radial_stress(r) == pytest.approx(sr_ref, rel=1e-12)
            assert cyl.hoop_stress(r) == pytest.approx(st_ref, rel=1e-12)

    def test_inner_hoop_is_max_tensile_for_internal_pressure(self):
        # 仅内压时，内壁环向应力为全场最大拉应力
        cyl = make_cylinder(a=100.0, b=200.0, p_i=10.0, p_o=0.0)
        inner_hoop = cyl.hoop_stress(cyl.inner_radius)
        for r in wall_radii(cyl):
            assert inner_hoop >= cyl.hoop_stress(r)
            assert inner_hoop >= cyl.radial_stress(r)
        assert inner_hoop > 0.0


class TestLinearity:
    """联动一：内压单独放大两倍，整个应力场按比例放大两倍。"""

    def test_doubling_internal_pressure_doubles_stress_field(self):
        a, b, p_i, p_o = 100.0, 200.0, 10.0, 3.0
        cyl1 = make_cylinder(a, b, p_i, p_o)
        cyl2 = make_cylinder(a, b, 2.0 * p_i, p_o)
        for r in wall_radii(cyl1):
            # 外压部分不动，故 (σ2 − σ1) 应等于 (σ1 中内压贡献)
            sr_i = cyl1.radial_stress(r) - reference_stresses(a, b, 0.0, p_o, r)[0]
            st_i = cyl1.hoop_stress(r) - reference_stresses(a, b, 0.0, p_o, r)[1]
            assert cyl2.radial_stress(r) - cyl1.radial_stress(r) == pytest.approx(
                sr_i, rel=1e-12
            )
            assert cyl2.hoop_stress(r) - cyl1.hoop_stress(r) == pytest.approx(
                st_i, rel=1e-12
            )

    def test_doubling_all_pressures_doubles_everything(self):
        a, b, p_i, p_o = 100.0, 200.0, 10.0, 3.0
        cyl1 = make_cylinder(a, b, p_i, p_o)
        cyl2 = make_cylinder(a, b, 2.0 * p_i, 2.0 * p_o)
        for r in wall_radii(cyl1):
            assert cyl2.radial_stress(r) == pytest.approx(
                2.0 * cyl1.radial_stress(r), rel=1e-12
            )
            assert cyl2.hoop_stress(r) == pytest.approx(
                2.0 * cyl1.hoop_stress(r), rel=1e-12
            )


class TestHydrostaticLimit:
    """联动二：叠加与内压相等的外压 → 均匀静水压状态，偏应力下降。"""

    def test_equal_pressures_give_uniform_hydrostatic_state(self):
        a, b, p = 100.0, 200.0, 10.0
        cyl = make_cylinder(a, b, p_i=p, p_o=p)
        for r in wall_radii(cyl):
            assert cyl.radial_stress(r) == pytest.approx(-p, abs=1e-10)
            assert cyl.hoop_stress(r) == pytest.approx(-p, abs=1e-10)

    def test_deviatoric_stress_drops_when_outer_pressure_added(self):
        a, b, p = 100.0, 200.0, 10.0
        cyl_internal_only = make_cylinder(a, b, p_i=p, p_o=0.0)
        cyl_hydrostatic = make_cylinder(a, b, p_i=p, p_o=p)
        for r in wall_radii(cyl_internal_only):
            dev_before = abs(
                cyl_internal_only.hoop_stress(r) - cyl_internal_only.radial_stress(r)
            )
            dev_after = abs(
                cyl_hydrostatic.hoop_stress(r) - cyl_hydrostatic.radial_stress(r)
            )
            assert dev_after < 1e-10
            assert dev_after < dev_before

    def test_von_mises_vanishes_in_hydrostatic_state(self):
        # 闭口筒：σz = A = −p，三个主应力均为 −p，von Mises 恰为零
        cyl = make_cylinder(a=100.0, b=200.0, p_i=10.0, p_o=10.0)
        for r in wall_radii(cyl):
            assert cyl.von_mises(r, EndCondition.CLOSED, 0.3) == pytest.approx(
                0.0, abs=1e-9
            )

    def test_von_mises_reflects_axial_condition_in_hydrostatic_state(self):
        # 同样的内外等压下，开口筒 σz=0、平面应变 σz=2νA=−6，
        # 应力状态并非静水，von Mises 必须分别给出各自的非零值
        p, nu = 10.0, 0.3
        cyl = make_cylinder(a=100.0, b=200.0, p_i=p, p_o=p)
        for r in wall_radii(cyl):
            assert cyl.von_mises(r, EndCondition.OPEN, nu) == pytest.approx(p)
            assert cyl.von_mises(r, EndCondition.PLANE_STRAIN, nu) == pytest.approx(
                p * (1.0 - 2.0 * nu)
            )


class TestGeometricSimilarity:
    """联动三：几何同系数缩放后，以 r/a 为横坐标的应力分布不变。"""

    @pytest.mark.parametrize("scale", [0.5, 2.0, 3.7])
    def test_scaling_invariance(self, scale):
        a, b, p_i, p_o = 100.0, 160.0, 7.0, 2.0
        cyl1 = make_cylinder(a, b, p_i, p_o)
        cyl2 = make_cylinder(scale * a, scale * b, p_i, p_o)
        for ratio in (1.0, 1.1, 1.3, 1.6):
            r1 = ratio * a
            r2 = ratio * (scale * a)
            assert cyl2.radial_stress(r2) == pytest.approx(
                cyl1.radial_stress(r1), rel=1e-12
            )
            assert cyl2.hoop_stress(r2) == pytest.approx(
                cyl1.hoop_stress(r1), rel=1e-12
            )


class TestAxialConditions:
    """三种轴向约束分别处理：闭口 σz=A，开口 σz=0，平面应变 σz=2νA。"""

    def test_closed_end_axial_stress(self):
        a, b, p_i, p_o = 100.0, 200.0, 10.0, 0.0
        cyl = make_cylinder(a, b, p_i, p_o)
        expected = (p_i * a**2 - p_o * b**2) / (b**2 - a**2)
        assert cyl.axial_stress(EndCondition.CLOSED) == pytest.approx(expected)
        assert cyl.axial_stress(EndCondition.CLOSED) == pytest.approx(cyl.A)

    def test_open_end_axial_stress_is_zero(self):
        cyl = make_cylinder()
        assert cyl.axial_stress(EndCondition.OPEN) == 0.0

    def test_plane_strain_axial_stress_includes_poisson_effect(self):
        nu = 0.3
        cyl = make_cylinder(a=100.0, b=200.0, p_i=10.0, p_o=0.0)
        assert cyl.axial_stress(EndCondition.PLANE_STRAIN, nu) == pytest.approx(
            2.0 * nu * cyl.A
        )
        # 平面应变 σz 必须随 ν 变化（泊松效应确实计入）
        assert cyl.axial_stress(EndCondition.PLANE_STRAIN, 0.4) != pytest.approx(
            cyl.axial_stress(EndCondition.PLANE_STRAIN, 0.2)
        )

    def test_three_conditions_give_distinct_axial_stresses(self):
        cyl = make_cylinder(a=100.0, b=200.0, p_i=10.0, p_o=0.0)
        values = {
            cyl.axial_stress(EndCondition.CLOSED, 0.3),
            cyl.axial_stress(EndCondition.OPEN, 0.3),
            cyl.axial_stress(EndCondition.PLANE_STRAIN, 0.3),
        }
        assert len(values) == 3


class TestSolveLameConstants:
    """系数 A、B 与边界条件自洽。"""

    def test_constants_reproduce_boundary_conditions(self):
        a, b, p_i, p_o = 100.0, 200.0, 10.0, 3.0
        A, B = solve_lame_constants(a, b, p_i, p_o)
        assert A - B / a**2 == pytest.approx(-p_i, rel=1e-12)
        assert A - B / b**2 == pytest.approx(-p_o, rel=1e-12)
