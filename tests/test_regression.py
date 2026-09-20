"""应力场联动关系的回归锚点：

A. 内压单独加倍 → 整个应力场按比例加倍（线性）；
B. 叠加等额外压 → 趋近均匀静水压，偏应力 σθ−σr 下降；
C. 几何整体缩放 → 以 r/a 为横坐标的无量纲应力分布不变。
"""

from __future__ import annotations

import pytest

from app.domain import EndCondition
from app.lame import deviatoric_diameter, stress_profile
from app.validation import validate_cylinder


class TestLinearity:
    @pytest.mark.parametrize("end", ["closed", "open", "plane_strain"])
    def test_double_internal_pressure_doubles_field(self, end):
        poisson = 0.3 if end == "plane_strain" else None
        c1 = validate_cylinder(100, 200, 100, 25, end, poisson)
        c2 = validate_cylinder(100, 200, 200, 25, end, poisson)
        p1 = stress_profile(c1, 31)
        p2 = stress_profile(c2, 31)
        # 2·σ(p_i=100) 应等于 σ(p_i=200) − σ(p_i=0) 的贡献；直接验证差值场
        c0 = validate_cylinder(100, 200, 0, 25, end, poisson)
        p0 = stress_profile(c0, 31)
        for x0, x1, x2 in zip(p0, p1, p2):
            for attr in ("sigma_r", "sigma_theta", "sigma_z"):
                doubled = 2.0 * getattr(x1, attr) - getattr(x0, attr)
                assert getattr(x2, attr) == pytest.approx(doubled, rel=1e-12)

    def test_pure_double_without_external(self):
        # 外压为零时更直接：全场每个分量都精确加倍。
        c1 = validate_cylinder(80, 160, 120, 0, "closed")
        c2 = validate_cylinder(80, 160, 240, 0, "closed")
        for x1, x2 in zip(stress_profile(c1, 25), stress_profile(c2, 25)):
            assert x2.sigma_r == pytest.approx(2 * x1.sigma_r)
            assert x2.sigma_theta == pytest.approx(2 * x1.sigma_theta)
            assert x2.sigma_z == pytest.approx(2 * x1.sigma_z)

    def test_A_B_scale_linearly_with_pressure(self):
        from app.lame import solve_lame_coefficients

        c1 = validate_cylinder(100, 200, 100, 40, "open")
        c2 = validate_cylinder(100, 200, 150, 90, "open")
        k1, k2 = solve_lame_coefficients(c1), solve_lame_coefficients(c2)
        # Δp_i=50, Δp_o=50 相当于叠加 50 的均匀外/内压 ⇒ A 减 50、B 不变。
        assert k2.B == pytest.approx(k1.B)
        assert k2.A == pytest.approx(k1.A - 50.0)


class TestHydrostaticSuperposition:
    def test_equal_pressures_is_pure_hydrostatic(self):
        # p_i = p_o = p ⇒ B = 0，σr = σθ = σz(闭口) = -p，偏应力处处为零。
        p = 80.0
        cyl = validate_cylinder(100, 200, p, p, "closed")
        pts = stress_profile(cyl, 21)
        for pt in pts:
            assert pt.sigma_r == pytest.approx(-p)
            assert pt.sigma_theta == pytest.approx(-p)
            assert pt.sigma_z == pytest.approx(-p)
            assert deviatoric_diameter(pt) == pytest.approx(0.0, abs=1e-12)

    def test_adding_equal_external_pressure_reduces_deviator(self):
        # 在 (p_i, p_o)=(100,0) 上叠加 100 外压 → (100,100)：
        # A 不变偏置下偏应力应大幅下降至零（叠加静水压不改变偏应力的严格
        # 说法针对“同时叠加相等内、外压”；此处叠加的是外压，B 直接归零）。
        c_before = validate_cylinder(100, 200, 100, 0, "open")
        c_after = validate_cylinder(100, 200, 100, 100, "open")
        before = stress_profile(c_before, 21)
        after = stress_profile(c_after, 21)
        for b, a in zip(before, after):
            assert abs(deviatoric_diameter(a)) < abs(deviatoric_diameter(b))
            assert abs(deviatoric_diameter(a)) < 1e-12

    def test_hydrostatic_superposition_principle(self):
        # 叠加原理：σ(p_i,p_o) = σ(p_i,0) + σ(0,p_o)。
        base = validate_cylinder(100, 200, 130, 60, "plane_strain", 0.3)
        part_i = validate_cylinder(100, 200, 130, 0, "plane_strain", 0.3)
        part_o = validate_cylinder(100, 200, 0, 60, "plane_strain", 0.3)
        for s, si, so in zip(
            stress_profile(base, 15),
            stress_profile(part_i, 15),
            stress_profile(part_o, 15),
        ):
            assert s.sigma_r == pytest.approx(si.sigma_r + so.sigma_r)
            assert s.sigma_theta == pytest.approx(si.sigma_theta + so.sigma_theta)
            assert s.sigma_z == pytest.approx(si.sigma_z + so.sigma_z)


class TestGeometricScaleInvariance:
    @pytest.mark.parametrize("scale", [0.1, 3.7, 100.0])
    def test_dimensionless_profile_invariant(self, scale):
        # 内外半径与壁厚同乘 s 后，以 x = r/a 为横坐标、σ 同量纲看分布不变。
        c1 = validate_cylinder(100, 200, 100, 20, "closed")
        a2, b2 = 100 * scale, 200 * scale
        c2 = validate_cylinder(a2, b2, 100, 20, "closed")
        n = 41
        p1, p2 = stress_profile(c1, n), stress_profile(c2, n)
        for x1, x2 in zip(p1, p2):
            x_coord_1 = x1.r / c1.a
            x_coord_2 = x2.r / c2.a
            assert x_coord_1 == pytest.approx(x_coord_2)
            assert x1.sigma_r == pytest.approx(x2.sigma_r)
            assert x1.sigma_theta == pytest.approx(x2.sigma_theta)
            assert x1.sigma_z == pytest.approx(x2.sigma_z)

    def test_coefficients_B_scales_as_length_squared(self):
        from app.lame import solve_lame_coefficients

        c1 = validate_cylinder(100, 200, 100, 10, "open")
        c2 = validate_cylinder(300, 600, 100, 10, "open")
        k1, k2 = solve_lame_coefficients(c1), solve_lame_coefficients(c2)
        assert k1.A == pytest.approx(k2.A)              # A 是应力量，与尺度无关
        assert k2.B == pytest.approx(k1.B * 9.0)       # B 带长度平方量纲
