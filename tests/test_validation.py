"""入口输入校验测试：坏数据必须在计算前被结构化挡回。"""

from __future__ import annotations

import math

import pytest

from app.validation import (
    POISSON_MAX,
    POISSON_MIN,
    ValidationError,
    validate_cylinder,
    validate_profile_points,
    validate_radius,
    validate_yield_strength,
)


class TestGeometry:
    @pytest.mark.parametrize("a,b", [(200, 100), (100, 100), (0, 100), (-5, 10)])
    def test_bad_radii(self, a, b):
        with pytest.raises(ValidationError) as exc:
            validate_cylinder(a, b, 10, 0, "open")
        assert exc.value.field in ("a", "b")
        assert exc.value.message

    @pytest.mark.parametrize("bad", [float("nan"), float("inf"), "-inf", None, "100", True])
    def test_non_numeric_radius(self, bad):
        with pytest.raises(ValidationError):
            validate_cylinder(bad, 200, 10, 0, "open")


class TestPressureSign:
    @pytest.mark.parametrize("p", [-1.0, -1e-9])
    def test_negative_pressure_rejected(self, p):
        # 约定受压为正；负值属于未按约定给符号，入口挡回。
        with pytest.raises(ValidationError) as exc:
            validate_cylinder(100, 200, p, 0, "open")
        assert exc.value.field == "p_i"
        with pytest.raises(ValidationError) as exc:
            validate_cylinder(100, 200, 10, p, "open")
        assert exc.value.field == "p_o"

    @pytest.mark.parametrize("bad", [float("nan"), float("inf"), None])
    def test_non_finite_pressure_rejected(self, bad):
        with pytest.raises(ValidationError):
            validate_cylinder(100, 200, bad, 0, "open")

    def test_zero_pressures_allowed(self):
        cyl = validate_cylinder(100, 200, 0, 0, "open")
        assert cyl.p_i == 0 and cyl.p_o == 0


class TestEndConditionAndPoisson:
    def test_unknown_end_condition(self):
        with pytest.raises(ValidationError) as exc:
            validate_cylinder(100, 200, 10, 0, "welded")
        assert exc.value.field == "end_condition"

    def test_plane_strain_requires_poisson(self):
        with pytest.raises(ValidationError) as exc:
            validate_cylinder(100, 200, 10, 0, "plane_strain", None)
        assert exc.value.field == "poisson"

    @pytest.mark.parametrize("nu", [0.5, -1.0, 0.6, -2.0, float("nan"), float("inf")])
    def test_unphysical_poisson_rejected(self, nu):
        with pytest.raises(ValidationError) as exc:
            validate_cylinder(100, 200, 10, 0, "plane_strain", nu)
        assert exc.value.field == "poisson"

    @pytest.mark.parametrize("nu", [-0.99, 0.0, 0.25, 0.4999])
    def test_physical_poisson_accepted(self, nu):
        cyl = validate_cylinder(100, 200, 10, 0, "plane_strain", nu)
        assert cyl.poisson == nu

    def test_poisson_optional_for_closed_and_open(self):
        assert validate_cylinder(100, 200, 10, 0, "closed").poisson is None
        assert validate_cylinder(100, 200, 10, 0, "open").poisson is None


class TestRadiusAndAuxiliary:
    def test_query_radius_outside_wall(self):
        cyl = validate_cylinder(100, 200, 10, 0, "open")
        with pytest.raises(ValidationError) as exc:
            validate_radius(201, cyl)
        assert exc.value.field == "r"
        with pytest.raises(ValidationError):
            validate_radius(99, cyl)

    @pytest.mark.parametrize("n", [1, 0, -3, 10001, 2.5])
    def test_bad_profile_points(self, n):
        with pytest.raises(ValidationError):
            validate_profile_points(n)

    def test_bad_yield_strength(self):
        with pytest.raises(ValidationError):
            validate_yield_strength(0)
        with pytest.raises(ValidationError):
            validate_yield_strength(-10)


def test_error_is_structured():
    try:
        validate_cylinder(10, 5, 1, 0, "open")
    except ValidationError as e:
        detail = e.to_detail()
        assert set(detail) == {"field", "reason"}
        assert isinstance(detail["reason"], str) and detail["reason"]
    else:
        pytest.fail("应当抛出 ValidationError")
