"""HTTP 接口测试：功能正确性、输入校验挡回、并发请求互不干扰。"""

import math
from concurrent.futures import ThreadPoolExecutor

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

BASE_PAYLOAD = {
    "inner_radius": 100.0,
    "outer_radius": 200.0,
    "internal_pressure": 60.0,
    "external_pressure": 0.0,
    "end_condition": "closed",
    "poisson_ratio": 0.3,
}


def analytic(a, b, p_i, p_o, r):
    """独立手算 Lamé 解，用于对照接口返回。"""
    A = (p_i * a**2 - p_o * b**2) / (b**2 - a**2)
    B = (p_i - p_o) * a**2 * b**2 / (b**2 - a**2)
    return A - B / r**2, A + B / r**2


class TestStressEndpoint:
    def test_inner_and_outer_hoop_stresses(self):
        resp = client.post("/stress", json={**BASE_PAYLOAD, "query_radius": 100.0})
        assert resp.status_code == 200
        body = resp.json()
        # 手算：A = 60·100²/(200²−100²) = 20，B = 60·100²·200²/30000 = 8e5
        # σθ(a) = 20 + 8e5/1e4 = 100；σθ(b) = 20 + 8e5/4e4 = 40
        assert body["inner_hoop_stress"] == pytest.approx(100.0)
        assert body["outer_hoop_stress"] == pytest.approx(40.0)
        assert body["lame_A"] == pytest.approx(20.0)
        assert body["lame_B"] == pytest.approx(8.0e5)
        # 闭口筒轴向应力 σz = A
        assert body["axial_stress"] == pytest.approx(20.0)

    def test_query_point_matches_lame_solution(self):
        for r in (100.0, 137.5, 163.0, 200.0):
            resp = client.post("/stress", json={**BASE_PAYLOAD, "query_radius": r})
            assert resp.status_code == 200
            q = resp.json()["query"]
            sr_ref, st_ref = analytic(100.0, 200.0, 60.0, 0.0, r)
            assert q["r"] == pytest.approx(r)
            assert q["sigma_r"] == pytest.approx(sr_ref, rel=1e-9)
            assert q["sigma_theta"] == pytest.approx(st_ref, rel=1e-9)

    def test_end_conditions_change_axial_stress(self):
        axial = {}
        for cond in ("closed", "open", "plane_strain"):
            resp = client.post(
                "/stress",
                json={**BASE_PAYLOAD, "end_condition": cond, "query_radius": 150.0},
            )
            assert resp.status_code == 200
            axial[cond] = resp.json()["axial_stress"]
        assert axial["closed"] == pytest.approx(20.0)
        assert axial["open"] == pytest.approx(0.0)
        assert axial["plane_strain"] == pytest.approx(2.0 * 0.3 * 20.0)

    def test_yield_check(self):
        # 内壁 von Mises 最大：σr=−60, σθ=100, σz=20 → vm = 80·√3 ≈ 138.6
        resp = client.post(
            "/stress",
            json={**BASE_PAYLOAD, "query_radius": 100.0, "yield_strength": 200.0},
        )
        assert resp.status_code == 200
        body = resp.json()
        assert body["yielded"] is False
        assert body["max_von_mises"] == pytest.approx(body["query"]["von_mises"], rel=1e-3)

        resp = client.post(
            "/stress",
            json={**BASE_PAYLOAD, "query_radius": 100.0, "yield_strength": 100.0},
        )
        assert resp.json()["yielded"] is True

    def test_yield_fields_null_without_yield_strength(self):
        resp = client.post("/stress", json={**BASE_PAYLOAD, "query_radius": 100.0})
        body = resp.json()
        assert body["yield_strength"] is None
        assert body["yielded"] is None


class TestProfileEndpoint:
    def test_profile_points_match_lame_solution_pointwise(self):
        resp = client.post("/profile", json={**BASE_PAYLOAD, "n_points": 41})
        assert resp.status_code == 200
        body = resp.json()
        points = body["points"]
        assert len(points) == 41
        # 等间距、覆盖整个壁厚
        assert points[0]["r"] == pytest.approx(100.0)
        assert points[-1]["r"] == pytest.approx(200.0)
        steps = {round(points[i + 1]["r"] - points[i]["r"], 9) for i in range(40)}
        assert len(steps) == 1
        # 每个点都必须等于 Lamé 解在该半径处的值（不许写死衰减形状）
        for p in points:
            sr_ref, st_ref = analytic(100.0, 200.0, 60.0, 0.0, p["r"])
            assert p["sigma_r"] == pytest.approx(sr_ref, rel=1e-9)
            assert p["sigma_theta"] == pytest.approx(st_ref, rel=1e-9)
        # 仅内压时环向应力沿壁厚严格单调下降
        hoops = [p["sigma_theta"] for p in points]
        assert all(h1 > h2 for h1, h2 in zip(hoops, hoops[1:]))

    def test_profile_reports_thin_wall_reference(self):
        resp = client.post("/profile", json={**BASE_PAYLOAD, "n_points": 10})
        body = resp.json()
        # 薄壁对照：p·r_m/t = 60·150/100 = 90
        assert body["thin_wall_hoop_stress"] == pytest.approx(90.0)
        # Lamé 内壁环向应力 100 与薄壁估计 90 明显不同 → 不是薄壁近似冒充
        assert body["points"][0]["sigma_theta"] == pytest.approx(100.0)


class TestInputValidation:
    def test_outer_radius_not_greater_than_inner(self):
        resp = client.post(
            "/stress",
            json={**BASE_PAYLOAD, "outer_radius": 100.0, "query_radius": 100.0},
        )
        assert resp.status_code == 400
        err = resp.json()["error"]
        assert err["code"] == "INVALID_GEOMETRY"
        assert "外半径" in err["message"]

    @pytest.mark.parametrize("field", ["inner_radius", "outer_radius"])
    def test_non_positive_radius_rejected(self, field):
        resp = client.post(
            "/stress", json={**BASE_PAYLOAD, field: -5.0, "query_radius": 100.0}
        )
        assert resp.status_code == 400
        assert resp.json()["error"]["code"] == "INVALID_GEOMETRY"

    @pytest.mark.parametrize("field", ["internal_pressure", "external_pressure"])
    def test_missing_pressure_rejected(self, field):
        payload = {**BASE_PAYLOAD, "query_radius": 100.0}
        del payload[field]
        resp = client.post("/stress", json=payload)
        assert resp.status_code == 422
        assert resp.json()["error"]["code"] == "INVALID_REQUEST"

    @pytest.mark.parametrize("nu", [0.6, 1.0, -1.5])
    def test_unphysical_poisson_ratio_rejected(self, nu):
        resp = client.post(
            "/stress",
            json={
                **BASE_PAYLOAD,
                "end_condition": "plane_strain",
                "poisson_ratio": nu,
                "query_radius": 150.0,
            },
        )
        assert resp.status_code == 400
        err = resp.json()["error"]
        assert err["code"] == "INVALID_POISSON_RATIO"
        assert "泊松比" in err["message"]

    def test_query_radius_out_of_wall_rejected(self):
        resp = client.post("/stress", json={**BASE_PAYLOAD, "query_radius": 250.0})
        assert resp.status_code == 400
        assert resp.json()["error"]["code"] == "QUERY_RADIUS_OUT_OF_RANGE"

    def test_error_response_is_structured(self):
        resp = client.post("/stress", json={**BASE_PAYLOAD, "inner_radius": 0.0,
                                            "query_radius": 100.0})
        assert resp.status_code == 400
        err = resp.json()["error"]
        assert set(err) >= {"code", "message"}
        assert err["message"]


class TestConcurrency:
    """并发请求各自独立：不同参数的响应互不串扰。"""

    def test_concurrent_requests_are_independent(self):
        cases = [
            (100.0, 200.0, 10.0 * k, 0.0)
            for k in range(1, 17)
        ]

        def run(case):
            a, b, p_i, p_o = case
            resp = client.post(
                "/stress",
                json={
                    "inner_radius": a,
                    "outer_radius": b,
                    "internal_pressure": p_i,
                    "external_pressure": p_o,
                    "end_condition": "closed",
                    "query_radius": a,
                },
            )
            assert resp.status_code == 200
            return p_i, resp.json()

        with ThreadPoolExecutor(max_workers=8) as pool:
            results = list(pool.map(run, cases))

        for p_i, body in results:
            # 仅内压、K=2：σθ(a) = p_i·(a²+b²)/(b²−a²) = p_i·5/3
            assert body["inner_hoop_stress"] == pytest.approx(p_i * 5.0 / 3.0)
            assert body["query"]["sigma_r"] == pytest.approx(-p_i)


class TestHealth:
    def test_health(self):
        resp = client.get("/health")
        assert resp.status_code == 200
        assert resp.json() == {"status": "ok"}
