"""HTTP 层端到端测试：接口契约、错误结构、屈服判定、并发隔离。"""

from __future__ import annotations

import concurrent.futures

import pytest

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


VALID = {
    "a": 100.0,
    "b": 200.0,
    "p_i": 100.0,
    "p_o": 0.0,
    "end_condition": "closed",
}


class TestStressEndpoint:
    def test_basic_response_contract(self):
        resp = client.post("/api/v1/stress", json=VALID)
        assert resp.status_code == 200
        body = resp.json()
        assert body["hoop_stress"]["inner"] == pytest.approx(166.6666667)
        assert body["hoop_stress"]["outer"] == pytest.approx(66.6666667)
        assert body["inner_wall"]["sigma_r"] == -100.0
        assert body["outer_wall"]["sigma_r"] == 0.0
        assert set(body["lame_coefficients"]) == {"A", "B"}
        assert body["inner_wall"]["sigma_z"] == pytest.approx(33.3333333)
        assert "thin_wall_reference" in body

    def test_at_radius_query(self):
        resp = client.post("/api/v1/stress", json={**VALID, "r": 150.0})
        assert resp.status_code == 200
        at = resp.json()["at_radius"]
        # A=33.333, B=1.3333e6, r=150: σr=33.333-59.259=-25.926,
        # σθ=33.333+59.259=92.593
        assert at["sigma_r"] == pytest.approx(-25.9259259)
        assert at["sigma_theta"] == pytest.approx(92.5925926)
        assert "von_mises" in at

    def test_radius_outside_wall_422(self):
        resp = client.post("/api/v1/stress", json={**VALID, "r": 999.0})
        assert resp.status_code == 422
        detail = resp.json()["details"]
        assert any(d["field"] == "r" for d in detail)
        assert all(d["reason"] for d in detail)

    def test_yield_check_not_yielded(self):
        resp = client.post(
            "/api/v1/stress", json={**VALID, "yield_strength": 250.0}
        )
        body = resp.json()["yield_check"]
        assert body["yielded"] is False
        assert body["criterion"] == "von_mises"
        assert body["max_von_mises"] > 0

    def test_yield_check_yielded(self):
        resp = client.post(
            "/api/v1/stress", json={**VALID, "yield_strength": 50.0}
        )
        body = resp.json()["yield_check"]
        assert body["yielded"] is True
        # 手工：内壁 (σr,σθ,σz)=(-100,166.67,33.33)
        # vm = sqrt(0.5[(-266.67)^2+(133.33)^2+(133.33)^2]) ≈ 152.75...
        assert body["max_von_mises"] > 150.0

    def test_open_end_no_axial(self):
        resp = client.post(
            "/api/v1/stress", json={**VALID, "end_condition": "open"}
        )
        assert resp.json()["inner_wall"]["sigma_z"] == 0.0

    def test_plane_strain_with_poisson(self):
        resp = client.post(
            "/api/v1/stress",
            json={**VALID, "end_condition": "plane_strain", "poisson": 0.3},
        )
        # σz = 2νA = 0.6·33.333 = 20.0
        assert resp.json()["inner_wall"]["sigma_z"] == 20.0

    def test_plane_strain_missing_poisson_422(self):
        resp = client.post(
            "/api/v1/stress", json={**VALID, "end_condition": "plane_strain"}
        )
        assert resp.status_code == 422
        assert resp.json()["error"] == "validation_failed"

    def test_bad_geometry_422(self):
        resp = client.post("/api/v1/stress", json={**VALID, "b": 50.0})
        assert resp.status_code == 422
        assert resp.json()["details"][0]["field"] == "b"

    def test_negative_pressure_422(self):
        resp = client.post("/api/v1/stress", json={**VALID, "p_i": -10.0})
        assert resp.status_code == 422

    def test_unknown_field_rejected(self):
        resp = client.post("/api/v1/stress", json={**VALID, "extra": 1})
        assert resp.status_code == 422


class TestProfileEndpoint:
    def test_profile_points(self):
        resp = client.post(
            "/api/v1/stress/profile", json={**VALID, "points": 11}
        )
        assert resp.status_code == 200
        body = resp.json()
        assert body["points"] == 11
        assert len(body["radii"]) == 11
        assert len(body["sigma_r"]) == 11 and len(body["sigma_theta"]) == 11
        assert body["radii"][0] == 100.0 and body["radii"][-1] == 200.0
        # 首末点必须满足边界条件，证明逐点来自 Lamé 解而非写死曲线。
        assert body["sigma_r"][0] == -100.0
        assert body["sigma_r"][-1] == 0.0
        assert body["sigma_theta"][0] == pytest.approx(166.6666667)
        assert body["sigma_theta"][-1] == pytest.approx(66.6666667)

    def test_profile_default_points(self):
        resp = client.post("/api/v1/stress/profile", json=VALID)
        assert len(resp.json()["radii"]) == 21

    def test_profile_invalid_points(self):
        resp = client.post(
            "/api/v1/stress/profile", json={**VALID, "points": 1}
        )
        assert resp.status_code == 422


class TestExampleEndpoint:
    def test_example_is_hand_checkable(self):
        resp = client.get("/api/v1/examples/internal-pressure")
        assert resp.status_code == 200
        body = resp.json()
        assert body["hoop_stress"]["inner"] == pytest.approx(166.6666667)
        assert body["hoop_stress"]["outer"] == pytest.approx(66.6666667)
        assert body["example_notes"]["max_tensile_location"] == "inner_wall hoop"
        assert body["yield_check"]["yielded"] is False


class TestHealthAndConcurrency:
    def test_health(self):
        assert client.get("/health").json() == {"status": "ok"}

    def test_concurrent_requests_are_isolated(self):
        # 不同几何/载荷的请求并发打入，结果必须各自正确，无中间量串扰。
        cases = [
            ({**VALID}, 166.6666667),
            ({**VALID, "a": 50.0, "b": 75.0}, 260.0),  # p(75²+50²)/(75²-50²)
            ({**VALID, "p_i": 200.0}, 333.3333333),
            ({**VALID, "p_o": 100.0, "end_condition": "open"}, -100.0),  # 等内外压=静水压
        ]

        def call(case):
            payload, expected_inner = case
            resp = client.post("/api/v1/stress", json=payload)
            assert resp.status_code == 200
            assert resp.json()["hoop_stress"]["inner"] == pytest.approx(expected_inner)
            return True

        with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:
            results = list(pool.map(call, cases * 12))
        assert all(results)

    def test_concurrent_profile_and_stress(self):
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            futs = [
                pool.submit(lambda: client.post("/api/v1/stress", json=VALID).json()),
                pool.submit(
                    lambda: client.post(
                        "/api/v1/stress/profile", json=VALID
                    ).json()
                ),
            ] * 10
            for f in futs:
                out = f.result()
                assert "hoop_stress" in out or "sigma_r" in out
