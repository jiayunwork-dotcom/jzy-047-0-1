"""示范算例：仅内压作用的闭口厚壁圆筒。

参数（单位：mm、MPa，手工可核对）：
    内半径 a = 100，外半径 b = 200（厚比 K = 2）
    内压 p_i = 60，外压 p_o = 0，闭口端

手算对照：
    A = p_i·a²/(b²−a²) = 60·10000/30000 = 20 MPa
    B = p_i·a²·b²/(b²−a²) = 60·10000·40000/30000 = 8×10⁵ MPa·mm²
    σθ(a) = A + B/a² = 20 + 80 = 100 MPa   ← 全场最大拉应力
    σθ(b) = A + B/b² = 20 + 20 = 40 MPa
    σr(a) = A − B/a² = −60 MPa（= −p_i）
    σr(b) = A − B/b² = 0 MPa（外壁自由）
    σz    = A = 20 MPa（闭口）

运行：python examples/demo_internal_pressure.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.lame import EndCondition, LameCylinder

A = 100.0  # 内半径 mm
B = 200.0  # 外半径 mm
P_I = 60.0  # 内压 MPa
P_O = 0.0  # 外压 MPa


def main() -> None:
    cylinder = LameCylinder.from_pressures(A, B, P_I, P_O)
    cond = EndCondition.CLOSED

    print("厚壁圆筒 Lamé 应力分析 —— 仅内压示范算例")
    print(f"  a = {A} mm, b = {B} mm, p_i = {P_I} MPa, p_o = {P_O} MPa, 闭口端")
    print(f"  Lamé 系数: A = {cylinder.A:g} MPa, B = {cylinder.B:g} MPa·mm²")
    print()
    print(f"  {'r (mm)':>8} {'σr (MPa)':>10} {'σθ (MPa)':>10} {'σz (MPa)':>10}")
    for i in range(11):
        r = A + (B - A) * i / 10
        s = cylinder.stress_at(r, cond)
        print(f"  {r:8.1f} {s.sigma_r:10.3f} {s.sigma_theta:10.3f} {s.sigma_z:10.3f}")
    print()

    inner_hoop = cylinder.hoop_stress(A)
    outer_hoop = cylinder.hoop_stress(B)
    print(f"  内壁环向应力 σθ(a) = {inner_hoop:g} MPa（手算 100）")
    print(f"  外壁环向应力 σθ(b) = {outer_hoop:g} MPa（手算 40）")
    print(f"  外壁径向应力 σr(b) = {cylinder.radial_stress(B):g} MPa（手算 0）")

    # 手工核对的关键断言：内壁环向应力为全场最大拉应力
    assert inner_hoop == 100.0
    assert outer_hoop == 40.0
    assert cylinder.radial_stress(B) == 0.0
    assert cylinder.radial_stress(A) == -P_I
    for i in range(101):
        r = A + (B - A) * i / 100
        assert inner_hoop >= cylinder.hoop_stress(r)
        assert inner_hoop >= cylinder.radial_stress(r)
    print("  手工核对通过：内壁环向应力 100 MPa 为全场最大拉应力。")


if __name__ == "__main__":
    main()
