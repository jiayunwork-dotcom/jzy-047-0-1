# 厚壁圆筒 Lamé 应力分析服务

给定筒体内/外半径、内/外压力与轴向约束方式，按 Lamé 经典弹性解计算
厚壁圆筒的应力场，经 HTTP 对外提供。纯计算、无状态，可被强度校核流程
反复并发调用。

## 物理模型

符号约定：**拉应力为正**；压力以"压为正"输入，故壁面径向应力等于负压力。

```
σr(r) = A − B/r²        径向（B 项取负号）
σθ(r) = A + B/r²        环向（B 项取正号，与径向相反）

A = (p_i·a² − p_o·b²) / (b² − a²)
B = (p_i − p_o)·a²·b² / (b² − a²)
```

系数推导见 `app/boundary.py` 模块 docstring（由 σr(a) = −p_i、σr(b) = −p_o
两个边界条件联立解出）。

轴向约束三种情形分别处理：

| `end_condition` | 含义 | 轴向应力 σz |
|---|---|---|
| `closed` | 闭口（端盖承压） | σz = A |
| `open` | 开口（轴向自由） | σz = 0 |
| `plane_strain` | 平面应变（εz = 0） | σz = 2νA（计入泊松效应） |

## 项目结构

```
app/
  boundary.py    # Lamé 系数 A、B 的边界条件求解（含推导说明）
  lame.py        # Lamé 应力场：σr、σθ、σz、von Mises
  thinwall.py    # 薄壁极限对照（仅作参照，不参与对外结果）
  validation.py  # 输入合法性校验（计算入口前挡回）
  models.py      # 请求/响应模型
  main.py        # FastAPI 装配
examples/
  demo_internal_pressure.py   # 仅内压示范算例（可手工核对）
tests/
  test_physics.py  # 物理判据与联动关系
  test_api.py      # HTTP 接口、校验挡回、并发独立性
```

## 接口

### `POST /stress`

单根筒的完整应力计算。

```json
{
  "inner_radius": 100.0,
  "outer_radius": 200.0,
  "internal_pressure": 60.0,
  "external_pressure": 0.0,
  "end_condition": "closed",
  "poisson_ratio": 0.3,
  "query_radius": 100.0,
  "yield_strength": 250.0
}
```

返回：内外壁环向应力、Lamé 系数、轴向应力、指定半径处的 σr/σθ/σz/von Mises、
全壁厚最大 von Mises，以及（给出 `yield_strength` 时）是否屈服。
`yield_strength` 可省略，省略时 `yielded` 为 `null`。

### `POST /profile`

沿壁厚等间距采样（`n_points`，默认 50），逐点返回 Lamé 解的
σr(r)、σθ(r)、σz、von Mises，供上游绘制应力沿壁分布曲线；
附薄壁公式对照值 `thin_wall_hoop_stress`。

### `GET /health`

健康检查。

### 错误格式

所有校验失败返回结构化错误，不带着坏数据进入计算：

```json
{ "error": { "code": "INVALID_GEOMETRY", "message": "外半径必须严格大于内半径：..." } }
```

| 场景 | HTTP | code |
|---|---|---|
| 外半径 ≤ 内半径、半径非正值 | 400 | `INVALID_GEOMETRY` |
| 压力缺失/非有限数 | 422 / 400 | `INVALID_REQUEST` / `INVALID_PRESSURE` |
| 泊松比超出 −1 < ν ≤ 0.5 | 400 | `INVALID_POISSON_RATIO` |
| 查询半径越出壁厚 | 400 | `QUERY_RADIUS_OUT_OF_RANGE` |
| 屈服强度非正值 | 400 | `INVALID_YIELD_STRENGTH` |

## 运行

### 容器（一步构建并运行）

```bash
docker build -t lame-service .
docker run --rm -p 8000:8000 lame-service
# 文档：http://localhost:8000/docs
```

### 本地开发（Python 3.11）

```bash
pip install -r requirements.txt -r requirements-dev.txt
uvicorn app.main:app --reload
```

### 示范算例

```bash
python examples/demo_internal_pressure.py
```

a=100 mm、b=200 mm、p_i=60 MPa、闭口：内壁环向应力恰为 100 MPa，
是全场最大拉应力，可手工核对。

## 测试

```bash
pytest
```

覆盖的物理判据（全部自动化，逐条独立）：

1. 径向平衡方程 dσr/dr = (σθ − σr)/r 在壁内处处成立（数值差分抽查）；
2. 薄壁极限：厚比 1.1 时环向应力贴近薄壁公式 p·r_m/t；厚比 2 时内壁
   环向应力显著高于薄壁估计（防止以薄壁近似冒充 Lamé 解）；
3. 仅内压时外壁径向应力恰为零；
4. 线性：内压放大两倍，应力场按比例放大两倍；
5. 静水压：叠加相等外压后趋于均匀静水状态，偏应力 (σθ−σr) 降为零；
6. 相似性：几何同系数缩放后，以 r/a 为横坐标的应力分布不变；
7. 三种轴向约束分别给出正确且互不相同的 σz。

接口测试另覆盖：profile 曲线逐点等于 Lamé 解（非写死形状）、屈服判定、
全部校验挡回分支，以及 16 路并发请求结果互不串扰（服务无共享可变状态，
每个请求独立构造自己的 Lamé 解）。
