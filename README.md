# 厚壁圆筒 Lamé 弹性应力计算服务

给定厚壁圆筒的内半径 `a`、外半径 `b`、内壁压力 `p_i`、外壁压力 `p_o`
以及两端轴向约束方式（闭口 / 开口 / 平面应变），计算沿壁厚任意半径处的
**径向应力 σr、环向应力 σθ、轴向应力 σz**，并单独报出内外壁环向应力。
服务只做这一件事：无持久化、无台账、无任何其他力学模型。

- 运行时：**Python 3.11**（容器基础镜像 `python:3.11-slim-bookworm`）
- Web 层：FastAPI + Uvicorn
- 约定：**压力受压为正**输入；**应力拉为正、压为负**输出；几何量同量纲即可

## 物理模型

轴对称弹性解取 Lamé 形式（`σr + σθ` 沿壁厚为常数）：

```
σr(r) = A − B/r²
σθ(r) = A + B/r²        ← B 项与径向符号相反，切勿写同号
```

边界条件 `σr(a) = −p_i`、`σr(b) = −p_o`（压力受压为正 → 边界压应力为负）：

```
A = (a²·p_i − b²·p_o) / (b² − a²)
B = a²·b²·(p_i − p_o) / (b² − a²)
```

轴向应力按端部约束分别处理：

| 约束 | σz |
|---|---|
| 闭口 `closed` | 端盖合力平衡：`σz = (p_i a² − p_o b²)/(b²−a²) = A` |
| 开口 `open` | 轴向自由：`σz = 0` |
| 平面应变 `plane_strain` | εz=0：`σz = ν(σr+σθ) = 2νA`（必须提供 `poisson`） |

σr、σθ 与端部条件无关；平面应变只通过 σz 体现泊松效应。

## 项目结构（模块按职责分开）

```
app/
  domain.py      领域数据模型（Cylinder / LaméCoeffs / StressPoint 等）
  validation.py  入口输入校验（几何、压力符号、泊松比、半径范围）
  lame.py        边界条件求解 A/B、Lamé 应力场、σz、von Mises（推导写在模块文档串）
  thin_wall.py   薄壁薄膜公式对照（p·r_m/t）与厚径比
  service.py     无状态计算编排（所有中间量均为请求内局部量）
  schemas.py     HTTP 请求模型
  main.py        FastAPI 路由、结构化错误、示范算例
tests/
  test_lame.py        系数、边界、符号、轴向三条件、手工算例
  test_physics.py     三条物理判据（平衡方程/薄壁极限/外壁 σr=0）
  test_regression.py  线性、叠加静水压、几何缩放不变性
  test_validation.py  非法输入结构化拦截
  test_api.py         HTTP 契约、屈服判定、并发隔离
```

## 一步构建并运行（容器）

```bash
# 方式一：Docker
docker build -t thick-cylinder-api .
docker run --rm -p 8000:8000 thick-cylinder-api

# 方式二：docker compose
docker compose up --build
```

打开交互文档：<http://localhost:8000/docs>

容器内运行测试（镜像内置 test 目标）：

```bash
docker build --target test -t lame-test .
docker run --rm lame-test
```

## 本地运行（Python 3.11）

```bash
pip install -r requirements-dev.txt
uvicorn app.main:app --port 8000
pytest -q          # 97 项测试
```

## HTTP 接口

### `POST /api/v1/stress` — 单根筒完整计算

请求体：

```json
{
  "a": 100.0,
  "b": 200.0,
  "p_i": 100.0,
  "p_o": 0.0,
  "end_condition": "closed",
  "poisson": null,
  "r": 150.0,
  "yield_strength": 250.0
}
```

`r` 与 `yield_strength` 可选；`end_condition` 取
`closed` / `open` / `plane_strain`；平面应变时 `poisson` 必填。

响应含 `inner_wall`、`outer_wall`（各自 σr/σθ/σz/von Mises）、
顶层 `hoop_stress.{inner,outer}`、Lamé 系数 A/B、`at_radius`、
`thin_wall_reference`，以及给定屈服强度时的 `yield_check`
（按 von Mises 等效应力判定，`yielded` 为布尔值）。

### `POST /api/v1/stress/profile` — 沿壁厚等间距采样

```json
{ "a": 100.0, "b": 200.0, "p_i": 100.0, "p_o": 0.0,
  "end_condition": "closed", "points": 21 }
```

返回等长点列 `radii / sigma_r / sigma_theta / sigma_z`。
每个点都由 `σr=A−B/r²、σθ=A+B/r²` 逐点实时计算，不是预设曲线。

### `GET /api/v1/examples/internal-pressure` — 仅内压示范算例

`a=100, b=200, p_i=100, p_o=0`（闭口），可手工核对：

```
σθ(a) = p_i(b²+a²)/(b²−a²) = 100·50000/30000 = 166.67  ← 全场最大拉应力
σθ(b) = 2 p_i a²/(b²−a²)                       =  66.67
σr(a) = −100,   σr(b) = 0,   σz = 33.33
```

### 错误响应（HTTP 422，结构化）

```json
{
  "error": "validation_failed",
  "details": [ { "field": "b", "reason": "外半径 b 必须严格大于内半径 a（收到 a=200.0, b=100.0）" } ]
}
```

在任何计算之前拦截：`b ≤ a`、半径非正/非数值、压力为负（符号约定）、
未知端部条件、平面应变缺泊松比或 ν 不在 `(−1, 0.5)`、查询半径超出
`[a,b]`、采样点数越界、屈服强度非正。

## 自动化测试守住的物理判据

1. **径向平衡方程**：`dσr/dr = (σθ−σr)/r`，对多种载荷组合在壁厚 5 个
   半径点做中心差分，残差低于特征量级的 1e-6。
2. **薄壁极限**：`b/a = 1.1` 时 Lamé 内壁环向与 `Δp·r_m/t` 的相对偏差
   < 0.5%（解析偏差 `2(K²+1)/(K+1)² − 1`）；`b/a = 2` 时 Lamé 值
   166.67 显著高于薄膜值 150（>10%），且偏差随 K→1 单调趋零。
3. **自由外表面**：仅内压时多组厚径比下 `σr(b)` 恒为零（abs < 1e-12）。

联动回归锚点：

- 内压加倍（外压不变）→ 全场三向应力按线性叠加严格加倍；
- 叠加与内压相等的外压 → `B = 0`，σr=σθ=σz=−p（纯静水压），
  偏应力 `σθ−σr` 处处降为 0；并验证叠加原理 σ(p_i,p_o)=σ(p_i,0)+σ(0,p_o)；
- 几何整体缩放 → 以 `r/a` 为横坐标的应力分布完全不变（A 不变，
  B 按长度平方缩放）。

此外含 48 次混合参数的线程池并发请求测试，验证请求间中间量互不干扰。
