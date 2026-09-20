"""HTTP 请求/响应的 Pydantic 模型。

压力字段在模型层就约束为非负（受压为正的约定），
几何与泊松比的细粒度规则交由 app.validation 统一执行，
以便模型校验与领域校验返回同构的结构化错误。
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from .domain import EndCondition

EndConditionLiteral = Literal["closed", "open", "plane_strain"]
# 与枚举保持单一事实来源，防止字面量与枚举漂移。
assert {v.value for v in EndCondition} == {"closed", "open", "plane_strain"}


class StressRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    a: float = Field(..., gt=0, description="内半径")
    b: float = Field(..., gt=0, description="外半径（必须 > a）")
    p_i: float = Field(..., ge=0, description="内壁压力，受压为正")
    p_o: float = Field(..., ge=0, description="外壁压力，受压为正")
    end_condition: EndConditionLiteral
    poisson: float | None = Field(None, description="泊松比，平面应变时必填")
    r: float | None = Field(None, gt=0, description="待查询的半径，需落在 [a,b]")
    yield_strength: float | None = Field(None, gt=0, description="屈服强度，可选")


class ProfileRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    a: float = Field(..., gt=0, description="内半径")
    b: float = Field(..., gt=0, description="外半径（必须 > a）")
    p_i: float = Field(..., ge=0, description="内壁压力，受压为正")
    p_o: float = Field(..., ge=0, description="外壁压力，受压为正")
    end_condition: EndConditionLiteral
    poisson: float | None = Field(None, description="泊松比，平面应变时必填")
    points: int = Field(21, ge=2, le=10_000, description="采样点数（含内外壁）")
