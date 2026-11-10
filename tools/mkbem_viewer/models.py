from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class MaterialProperties:
    mu: float
    nu: float
    k_t: float
    k_p: float
    beta_ts: float
    alpha_p: float
    initial_temperature: float
    eta: float
    rho: float
    cp: float
    bulk_modulus_input: float
    skempton_b: float


@dataclass(frozen=True)
class BoundaryCondition:
    ux_code: int
    ux_value: float
    uy_code: int
    uy_value: float
    temp_code: int
    temp_value: float
    pressure_code: int
    pressure_value: float


@dataclass(frozen=True)
class Source:
    kind: str
    strength: float
    x: float
    y: float


@dataclass(frozen=True)
class InternalElement:
    nodes: tuple[int, int, int, int]
    original_volumetric_strain: float
    original_temperature: float
    original_pressure: float


@dataclass(frozen=True)
class BemInput:
    path: Path
    title: str
    boundary_count: int
    interior_point_count: int
    boundary_group_count: int
    source_count: int
    internal_node_count: int
    internal_element_count: int
    time_step_count: int
    material: MaterialProperties
    boundary_nodes: list[tuple[float, float]]
    boundary_conditions: list[BoundaryCondition]
    interior_points: list[tuple[float, float]]
    time_steps: list[float]
    sources: list[Source]
    boundary_group_last_nodes: list[int] = field(default_factory=list)
    internal_nodes: list[tuple[float, float]] = field(default_factory=list)
    internal_elements: list[InternalElement] = field(default_factory=list)


@dataclass(frozen=True)
class BoundaryResult:
    x: float
    y: float
    ux: float
    uy: float
    temperature: float
    pressure: float
    traction_x: float
    traction_y: float
    thermal_flux: float
    flow_rate: float


@dataclass(frozen=True)
class InteriorPrimaryResult:
    x: float
    y: float
    ux: float
    uy: float
    temperature: float
    pressure: float


@dataclass(frozen=True)
class InteriorStressResult:
    x: float
    y: float
    stress_x: float
    stress_y: float
    stress_xy: float
    thermal_flux_x: float
    thermal_flux_y: float
    flow_rate_x: float
    flow_rate_y: float


@dataclass(frozen=True)
class ResultStep:
    time_year: float
    boundary: list[BoundaryResult]
    interior_primary: list[InteriorPrimaryResult]
    interior_stress: list[InteriorStressResult]


@dataclass(frozen=True)
class BemResult:
    path: Path
    steps: list[ResultStep]

    @property
    def time_step_count(self) -> int:
        return len(self.steps)
