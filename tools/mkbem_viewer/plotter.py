from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pyvista as pv

from .models import BemInput, ResultStep


BOUNDARY_FIELDS = {
    "ux": "ux",
    "uy": "uy",
    "temperature": "temperature",
    "pressure": "pressure",
    "traction_x": "traction_x",
    "traction_y": "traction_y",
    "thermal_flux": "thermal_flux",
    "flow_rate": "flow_rate",
}

INTERIOR_FIELDS = {
    "ux": "ux",
    "uy": "uy",
    "temperature": "temperature",
    "pressure": "pressure",
    "stress_x": "stress_x",
    "stress_y": "stress_y",
    "stress_xy": "stress_xy",
    "thermal_flux_x": "thermal_flux_x",
    "thermal_flux_y": "thermal_flux_y",
    "flow_rate_x": "flow_rate_x",
    "flow_rate_y": "flow_rate_y",
}


def _points_2d(points: Iterable[tuple[float, float]]) -> np.ndarray:
    return np.array([(x, y, 0.0) for x, y in points], dtype=float)


def make_boundary_polyline(model: BemInput) -> pv.PolyData:
    points = _points_2d(model.boundary_nodes)
    mesh = pv.PolyData(points)
    if len(points) > 1:
        closed_indices = list(range(len(points))) + [0]
        mesh.lines = np.array([len(closed_indices), *closed_indices], dtype=np.int_)
    return mesh


def make_boundary_segment(model: BemInput, index: int) -> pv.PolyData:
    if not model.boundary_nodes:
        return pv.PolyData()
    if index < 0 or index >= len(model.boundary_nodes):
        return pv.PolyData()

    next_index = 0 if index == len(model.boundary_nodes) - 1 else index + 1
    mesh = pv.PolyData(_points_2d([model.boundary_nodes[index], model.boundary_nodes[next_index]]))
    mesh.lines = np.array([2, 0, 1], dtype=np.int_)
    return mesh


def make_interior_points(model: BemInput) -> pv.PolyData:
    return pv.PolyData(_points_2d(model.interior_points))


def make_boundary_results(step: ResultStep, field: str) -> pv.PolyData:
    if field not in BOUNDARY_FIELDS:
        raise ValueError(f"Unknown boundary field: {field}")
    points = _points_2d((row.x, row.y) for row in step.boundary)
    mesh = pv.PolyData(points)
    mesh[field] = np.array([getattr(row, BOUNDARY_FIELDS[field]) for row in step.boundary], dtype=float)
    mesh["displacement"] = np.array([(row.ux, row.uy, 0.0) for row in step.boundary], dtype=float)
    return mesh


def make_interior_results(step: ResultStep, field: str) -> pv.PolyData:
    if field not in INTERIOR_FIELDS:
        raise ValueError(f"Unknown interior field: {field}")

    if field in {"stress_x", "stress_y", "stress_xy", "thermal_flux_x", "thermal_flux_y", "flow_rate_x", "flow_rate_y"}:
        rows = step.interior_stress
    else:
        rows = step.interior_primary

    points = _points_2d((row.x, row.y) for row in rows)
    mesh = pv.PolyData(points)
    mesh[field] = np.array([getattr(row, INTERIOR_FIELDS[field]) for row in rows], dtype=float)

    if step.interior_primary and len(step.interior_primary) == len(rows):
        mesh["displacement"] = np.array([(row.ux, row.uy, 0.0) for row in step.interior_primary], dtype=float)
    return mesh
