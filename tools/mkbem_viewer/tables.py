from __future__ import annotations

import csv
from pathlib import Path

from .models import ResultStep


BOUNDARY_RESULT_HEADERS = [
    "index", "time_year", "x", "y", "ux", "uy", "temperature", "pressure",
    "traction_x", "traction_y", "thermal_flux", "flow_rate",
]

INTERIOR_RESULT_HEADERS = [
    "index", "time_year", "x", "y", "ux", "uy", "temperature", "pressure",
    "stress_x", "stress_y", "stress_xy", "thermal_flux_x", "thermal_flux_y", "flow_rate_x", "flow_rate_y",
]


def boundary_result_rows(step: ResultStep) -> list[list[float | int]]:
    rows: list[list[float | int]] = []
    for index, row in enumerate(step.boundary):
        rows.append([
            index,
            step.time_year,
            row.x,
            row.y,
            row.ux,
            row.uy,
            row.temperature,
            row.pressure,
            row.traction_x,
            row.traction_y,
            row.thermal_flux,
            row.flow_rate,
        ])
    return rows


def interior_result_rows(step: ResultStep) -> list[list[float | int | str]]:
    rows: list[list[float | int | str]] = []
    stress_by_index = {index: row for index, row in enumerate(step.interior_stress)}
    for index, row in enumerate(step.interior_primary):
        stress = stress_by_index.get(index)
        rows.append([
            index,
            step.time_year,
            row.x,
            row.y,
            row.ux,
            row.uy,
            row.temperature,
            row.pressure,
            stress.stress_x if stress else "",
            stress.stress_y if stress else "",
            stress.stress_xy if stress else "",
            stress.thermal_flux_x if stress else "",
            stress.thermal_flux_y if stress else "",
            stress.flow_rate_x if stress else "",
            stress.flow_rate_y if stress else "",
        ])
    return rows


def result_table(dataset: str, step: ResultStep) -> tuple[list[str], list[list[float | int | str]]]:
    if dataset == "boundary":
        return BOUNDARY_RESULT_HEADERS, boundary_result_rows(step)
    if dataset == "interior":
        return INTERIOR_RESULT_HEADERS, interior_result_rows(step)
    raise ValueError(f"Unknown result dataset: {dataset}")


def write_result_csv(path: str | Path, dataset: str, step: ResultStep) -> None:
    headers, rows = result_table(dataset, step)
    with Path(path).open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        writer.writerows(rows)
