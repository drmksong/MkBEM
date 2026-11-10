from __future__ import annotations

import math
import re
from pathlib import Path

from .models import BemResult, BoundaryResult, InteriorPrimaryResult, InteriorStressResult, ResultStep

_TIME_RE = re.compile(r"^\s*Time\s*=\s*([-+0-9.eE]+)")
_NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def _numbers(line: str) -> list[float]:
    return [float(value) for value in _NUMBER_RE.findall(line)]


def _has_nonfinite(values: list[float]) -> bool:
    return any(not math.isfinite(value) for value in values)


def parse_result_file(path: str | Path, boundary_count: int | None = None, interior_count: int | None = None) -> BemResult:
    result_path = Path(path)
    lines = result_path.read_text(encoding="utf-8", errors="replace").splitlines()
    steps: list[ResultStep] = []
    index = 0

    while index < len(lines):
        time_match = _TIME_RE.match(lines[index])
        if not time_match:
            index += 1
            continue

        time_year = float(time_match.group(1))
        index += 1

        boundary: list[BoundaryResult] = []
        while index < len(lines) and "Interior Points Displacement" not in lines[index] and not _TIME_RE.match(lines[index]):
            values = _numbers(lines[index])
            if len(values) == 10:
                if _has_nonfinite(values):
                    raise ValueError(f"Non-finite boundary result near line {index + 1}")
                boundary.append(BoundaryResult(*values))
            index += 1

        if boundary_count is not None and len(boundary) != boundary_count:
            raise ValueError(f"Expected {boundary_count} boundary rows at time {time_year}, got {len(boundary)}")

        interior_primary: list[InteriorPrimaryResult] = []
        if index < len(lines) and "Interior Points Displacement" in lines[index]:
            index += 1
            if index < len(lines):
                index += 1
            while index < len(lines) and "Interior Points Stress" not in lines[index] and not _TIME_RE.match(lines[index]):
                values = _numbers(lines[index])
                if len(values) == 6:
                    if _has_nonfinite(values):
                        raise ValueError(f"Non-finite interior primary result near line {index + 1}")
                    interior_primary.append(InteriorPrimaryResult(*values))
                index += 1

        if interior_count is not None and len(interior_primary) != interior_count:
            raise ValueError(f"Expected {interior_count} interior primary rows at time {time_year}, got {len(interior_primary)}")

        interior_stress: list[InteriorStressResult] = []
        if index < len(lines) and "Interior Points Stress" in lines[index]:
            index += 1
            if index < len(lines):
                index += 1
            while index < len(lines) and not _TIME_RE.match(lines[index]):
                values = _numbers(lines[index])
                if len(values) == 9:
                    if _has_nonfinite(values):
                        raise ValueError(f"Non-finite interior stress result near line {index + 1}")
                    interior_stress.append(InteriorStressResult(*values))
                index += 1

        if interior_count is not None and len(interior_stress) != interior_count:
            raise ValueError(f"Expected {interior_count} interior stress rows at time {time_year}, got {len(interior_stress)}")

        steps.append(ResultStep(time_year=time_year, boundary=boundary, interior_primary=interior_primary, interior_stress=interior_stress))

    return BemResult(path=result_path, steps=steps)
