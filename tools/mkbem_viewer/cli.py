from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from .input_parser import parse_input_file
from .result_parser import parse_result_file
from .tables import write_result_csv


def _default_result_path(input_path: Path) -> Path:
    return input_path.with_suffix(".out")


def build_summary(input_path: Path, result_path: Path | None) -> dict[str, object]:
    model = parse_input_file(input_path)
    summary: dict[str, object] = {
        "title": model.title,
        "input": str(model.path),
        "boundary_nodes": model.boundary_count,
        "interior_points": model.interior_point_count,
        "internal_nodes": model.internal_node_count,
        "internal_elements": model.internal_element_count,
        "input_time_steps": model.time_step_count,
        "sources": [asdict(source) for source in model.sources],
        "material": asdict(model.material),
    }

    if result_path is None:
        result_path = _default_result_path(input_path)

    if result_path.exists():
        result = parse_result_file(result_path, boundary_count=model.boundary_count, interior_count=model.interior_point_count)
        summary.update(
            {
                "result": str(result.path),
                "result_time_steps": result.time_step_count,
                "time_years": [step.time_year for step in result.steps],
                "boundary_rows_per_step": [len(step.boundary) for step in result.steps],
                "interior_rows_per_step": [len(step.interior_primary) for step in result.steps],
            }
        )
    else:
        summary["result"] = None

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Inspect MkBEM input and result files.")
    parser.add_argument("--input", required=True, type=Path, help="MkBEM .dat input file")
    parser.add_argument("--result", type=Path, help="MkBEM .out result file. Defaults to input path with .out extension.")
    parser.add_argument("--json", action="store_true", help="Print machine-readable JSON")
    parser.add_argument("--view", action="store_true", help="Open the PyQt/PyVista viewer after validating the files")
    parser.add_argument("--export-csv", type=Path, help="Export one result time step to CSV")
    parser.add_argument("--dataset", choices=["boundary", "interior"], default="boundary", help="Dataset to export when --export-csv is used")
    parser.add_argument("--step", type=int, default=-1, help="Result time-step index to export. Defaults to the last step")
    args = parser.parse_args()

    summary = build_summary(args.input, args.result)

    if args.export_csv is not None:
        result_path = args.result or _default_result_path(args.input)
        model = parse_input_file(args.input)
        result = parse_result_file(result_path, boundary_count=model.boundary_count, interior_count=model.interior_point_count)
        step_index = args.step if args.step >= 0 else result.time_step_count - 1
        if step_index < 0 or step_index >= result.time_step_count:
            raise SystemExit(f"Invalid step index {step_index}; result has {result.time_step_count} steps")
        write_result_csv(args.export_csv, args.dataset, result.steps[step_index])
        print(f"Exported {args.dataset} step {step_index} to {args.export_csv}")
        return 0

    if args.view:
        from .viewer import main as viewer_main

        viewer_args = ["--input", str(args.input)]
        if args.result is not None:
            viewer_args.extend(["--result", str(args.result)])
        import sys

        original_argv = sys.argv
        try:
            sys.argv = [original_argv[0], *viewer_args]
            return viewer_main()
        finally:
            sys.argv = original_argv

    if args.json:
        print(json.dumps(summary, indent=2))
        return 0

    print(f"Case: {summary['title']}")
    print(f"Input: {summary['input']}")
    print(f"Boundary nodes: {summary['boundary_nodes']}")
    print(f"Interior points: {summary['interior_points']}")
    print(f"Input time steps: {summary['input_time_steps']}")
    print(f"Sources: {len(summary['sources'])}")
    if summary.get("result"):
        print(f"Result: {summary['result']}")
        print(f"Result time steps: {summary['result_time_steps']}")
        print(f"Time years: {summary['time_years']}")
    else:
        print("Result: not found")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
