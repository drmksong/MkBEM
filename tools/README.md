# MkBEM Python Utilities

Stage 1 provides parsers and a command-line summary tool for MkBEM input and result files. Stage 2 adds a PyQt/PyVista viewer for model geometry and time-step result fields.

## Environment

The requested conda environment is named `bem` and uses Python 3.11.

```sh
conda activate bem
conda env update -n bem -f tools/environment.yml
```

For a pip-only setup inside an existing Python environment:

```sh
python -m pip install -r tools/requirements.txt
```

## CLI Summary

From the `02_MkBEM` directory:

```sh
python -m tools.mkbem_viewer.cli --input tests/THM_TeST2.dat --result tests/THM_TeST2.out
```

JSON output is also available:

```sh
python -m tools.mkbem_viewer.cli --input tests/THM_TeST2.dat --result tests/THM_TeST2.out --json
```

Export the current benchmark result table to CSV from the command line:

```sh
python -m tools.mkbem_viewer.cli --input tests/THM_TeST2.dat --result tests/THM_TeST2.out --dataset interior --step 4 --export-csv interior_step_004.csv
```

## Visualization Direction

Use PyQt with PyVista/PyVistaQt for the next stage. The current `bem` conda environment provides PyQt5 through conda-forge's `pyqt` package. PyVista is a good fit because MkBEM data can be represented as boundary polylines, interior point clouds, scalar fields, vector glyphs, and later VTK exports for ParaView.

## Viewer

From the `02_MkBEM` directory:

```sh
conda activate bem
python -m tools.mkbem_viewer.viewer --input tests/THM_TeST2.dat --result tests/THM_TeST2.out
```

The viewer currently supports:

- boundary geometry and interior point display
- boundary coordinate table with node index, coordinate, next-node connectivity, element midpoint, and boundary-condition code/value columns
- selected boundary table rows highlighted in the model view
- result table for the selected time step and dataset
- CSV export for the selected result table
- time-step selection
- boundary fields: displacement, temperature, pressure, traction, thermal flux, flow rate
- interior fields: displacement, temperature, pressure, stress, thermal flux, flow rate
- displacement vector overlays
