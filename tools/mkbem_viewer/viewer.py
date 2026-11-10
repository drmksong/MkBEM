from __future__ import annotations

import argparse
import os
from dataclasses import asdict
from pathlib import Path

from PyQt5 import QtCore, QtWidgets

from .input_parser import parse_input_file
from .models import BemInput, BemResult
from .plotter import BOUNDARY_FIELDS, INTERIOR_FIELDS, make_boundary_polyline, make_boundary_results, make_boundary_segment, make_interior_points, make_interior_results
from .result_parser import parse_result_file
from .tables import result_table, write_result_csv


class PlaceholderPlotter:
    def __init__(self) -> None:
        self.interactor = QtWidgets.QLabel("PyVista rendering is disabled in this environment.")
        self.interactor.setAlignment(QtCore.Qt.AlignCenter)

    def clear(self) -> None:
        pass

    def set_background(self, color: str) -> None:
        self.interactor.setStyleSheet(f"background: {color}; color: #222;")

    def add_text(self, text: str, **kwargs) -> None:
        self.interactor.setText(text)

    def add_mesh(self, *args, **kwargs) -> None:
        pass

    def add_arrows(self, *args, **kwargs) -> None:
        pass

    def view_xy(self) -> None:
        pass

    def reset_camera(self) -> None:
        pass

    def close(self) -> None:
        pass


def _create_plotter(parent: QtWidgets.QWidget):
    if os.environ.get("QT_QPA_PLATFORM") == "offscreen":
        return PlaceholderPlotter()

    from pyvistaqt import QtInteractor

    return QtInteractor(parent)


class BemViewer(QtWidgets.QMainWindow):
    def __init__(self, input_path: Path | None = None, result_path: Path | None = None) -> None:
        super().__init__()
        self.setWindowTitle("MkBEM Result Viewer")
        self.resize(1280, 820)

        self.model: BemInput | None = None
        self.result: BemResult | None = None
        self.input_path: Path | None = None
        self.result_path: Path | None = None
        self.selected_boundary_index: int | None = None

        self.plotter = _create_plotter(self)
        self.summary = QtWidgets.QPlainTextEdit()
        self.summary.setReadOnly(True)
        self.summary.setMinimumWidth(360)
        self.boundary_table = QtWidgets.QTableWidget()
        self.boundary_table.setColumnCount(14)
        self.boundary_table.setHorizontalHeaderLabels([
            "node", "x", "y", "next", "mid x", "mid y",
            "ux code", "ux", "uy code", "uy", "T code", "T", "P code", "P",
        ])
        self.boundary_table.horizontalHeader().setStretchLastSection(True)
        self.boundary_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.boundary_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.boundary_table.setAlternatingRowColors(True)
        self.boundary_table.itemSelectionChanged.connect(self._boundary_selection_changed)
        self.result_table = QtWidgets.QTableWidget()
        self.result_table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.result_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.result_table.setAlternatingRowColors(True)

        self.time_combo = QtWidgets.QComboBox()
        self.time_combo.currentIndexChanged.connect(self._selection_changed)
        self.dataset_combo = QtWidgets.QComboBox()
        self.dataset_combo.addItems(["boundary", "interior"])
        self.dataset_combo.currentTextChanged.connect(self._dataset_changed)
        self.field_combo = QtWidgets.QComboBox()
        self.field_combo.currentTextChanged.connect(self._selection_changed)
        self.deformed_check = QtWidgets.QCheckBox("Show displacement vectors")
        self.deformed_check.setChecked(True)
        self.deformed_check.stateChanged.connect(self.refresh_scene)

        controls = QtWidgets.QWidget()
        controls_layout = QtWidgets.QVBoxLayout(controls)
        controls_layout.addWidget(self._button("Open Input", self.open_input))
        controls_layout.addWidget(self._button("Open Result", self.open_result))
        controls_layout.addWidget(self._button("Export Current Results CSV", self.export_current_results))
        controls_layout.addSpacing(12)
        controls_layout.addWidget(QtWidgets.QLabel("Time Step"))
        controls_layout.addWidget(self.time_combo)
        controls_layout.addWidget(QtWidgets.QLabel("Dataset"))
        controls_layout.addWidget(self.dataset_combo)
        controls_layout.addWidget(QtWidgets.QLabel("Field"))
        controls_layout.addWidget(self.field_combo)
        controls_layout.addWidget(self.deformed_check)
        controls_layout.addSpacing(12)
        details = QtWidgets.QTabWidget()
        details.addTab(self.summary, "Summary")
        details.addTab(self.boundary_table, "Boundary Coordinates")
        details.addTab(self.result_table, "Results")
        controls_layout.addWidget(details, stretch=1)

        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        splitter.addWidget(controls)
        splitter.addWidget(getattr(self.plotter, "interactor", self.plotter))
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        self.setCentralWidget(splitter)

        self._update_field_options()
        if input_path is not None:
            self.load_input(input_path)
        if result_path is not None:
            self.load_result(result_path)

    def closeEvent(self, event) -> None:  # noqa: N802
        self.plotter.close()
        super().closeEvent(event)

    def _button(self, text: str, slot) -> QtWidgets.QPushButton:
        button = QtWidgets.QPushButton(text)
        button.clicked.connect(slot)
        return button

    def open_input(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open MkBEM input", str(Path.cwd()), "MkBEM input (*.dat);;All files (*)")
        if path:
            self.load_input(Path(path))
            candidate = Path(path).with_suffix(".out")
            if candidate.exists():
                self.load_result(candidate)

    def open_result(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open MkBEM result", str(Path.cwd()), "MkBEM result (*.out);;All files (*)")
        if path:
            self.load_result(Path(path))

    def load_input(self, path: Path) -> None:
        try:
            self.model = parse_input_file(path)
            self.input_path = path
            self.selected_boundary_index = None
            self._refresh_summary()
            self._refresh_boundary_table()
            self.refresh_scene()
        except Exception as error:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "Input error", str(error))

    def load_result(self, path: Path) -> None:
        try:
            if self.model is None:
                self.result = parse_result_file(path)
            else:
                self.result = parse_result_file(path, boundary_count=self.model.boundary_count, interior_count=self.model.interior_point_count)
            self.result_path = path
            self._update_time_options()
            self._refresh_summary()
            self._refresh_result_table()
            self.refresh_scene()
        except Exception as error:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "Result error", str(error))

    def _update_time_options(self) -> None:
        self.time_combo.blockSignals(True)
        self.time_combo.clear()
        if self.result is not None:
            for index, step in enumerate(self.result.steps):
                self.time_combo.addItem(f"{index}: {step.time_year:g} yr", index)
        self.time_combo.blockSignals(False)

    def _selection_changed(self) -> None:
        self._refresh_result_table()
        self.refresh_scene()

    def _dataset_changed(self) -> None:
        self._update_field_options()
        self._refresh_result_table()

    def _update_field_options(self) -> None:
        self.field_combo.blockSignals(True)
        self.field_combo.clear()
        fields = BOUNDARY_FIELDS if self.dataset_combo.currentText() == "boundary" else INTERIOR_FIELDS
        self.field_combo.addItems(fields.keys())
        self.field_combo.blockSignals(False)
        self.refresh_scene()

    def _refresh_summary(self) -> None:
        if self.model is None:
            self.summary.setPlainText("Open a MkBEM .dat file to begin.")
            self.boundary_table.setRowCount(0)
            self.result_table.setRowCount(0)
            return

        lines = [
            f"Case: {self.model.title}",
            f"Input: {self.model.path}",
            f"Boundary nodes: {self.model.boundary_count}",
            f"Interior points: {self.model.interior_point_count}",
            f"Internal elements: {self.model.internal_element_count}",
            f"Time steps: {self.model.time_step_count}",
            f"Sources: {self.model.source_count}",
            "",
            "Material",
        ]
        for key, value in asdict(self.model.material).items():
            lines.append(f"  {key}: {value:g}")

        if self.result is not None:
            lines.extend([
                "",
                f"Result: {self.result.path}",
                f"Result time steps: {self.result.time_step_count}",
            ])
        self.summary.setPlainText("\n".join(lines))

    def _refresh_boundary_table(self) -> None:
        if self.model is None:
            self.boundary_table.setRowCount(0)
            return

        self.boundary_table.blockSignals(True)
        self.boundary_table.setRowCount(self.model.boundary_count)
        for index, (x, y) in enumerate(self.model.boundary_nodes):
            next_index = 0 if index == self.model.boundary_count - 1 else index + 1
            next_x, next_y = self.model.boundary_nodes[next_index]
            bc = self.model.boundary_conditions[index]
            values = [
                str(index),
                f"{x:.6g}",
                f"{y:.6g}",
                str(next_index),
                f"{(x + next_x) / 2.0:.6g}",
                f"{(y + next_y) / 2.0:.6g}",
                str(bc.ux_code),
                f"{bc.ux_value:.6g}",
                str(bc.uy_code),
                f"{bc.uy_value:.6g}",
                str(bc.temp_code),
                f"{bc.temp_value:.6g}",
                str(bc.pressure_code),
                f"{bc.pressure_value:.6g}",
            ]
            for column, value in enumerate(values):
                self.boundary_table.setItem(index, column, QtWidgets.QTableWidgetItem(value))

        self.boundary_table.resizeColumnsToContents()
        self.boundary_table.blockSignals(False)

    def _current_step(self):
        if self.result is None or not self.result.steps:
            return None
        index = self.time_combo.currentIndex()
        if index < 0:
            index = 0
        return self.result.steps[index]

    def _refresh_result_table(self) -> None:
        step = self._current_step()
        if step is None:
            self.result_table.setRowCount(0)
            self.result_table.setColumnCount(0)
            return

        headers, rows = result_table(self.dataset_combo.currentText(), step)
        self.result_table.setColumnCount(len(headers))
        self.result_table.setHorizontalHeaderLabels(headers)
        self.result_table.setRowCount(len(rows))
        for row_index, row in enumerate(rows):
            for column_index, value in enumerate(row):
                if isinstance(value, float):
                    text = f"{value:.8g}"
                else:
                    text = str(value)
                self.result_table.setItem(row_index, column_index, QtWidgets.QTableWidgetItem(text))
        self.result_table.resizeColumnsToContents()

    def export_current_results(self) -> None:
        step = self._current_step()
        if step is None:
            QtWidgets.QMessageBox.information(self, "Export results", "Open a result file before exporting.")
            return

        dataset = self.dataset_combo.currentText()
        default_name = f"mkbem_{dataset}_step_{self.time_combo.currentIndex():03d}.csv"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export current results", default_name, "CSV files (*.csv);;All files (*)")
        if not path:
            return
        try:
            write_result_csv(path, dataset, step)
        except Exception as error:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "Export error", str(error))

    def _boundary_selection_changed(self) -> None:
        selected = self.boundary_table.selectionModel().selectedRows()
        self.selected_boundary_index = selected[0].row() if selected else None
        self.refresh_scene()

    def refresh_scene(self) -> None:
        self.plotter.clear()
        self.plotter.set_background("white")

        if self.model is None:
            self.plotter.add_text("Open a MkBEM input file", color="black", font_size=12)
            return

        boundary = make_boundary_polyline(self.model)
        self.plotter.add_mesh(boundary, color="black", line_width=3, label="boundary")
        if self.selected_boundary_index is not None:
            selected = make_boundary_segment(self.model, self.selected_boundary_index)
            self.plotter.add_mesh(selected, color="#d55e00", line_width=8, label="selected boundary element")

        if self.model.interior_points:
            interior = make_interior_points(self.model)
            self.plotter.add_mesh(interior, color="#666666", point_size=9, render_points_as_spheres=True, label="interior points")

        if self.result is not None and self.result.steps:
            step = self._current_step()
            if step is None:
                return
            field = self.field_combo.currentText()
            if self.dataset_combo.currentText() == "boundary" and step.boundary:
                mesh = make_boundary_results(step, field)
                self.plotter.add_mesh(mesh, scalars=field, point_size=14, render_points_as_spheres=True, cmap="viridis", scalar_bar_args={"title": field})
                if self.deformed_check.isChecked() and "displacement" in mesh.point_data:
                    self.plotter.add_arrows(mesh.points, mesh["displacement"], mag=0.35, color="#d55e00")
            elif step.interior_primary:
                mesh = make_interior_results(step, field)
                self.plotter.add_mesh(mesh, scalars=field, point_size=15, render_points_as_spheres=True, cmap="viridis", scalar_bar_args={"title": field})
                if self.deformed_check.isChecked() and "displacement" in mesh.point_data:
                    self.plotter.add_arrows(mesh.points, mesh["displacement"], mag=0.35, color="#d55e00")

        self.plotter.view_xy()
        self.plotter.reset_camera()


def main() -> int:
    parser = argparse.ArgumentParser(description="Open a PyQt/PyVista viewer for MkBEM files.")
    parser.add_argument("--input", type=Path, help="MkBEM .dat input file")
    parser.add_argument("--result", type=Path, help="MkBEM .out result file")
    args = parser.parse_args()

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    viewer = BemViewer(args.input, args.result)
    viewer.show()
    return app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())