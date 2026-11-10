from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from .models import BemInput, BoundaryCondition, InternalElement, MaterialProperties, Source


class TokenReader:
    def __init__(self, tokens: list[str]) -> None:
        self._tokens = tokens
        self._index = 0

    def _next(self) -> str:
        if self._index >= len(self._tokens):
            raise ValueError("Unexpected end of MkBEM input file")
        token = self._tokens[self._index]
        self._index += 1
        return token

    def integer(self, name: str) -> int:
        token = self._next()
        try:
            return int(token)
        except ValueError as error:
            raise ValueError(f"Expected integer for {name}, got {token!r}") from error

    def real(self, name: str) -> float:
        token = self._next()
        try:
            return float(token)
        except ValueError as error:
            raise ValueError(f"Expected number for {name}, got {token!r}") from error

    def word(self, name: str) -> str:
        return self._next()


def _strip_comment(line: str) -> str:
    return line.split(";", 1)[0]


def _tokens(lines: Iterator[str]) -> list[str]:
    return "\n".join(_strip_comment(line) for line in lines).split()


def parse_input_file(path: str | Path) -> BemInput:
    input_path = Path(path)
    with input_path.open("r", encoding="utf-8", errors="replace") as file:
        title = file.readline().strip()
        reader = TokenReader(_tokens(file))

    boundary_count = reader.integer("N")
    interior_point_count = reader.integer("L")
    boundary_group_count = reader.integer("M")
    source_count = reader.integer("NS")

    boundary_group_last_nodes = [reader.integer(f"Last[{index}]") for index in range(boundary_group_count)]

    internal_node_count = reader.integer("NI")
    internal_element_count = reader.integer("NE")
    time_step_count = reader.integer("NT")

    material = MaterialProperties(
        mu=reader.real("Mu"),
        nu=reader.real("Nu"),
        k_t=reader.real("K_T"),
        k_p=reader.real("K_P"),
        beta_ts=reader.real("Beta_TS"),
        alpha_p=reader.real("Alpha_P"),
        initial_temperature=reader.real("T_0"),
        eta=reader.real("Eta"),
        rho=reader.real("Rho"),
        cp=reader.real("Cp"),
        bulk_modulus_input=reader.real("K"),
        skempton_b=reader.real("B"),
    )

    boundary_nodes = [(reader.real(f"X[{index}]"), reader.real(f"Y[{index}]")) for index in range(boundary_count)]

    boundary_conditions = []
    for index in range(boundary_count):
        boundary_conditions.append(
            BoundaryCondition(
                ux_code=reader.integer(f"Code[{index}, ux]"),
                ux_value=reader.real(f"Bc[{index}, ux]"),
                uy_code=reader.integer(f"Code[{index}, uy]"),
                uy_value=reader.real(f"Bc[{index}, uy]"),
                temp_code=reader.integer(f"Code[{index}, temp]"),
                temp_value=reader.real(f"Bc[{index}, temp]"),
                pressure_code=reader.integer(f"Code[{index}, pressure]"),
                pressure_value=reader.real(f"Bc[{index}, pressure]"),
            )
        )

    interior_points = [(reader.real(f"Xi[{index}]"), reader.real(f"Yi[{index}]")) for index in range(interior_point_count)]
    time_steps = [reader.real(f"Time[{index}]") for index in range(time_step_count)]

    sources = []
    for index in range(source_count):
        sources.append(
            Source(
                kind=reader.word(f"SourceType[{index}]").upper(),
                strength=reader.real(f"Source[{index}]") ,
                x=reader.real(f"SourceX[{index}]"),
                y=reader.real(f"SourceY[{index}]"),
            )
        )

    internal_nodes = [(reader.real(f"NodeX[{index}]"), reader.real(f"NodeY[{index}]")) for index in range(internal_node_count)]
    internal_elements = []
    for index in range(internal_element_count):
        internal_elements.append(
            InternalElement(
                nodes=(
                    reader.integer(f"Elem[{index}, 0]"),
                    reader.integer(f"Elem[{index}, 1]"),
                    reader.integer(f"Elem[{index}, 2]"),
                    reader.integer(f"Elem[{index}, 3]"),
                ),
                original_volumetric_strain=reader.real(f"OriVolStrain[{index}]"),
                original_temperature=reader.real(f"OriTemp[{index}]"),
                original_pressure=reader.real(f"OriPress[{index}]"),
            )
        )

    return BemInput(
        path=input_path,
        title=title,
        boundary_count=boundary_count,
        interior_point_count=interior_point_count,
        boundary_group_count=boundary_group_count,
        source_count=source_count,
        internal_node_count=internal_node_count,
        internal_element_count=internal_element_count,
        time_step_count=time_step_count,
        material=material,
        boundary_nodes=boundary_nodes,
        boundary_conditions=boundary_conditions,
        interior_points=interior_points,
        time_steps=time_steps,
        sources=sources,
        boundary_group_last_nodes=boundary_group_last_nodes,
        internal_nodes=internal_nodes,
        internal_elements=internal_elements,
    )
