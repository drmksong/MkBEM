"""Utilities for reading and inspecting MkBEM benchmark data."""

from .input_parser import parse_input_file
from .result_parser import parse_result_file

__all__ = ["parse_input_file", "parse_result_file"]
