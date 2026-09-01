"""Shared pytest fixtures for the sedumipy test suite.

Currently just the benchmark-metrics collector used by
test_benchmarks.py (SDPLIB/DIMACS problems, see that file's own
docstring) -- a session-scoped list the benchmark tests append per-
problem timing/accuracy results to, printed as a summary table and
written to benchmark_results.csv once the whole run finishes.
"""

from __future__ import annotations

import csv
import pathlib
from typing import Any

import pytest

_METRICS_FIELDS = ["source", "name", "time_s", "iter", "pobj", "numerr", "status"]
_CSV_PATH = pathlib.Path(__file__).parent.parent / "benchmark_results.csv"

# Module-level store shared between the fixture and the hook (a
# pytest_terminal_summary hook has no access to fixtures).
_benchmark_results: list[dict[str, Any]] = []


@pytest.fixture(scope="session")
def benchmark_collector() -> list[dict[str, Any]]:
    """Session-scoped list that benchmark tests append result dicts to.

    Printed in pytest_terminal_summary rather than via plain stdout so
    it survives pytest's output capturing regardless of -s/--capture.
    """
    return _benchmark_results


def pytest_terminal_summary(terminalreporter, exitstatus, config) -> None:  # noqa: ARG001
    results = _benchmark_results
    if not results:
        return

    results.sort(key=lambda r: (r["source"], r["name"]))

    header = f"{'Source':<7} {'Name':<20} {'Time(s)':>8} {'Iter':>5} {'pobj':>16} {'numerr':>7} {'status'}"
    sep = "-" * len(header)
    terminalreporter.write_sep("=", "BENCHMARK METRICS SUMMARY")
    terminalreporter.write_line(header)
    terminalreporter.write_line(sep)
    prev_source = None
    for r in results:
        if prev_source and r["source"] != prev_source:
            terminalreporter.write_line(sep)
        prev_source = r["source"]
        terminalreporter.write_line(
            f"{r['source']:<7} {r['name']:<20} {r['time_s']:>8.3f} {r['iter']:>5d} "
            f"{r['pobj']:>16.6g} {r['numerr']:>7d} {r['status']}"
        )
    terminalreporter.write_sep("=", f"{len(results)} problems")

    with open(_CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_METRICS_FIELDS)
        writer.writeheader()
        writer.writerows(results)
    terminalreporter.write_line(f"Benchmark results written to: {_CSV_PATH}")
