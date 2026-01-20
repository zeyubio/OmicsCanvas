#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""OmicsCanvas CLI (pip package entrypoint)

Design goals:
- Do NOT modify any step script logic.
- Provide a stable command `omicscanvas` to list and run numbered step scripts.

Examples:
  omicscanvas list
  omicscanvas run 10 -- --help
  omicscanvas run 10 -- --mode 3d --gene-type gene --dir-chip caculate_matrix

Notes:
- Arguments after `--` are passed to the step script unchanged.
"""

from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path


def _steps_dir() -> Path:
    # Installed package path
    return (Path(__file__).resolve().parent / "steps").resolve()


def list_steps() -> list[Path]:
    sd = _steps_dir()
    return sorted(sd.glob("[0-9][0-9]_*.py"))


def _find_step(step: str) -> Path:
    sd = _steps_dir()
    step = step.strip()
    if len(step) == 1:
        step = f"0{step}"
    matches = sorted(sd.glob(f"{step}_*.py"))
    if not matches:
        raise FileNotFoundError(f"Step {step} not found under: {sd}")
    if len(matches) > 1:
        raise RuntimeError(
            "Step matched multiple scripts: " + ", ".join([m.name for m in matches])
        )
    return matches[0]


def cmd_list(_: argparse.Namespace) -> int:
    for fp in list_steps():
        print(fp.name)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    target = _find_step(args.step)

    passthrough = args.passthrough
    if passthrough and passthrough[0] == "--":
        passthrough = passthrough[1:]

    # Replace argv for the executed script
    sys.argv = [str(target)] + passthrough
    runpy.run_path(str(target), run_name="__main__")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="omicscanvas",
        description="Run OmicsCanvas step scripts via a stable CLI.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    sub = p.add_subparsers(dest="cmd", required=True)

    p_list = sub.add_parser("list", help="List available step scripts.")
    p_list.set_defaults(func=cmd_list)

    p_run = sub.add_parser("run", help="Run a step script by number (01..16).")
    p_run.add_argument("step", help="Step number, e.g. 10 or 01")
    p_run.add_argument(
        "passthrough",
        nargs=argparse.REMAINDER,
        help="Arguments passed to the step script (use `--` to separate).",
    )
    p_run.set_defaults(func=cmd_run)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    ns = parser.parse_args(argv)
    return int(ns.func(ns))


if __name__ == "__main__":
    raise SystemExit(main())
