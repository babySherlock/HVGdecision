"""Command-line entry point for HVGDecision 0.9."""

from __future__ import annotations

import argparse
import json

from . import __version__
from .config import load_yaml
from .inspect import inspect_dataset
from .modes import normalize_mode
from .refine import refine_panel


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        prog="hvgdecision",
        description=(
            "Integration-aware HVG refinement. Select within_domain or "
            "cross_domain explicitly."
        ),
    )
    command.add_argument("--version", action="version", version=__version__)
    subcommands = command.add_subparsers(dest="command", required=True)

    inspect_cmd = subcommands.add_parser("inspect")
    inspect_cmd.add_argument("--config", required=True)

    refine_cmd = subcommands.add_parser("refine")
    refine_cmd.add_argument("--config", required=True)
    refine_cmd.add_argument(
        "--mode",
        choices=["within_domain", "cross_domain"],
        help="Algorithm mode. May also be set as mode: in YAML.",
    )
    refine_cmd.add_argument(
        "--delete-budget",
        type=int,
        default=None,
        help="Cross-domain fixed removal budget (commonly 5, 10, or 20).",
    )
    return command


def main() -> int:
    args = parser().parse_args()
    config = load_yaml(args.config)
    if args.command == "inspect":
        result = inspect_dataset(config)
    else:
        if args.mode is not None:
            config["mode"] = args.mode
        config["mode"] = normalize_mode(config.get("mode"))
        if args.delete_budget is not None:
            config["cross_domain_delete_budget"] = int(args.delete_budget)
        result = refine_panel(config)
    print(json.dumps(result, ensure_ascii=False, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
