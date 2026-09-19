"""Configuration loading and validation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_yaml(path: str | Path) -> dict[str, Any]:
    path = Path(path).expanduser().resolve()
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise TypeError(f"Configuration must be a mapping: {path}")
    data["_config_path"] = str(path)
    return data


def require(config: dict[str, Any], *keys: str) -> None:
    missing = [key for key in keys if key not in config or config[key] in (None, "")]
    if missing:
        raise KeyError(f"Configuration is missing required keys: {missing}")


def resolve_output(config: dict[str, Any], override: str | None = None) -> Path:
    value = override or config.get("output_dir")
    if not value:
        raise KeyError("output_dir is required")
    path = Path(value).expanduser().resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path
