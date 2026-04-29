from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from isograph_benchmark.paths import rel


def load_yaml(path: str | Path) -> dict[str, Any]:
    p = Path(path)
    if not p.is_absolute():
        p = rel(str(p))
    return yaml.safe_load(p.read_text()) or {}

