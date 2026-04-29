from __future__ import annotations

from pathlib import Path

try:
    from pyhere import here as _pyhere
except ImportError:  # pragma: no cover - local bootstrap before env creation
    _pyhere = None


def root() -> Path:
    if _pyhere is not None:
        return Path(_pyhere()).resolve()
    current = Path.cwd().resolve()
    for candidate in [current, *current.parents]:
        if (candidate / ".here").exists():
            return candidate
    return current


def rel(*parts: str) -> Path:
    return root().joinpath(*parts)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path
