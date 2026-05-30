"""Public Python package for Cell-GPS.

The original implementation remains under ``sfplot`` for backward
compatibility; new code should import from ``cellgps``.
"""

from __future__ import annotations

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from typing import Any

import sfplot as _legacy

__author__ = _legacy.__author__
__email__ = _legacy.__email__

try:
    __version__ = version("Cell-GPS")
except PackageNotFoundError:
    __version__ = getattr(_legacy, "__version__", "0.0.5")

_SUBPACKAGES = {"analysis", "preprocessing", "plotting", "gui"}

__all__ = sorted(set(_legacy.__all__) | _SUBPACKAGES)


def __getattr__(name: str) -> Any:
    if name in _SUBPACKAGES:
        module = import_module(f".{name}", __name__)
        globals()[name] = module
        return module

    value = getattr(_legacy, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
