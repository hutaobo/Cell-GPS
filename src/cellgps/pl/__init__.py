"""Plotting API following the scverse ``pl`` namespace convention."""

from __future__ import annotations

from typing import Any

import cellgps.plotting as _plotting

__all__ = list(_plotting.__all__)


def __getattr__(name: str) -> Any:
    value = getattr(_plotting, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
