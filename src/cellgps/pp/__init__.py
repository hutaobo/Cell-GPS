"""Preprocessing API following the scverse ``pp`` namespace convention."""

from __future__ import annotations

from typing import Any

import cellgps.preprocessing as _preprocessing

__all__ = list(_preprocessing.__all__)


def __getattr__(name: str) -> Any:
    value = getattr(_preprocessing, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
