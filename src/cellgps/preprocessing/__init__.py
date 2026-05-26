"""Cell-GPS preprocessing namespace."""

from __future__ import annotations

from typing import Any

import sfplot.preprocessing as _legacy

__path__ = list(_legacy.__path__)
__all__ = list(_legacy.__all__)


def __getattr__(name: str) -> Any:
    value = getattr(_legacy, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
