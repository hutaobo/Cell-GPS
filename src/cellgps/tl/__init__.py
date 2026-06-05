"""Analysis API following the scverse ``tl`` namespace convention."""

from __future__ import annotations

from typing import Any

import cellgps.analysis as _analysis

__all__ = list(_analysis.__all__)


def __getattr__(name: str) -> Any:
    value = getattr(_analysis, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
