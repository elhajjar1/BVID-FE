"""Small input-validation helpers shared by the BVID-FE dataclasses.

``NaN <= 0`` and ``inf > 0`` are both ``False``/``True`` in the wrong way for
a bare ``x <= 0`` guard, so every numeric field goes through
:func:`require_finite` (rejects ``bool``, strings, NaN and +/-inf) before its
range is checked. All failures raise :class:`ValueError` naming the field.
"""

from __future__ import annotations

import dataclasses
import math
import numbers
from collections.abc import Mapping
from typing import Any


def require_finite(value: Any, name: str) -> float:
    """Return *value* as a float if it is a finite real number; else raise.

    ``bool`` is rejected even though it is an ``int`` subclass, as are
    strings, so ``"0.152"`` or ``True`` cannot slip through as a dimension.
    """
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ValueError(f"{name} must be a number (got {type(value).__name__}: {value!r})")
    f = float(value)
    if not math.isfinite(f):
        raise ValueError(f"{name} must be finite (got {value!r})")
    return f


def require_positive(value: Any, name: str) -> float:
    """Finite and ``> 0``."""
    f = require_finite(value, name)
    if f <= 0:
        raise ValueError(f"{name} must be > 0 (got {value!r})")
    return f


def require_non_negative(value: Any, name: str) -> float:
    """Finite and ``>= 0``."""
    f = require_finite(value, name)
    if f < 0:
        raise ValueError(f"{name} must be >= 0 (got {value!r})")
    return f


def require_xy(value: Any, name: str) -> tuple:
    """A pair of finite numbers, returned as a tuple (lists from JSON are
    accepted). The element objects are kept as given."""
    if isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an (x, y) pair of numbers (got {value!r})")
    try:
        pair = tuple(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an (x, y) pair of numbers (got {value!r})") from exc
    if len(pair) != 2:
        raise ValueError(f"{name} must be an (x, y) pair of numbers (got {value!r})")
    require_finite(pair[0], f"{name}[0]")
    require_finite(pair[1], f"{name}[1]")
    return pair


def dataclass_kwargs(cls: Any, data: Any, name: str) -> dict:
    """Check that *data* is a mapping of keyword arguments for dataclass
    *cls* and return it as a ``dict``, for the ``from_dict`` loaders.

    Raises ``ValueError`` naming *name* for a non-mapping, a key that is
    not a field of *cls* (a typo would otherwise fall back to the default
    silently), or a missing field that has no default.
    """
    if not isinstance(data, Mapping):
        raise ValueError(f"{name} must be an object/dict (got {type(data).__name__}: {data!r})")
    init_fields = [f for f in dataclasses.fields(cls) if f.init]
    names = [f.name for f in init_fields]
    unknown = [k for k in data if k not in names]
    if unknown:
        raise ValueError(f"{name} has unknown field(s) {unknown}; {cls.__name__} fields: {names}")
    missing = [
        f.name
        for f in init_fields
        if f.name not in data
        and f.default is dataclasses.MISSING
        and f.default_factory is dataclasses.MISSING
    ]
    if missing:
        raise ValueError(f"{name} is missing required field(s) {missing}")
    return dict(data)
