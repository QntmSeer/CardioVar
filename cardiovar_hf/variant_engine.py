"""Hugging Face Space variant_engine — delegates to real Enformer inference with no synthetic fallback."""

import sys
from pathlib import Path

_ROOT = str(Path(__file__).resolve().parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from variant_engine import compute_variant_impact  # noqa: E402, F401
