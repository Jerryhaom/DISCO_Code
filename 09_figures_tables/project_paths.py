"""Portable paths shared by figure utilities (standard library only)."""
import os
from pathlib import Path
CODE_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get("DISCO_DATA_ROOT", CODE_ROOT / "data")).expanduser().resolve()
OUTPUT_ROOT = Path(os.environ.get("DISCO_OUTPUT_ROOT", CODE_ROOT / "results")).expanduser().resolve()
