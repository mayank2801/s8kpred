"""
s8kpred/config.py
-----------------
Centralised configuration: all hard-coded paths live here so every other
module can import from a single source of truth.
"""
import os
from pathlib import Path

# ------------------------------------------------------------------
# Package data directory (ships with the wheel / sdist)
# ------------------------------------------------------------------
_PKG_DIR = Path(__file__).parent
DATA_DIR = _PKG_DIR / "data"

# Built-in model & lookup-table paths
TRIPEPTIDE_PROPENSITY_3STATE = DATA_DIR / "TriPeptidePropensityThreeStateSecStructure2AND.csv"
TRIPEPTIDE_PROPENSITY_8STATE = DATA_DIR / "TriPeptidePropensityEightStateSecStructure.csv"
TRIPEPTIDE_BINARY_TABLE      = DATA_DIR / "TripeptideBinaryTable_60.csv"
MODEL_3STATE                 = DATA_DIR / "model_3state.json"
MODEL_8STATE                 = DATA_DIR / "model_8state.ubj"

# ------------------------------------------------------------------
# External tools — overridable via environment variables or CLI flags
# ------------------------------------------------------------------
PSIBLAST_DEFAULT = os.environ.get(
    "S8KPRED_PSIBLAST",
    "psiblast"          # assumes psiblast is on PATH by default
)
BLASTDB_DEFAULT = os.environ.get(
    "S8KPRED_BLASTDB",
    ""                  # must be supplied by the user
)

PSIBLAST_ITERATIONS = int(os.environ.get("S8KPRED_ITERATIONS", "3"))

# ------------------------------------------------------------------
# Secondary structure mappings
# ------------------------------------------------------------------
IDX2CHAR_3STATE = {0: "E", 1: "H", 2: "L"}
IDX2CHAR_8STATE = {0: "B", 1: "E", 2: "G", 3: "H", 4: "I", 5: "L", 6: "S", 7: "T"}

# Map 8-state labels → simplified 3-state (H/E/L) for cartoon plots
MAP_8TO3 = {"H": "H", "G": "H", "I": "H",
            "E": "E", "B": "E",
            "T": "L", "S": "L", "L": "L"}

WINDOW_SIZE      = 17
PADDING_RESIDUE  = "G"
PADDING_LEN      = 8
