"""
s8kpred/features/propensity.py
-------------------------------
Build tripeptide-propensity feature matrices for 3-state and 8-state models.
All logic extracted from SecondaryStructurePredictionThreeState.py and
SecondaryStructurePredictionEightState.py, made reusable.
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

from s8kpred.config import (
    TRIPEPTIDE_PROPENSITY_3STATE,
    TRIPEPTIDE_PROPENSITY_8STATE,
    TRIPEPTIDE_BINARY_TABLE,
    WINDOW_SIZE,
    PADDING_LEN,
    PADDING_RESIDUE,
)
from s8kpred.utils.fasta import SeqRecord


def _pad(seq: str) -> str:
    pad = PADDING_RESIDUE * PADDING_LEN
    return (pad + seq.upper() + pad)


def _safe_tripeptide(seq: str, k: int) -> str:
    """Return GXG-safe tripeptide centred at position k."""
    n   = len(seq)
    s1  = seq[k - 1] if k > 0 else PADDING_RESIDUE
    s2  = seq[k]
    s3  = seq[k + 1] if k + 1 < n else PADDING_RESIDUE
    return (s1 + s2 + s3).replace("X", PADDING_RESIDUE)


def _lookup_propensity(propdf: pd.DataFrame, tripeptide: str, n_cols: int) -> List[float]:
    """Return propensity values for a tripeptide; fall back to zeros."""
    idx = propdf.index[propdf["TriPeptide"] == tripeptide].tolist()
    if not idx:
        return [0.0] * n_cols
    return propdf.iloc[idx[0], 2: 2 + n_cols].tolist()


# ── 3-state feature builder ────────────────────────────────────────────────

def build_3state_features(
    records: List[SeqRecord],
    propensity_csv: Path = TRIPEPTIDE_PROPENSITY_3STATE,
    binary_csv: Path = TRIPEPTIDE_BINARY_TABLE,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build the propensity + PSSM input DataFrame for the 3-state model.

    Returns
    -------
    inputdf  : propensity feature rows (one per window)
    binarydf : the full binary lookup table (for later merge)
    """
    prop_df   = pd.read_csv(propensity_csv)
    binary_df = pd.read_csv(binary_csv)

    PROP_COLS = 3  # H, E, L propensities per residue

    col_names = ["ID", "sequence", "Residue9th"]
    for r in range(1, WINDOW_SIZE + 1):
        for s in ("H", "E", "L"):
            col_names.append(f"{s}{r}")

    rows = []
    for seq_idx, record in enumerate(records, start=1):
        seq_id  = f"Seq_{seq_idx}"
        padded  = _pad(record.seq)
        seq_len = len(padded)

        for j in range(seq_len - WINDOW_SIZE + 1):
            window = padded[j: j + WINDOW_SIZE]
            if len(window) != WINDOW_SIZE:
                continue

            residue9 = padded[j + 8]
            propensities: List[float] = []

            for k in range(j, j + WINDOW_SIZE):
                tri   = _safe_tripeptide(padded, k)
                propensities.extend(_lookup_propensity(prop_df, tri, PROP_COLS))

            row = [seq_id, window, residue9] + propensities
            rows.append(row)

    inputdf = pd.DataFrame(rows, columns=col_names)
    return inputdf, binary_df


# ── 8-state feature builder ────────────────────────────────────────────────

def build_8state_features(
    records: List[SeqRecord],
    propensity_csv: Path = TRIPEPTIDE_PROPENSITY_8STATE,
    binary_csv: Path = TRIPEPTIDE_BINARY_TABLE,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build the propensity + binary feature DataFrame for the 8-state model.

    Returns
    -------
    inputdf  : propensity + binary feature rows (one per window)
    binarydf : the full binary lookup table (for reference)
    """
    prop_df   = pd.read_csv(propensity_csv)
    binary_df = pd.read_csv(binary_csv)

    PROP_COLS = 8  # B, E, G, H, I, L, S, T propensities per residue

    prop_col_names = []
    for r in range(1, WINDOW_SIZE + 1):
        for s in ("B", "E", "G", "H", "I", "L", "S", "T"):
            prop_col_names.append(f"PR{r}_{s}")

    binary_col_names = [f"R{r}_{aa}" for r in range(1, 4) for aa in "ACDEFGHIKLMNPQRSTVWY"]

    col_names = ["ID", "sequence", "Residue9th"] + prop_col_names + binary_col_names

    rows = []
    for seq_idx, record in enumerate(records, start=1):
        seq_id  = f"Seq_{seq_idx}"
        padded  = _pad(record.seq)
        seq_len = len(padded)

        for j in range(seq_len - WINDOW_SIZE + 1):
            window = padded[j: j + WINDOW_SIZE]
            if len(window) != WINDOW_SIZE:
                continue

            residue9 = padded[j + 8]
            propensities: List[float] = []

            for k in range(j, j + WINDOW_SIZE):
                tri = _safe_tripeptide(padded, k)
                propensities.extend(_lookup_propensity(prop_df, tri, PROP_COLS))

            # Central tripeptide binary encoding (positions 7-9)
            mid_tri  = (padded[j + 7] + padded[j + 8] + padded[j + 9]).replace("X", PADDING_RESIDUE)
            bin_idx  = binary_df.index[binary_df["Tripeptide"] == mid_tri].tolist()
            if bin_idx:
                bin_vals = binary_df.iloc[bin_idx[0], 1:61].tolist()
            else:
                bin_vals = [0] * 60

            row = [seq_id, window, residue9] + propensities + bin_vals
            rows.append(row)

    inputdf = pd.DataFrame(rows, columns=col_names)
    return inputdf, binary_df
