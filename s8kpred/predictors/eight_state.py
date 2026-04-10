"""
s8kpred/predictors/eight_state.py
----------------------------------
8-state (H / G / I / E / B / T / S / L) secondary structure predictor.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from s8kpred.config import IDX2CHAR_8STATE, MODEL_8STATE
from s8kpred.features.propensity import build_8state_features
from s8kpred.utils.fasta import SeqRecord
from s8kpred.utils.io import write_ss2, write_horiz, write_fasta_result


def predict_8state(
    records: List[SeqRecord],
    pssm_csv: Path,
    job_dir: Path,
    model_path: Path = MODEL_8STATE,
    verbose: bool = True,
) -> Dict[str, str]:
    """
    Run 8-state secondary structure prediction.

    Parameters
    ----------
    records    : list of SeqRecord (from FASTA parser)
    pssm_csv   : path to the PSSM feature CSV produced by extract_pssm_features()
    job_dir    : directory where result files will be written
    model_path : path to the XGBoost .ubj model file
    verbose    : print per-sequence summaries

    Returns
    -------
    dict mapping sequence ID → predicted secondary structure string
    """
    import xgboost as xgb

    # ── Feature construction ─────────────────────────────────────────────
    inputdf, _ = build_8state_features(records)
    pssm_df    = pd.read_csv(pssm_csv)

    merged = pd.merge(inputdf, pssm_df, on=["ID", "sequence"], how="inner")

    if len(merged) < 2:
        raise ValueError("Not enough merged feature rows for 8-state prediction.")

    # ── Load model & predict ─────────────────────────────────────────────
    model = xgb.XGBClassifier()
    model.load_model(str(model_path))

    X             = merged.iloc[:, 3:]
    y_idx         = model.predict(X)
    probabilities = model.predict_proba(X).astype(float)

    vec_map = np.vectorize(lambda x: IDX2CHAR_8STATE.get(int(x), "L"))
    y_pred  = vec_map(y_idx)

    del model

    # ── Build per-sequence results ───────────────────────────────────────
    results: Dict[str, str] = {}
    seq_ids = merged["ID"].unique()

    for seq_id in seq_ids:
        mask      = merged["ID"] == seq_id
        residues  = merged.loc[mask, "Residue9th"].tolist()
        pred_ss   = y_pred[mask]
        probs     = probabilities[mask]
        record_id = records[int(seq_id.split("_")[1]) - 1].id

        predicted_str = "".join(pred_ss)
        results[record_id] = predicted_str

        if verbose:
            print(f"\n  [8-state] {record_id}")
            print(f"    Sequence : {''.join(residues)}")
            print(f"    Predicted: {predicted_str}")

        probs_arr = np.array(probs)

        write_ss2(
            job_dir / "ResultEightState.ss2",
            residues, pred_ss, probs_arr,
            header=record_id, n_classes=8,
        )
        write_horiz(
            job_dir / "ResultEightState.horiz",
            residues, pred_ss, probs_arr,
            header=record_id,
        )
        write_fasta_result(
            job_dir / "ResultEightState.fas",
            residues, pred_ss,
            header=record_id,
        )

        df_out = pd.DataFrame({
            "Residue":             residues,
            "Secondary Structure": pred_ss,
        })
        prob_cols = pd.DataFrame(probs_arr, columns=["B", "E", "G", "H", "I", "L", "S", "T"])
        df_out = pd.concat([df_out, prob_cols], axis=1)
        df_out.to_csv(job_dir / "ResultEightState.csv", mode="a", index=True)

    return results
