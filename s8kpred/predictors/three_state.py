"""
s8kpred/predictors/three_state.py
----------------------------------
3-state (H / E / L) secondary structure predictor.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from s8kpred.config import IDX2CHAR_3STATE, MODEL_3STATE
from s8kpred.features.propensity import build_3state_features
from s8kpred.utils.fasta import SeqRecord
from s8kpred.utils.io import write_ss2, write_horiz, write_fasta_result


def predict_3state(
    records: List[SeqRecord],
    pssm_csv: Path,
    job_dir: Path,
    model_path: Path = MODEL_3STATE,
    verbose: bool = True,
) -> Dict[str, str]:
    """
    Run 3-state secondary structure prediction.

    Parameters
    ----------
    records    : list of SeqRecord (from FASTA parser)
    pssm_csv   : path to the PSSM feature CSV produced by extract_pssm_features()
    job_dir    : directory where result files will be written
    model_path : path to the XGBoost .json model file
    verbose    : print per-sequence summaries

    Returns
    -------
    dict mapping sequence ID → predicted secondary structure string
    """
    import xgboost as xgb

    # ── Feature construction ─────────────────────────────────────────────
    inputdf, binary_df = build_3state_features(records)
    pssm_df  = pd.read_csv(pssm_csv)

    merged = pd.merge(inputdf, pssm_df, on=["ID", "sequence"], how="inner")

    # Extract central tripeptide for binary merge
    merged["Tripeptide"] = (
        merged["sequence"].astype(str).str.replace("X", "G").str[7:10]
    )
    merged = pd.merge(merged, binary_df, how="left", on="Tripeptide")
    merged = merged.drop(columns=["Tripeptide"])
    merged = merged.dropna()

    if len(merged) < 2:
        raise ValueError("Not enough merged feature rows for 3-state prediction.")

    # ── Load model & predict ─────────────────────────────────────────────
    model = xgb.XGBClassifier()
    model.load_model(str(model_path))

    X            = merged.iloc[:, 3:]
    y_idx        = model.predict(X)
    probabilities = model.predict_proba(X).astype(float)

    vec_map = np.vectorize(lambda x: IDX2CHAR_3STATE.get(int(x), "L"))
    y_pred  = vec_map(y_idx)

    del model  # free memory

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
            print(f"\n  [3-state] {record_id}")
            print(f"    Sequence : {''.join(residues)}")
            print(f"    Predicted: {predicted_str}")

        # ── Write output files ───────────────────────────────────────────
        probs_arr = np.array(probs)

        write_ss2(
            job_dir / "ResultThreeState.ss2",
            residues, pred_ss, probs_arr,
            header=record_id, n_classes=3,
        )
        write_horiz(
            job_dir / "ResultThreeState.horiz",
            residues, pred_ss, probs_arr,
            header=record_id,
        )
        write_fasta_result(
            job_dir / "ResultThreeState.fas",
            residues, pred_ss,
            header=record_id,
        )

        # DataFrame CSV
        df_out = pd.DataFrame({
            "Residue":            residues,
            "Secondary Structure": pred_ss,
            "ID":                 seq_id,
        })
        prob_cols = pd.DataFrame(probs_arr, columns=["E", "H", "L"])
        df_out = pd.concat([df_out, prob_cols], axis=1)
        df_out.to_csv(job_dir / "ResultThreeState.csv", mode="a",
                      index=True, header=not (job_dir / "ResultThreeState.csv").stat().st_size
                      if (job_dir / "ResultThreeState.csv").exists() else True)

    return results
