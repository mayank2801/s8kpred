"""
s8kpred/api.py
--------------
High-level public API.  This is what end-users and the CLI call.
"""
from __future__ import annotations

import io
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional

from s8kpred.config import (
    PSIBLAST_DEFAULT,
    BLASTDB_DEFAULT,
    PSIBLAST_ITERATIONS,
    MAP_8TO3,
)
from s8kpred.features.pssm import generate_pssm_files, extract_pssm_features
from s8kpred.predictors.three_state import predict_3state
from s8kpred.predictors.eight_state import predict_8state
from s8kpred.utils.fasta import parse_fasta, SeqRecord, validate_sequence
from s8kpred.utils.io import make_job_dir


# ── Result container ───────────────────────────────────────────────────────

class PredictionResult:
    """
    Container returned by :func:`predict` / :func:`predict_file`.

    Attributes
    ----------
    results_3state : dict   {seq_id: secondary_structure_string (H/E/L)}
    results_8state : dict   {seq_id: secondary_structure_string (H/G/I/E/B/T/S/L)}
    job_dir        : Path   directory containing all output files
    """

    def __init__(
        self,
        results_3state: Dict[str, str],
        results_8state: Dict[str, str],
        job_dir: Path,
    ):
        self.results_3state = results_3state
        self.results_8state = results_8state
        self.job_dir        = job_dir

    def __repr__(self) -> str:
        n = len(self.results_3state)
        return (
            f"PredictionResult({n} sequence(s), job_dir={str(self.job_dir)!r})"
        )

    def summary(self) -> str:
        lines = [f"S8kPred results  —  job directory: {self.job_dir}\n"]
        for seq_id, ss3 in self.results_3state.items():
            ss8 = self.results_8state.get(seq_id, "")
            lines.append(f"  Sequence : {seq_id}  ({len(ss3)} residues)")
            lines.append(f"  3-state  : {ss3}")
            lines.append(f"  8-state  : {ss8}")
            lines.append("")
        return "\n".join(lines)


# ── Core prediction pipeline ───────────────────────────────────────────────

def _run_pipeline(
    records: List[SeqRecord],
    output_dir: str | Path,
    job_id: Optional[str],
    blastdb: str,
    psiblast_bin: str,
    num_iterations: int,
    num_threads: Optional[int],
    model_3state: Optional[Path],
    model_8state: Optional[Path],
    run_3state: bool,
    run_8state: bool,
    make_plot: bool,
    verbose: bool,
) -> PredictionResult:

    job_dir = make_job_dir(base=output_dir, job_id=job_id)

    # Write combined FASTA for PSSM generation
    fasta_path = job_dir / "FASTA" / "input_sequence.fasta"
    with open(fasta_path, "w") as fh:
        for rec in records:
            fh.write(f">{rec.id} {rec.description}\n{rec.seq}\n")

    # ── 1. PSSM generation ─────────────────────────────────────────────
    if verbose:
        print(f"\n[1/3] Generating PSSM files  (job: {job_dir.name})")

    if not blastdb:
        raise ValueError(
            "A BLAST database path is required (--blastdb / S8KPRED_BLASTDB).\n"
            "Download UniRef50 and point s8kpred at it."
        )

    generate_pssm_files(
        job_dir=job_dir,
        fasta_path=fasta_path,
        database=blastdb,
        psiblast_bin=psiblast_bin,
        num_iterations=num_iterations,
        num_threads=num_threads,
        verbose=verbose,
    )

    if verbose:
        print("\n[2/3] Extracting PSSM features")

    pssm_csv = extract_pssm_features(job_dir=job_dir, verbose=verbose)

    # ── 2. Predictions ─────────────────────────────────────────────────
    if verbose:
        print("\n[3/3] Running secondary structure predictions")

    kwargs3 = {"model_path": model_3state} if model_3state else {}
    kwargs8 = {"model_path": model_8state} if model_8state else {}

    results_3 = predict_3state(records, pssm_csv, job_dir, verbose=verbose, **kwargs3) if run_3state else {}
    results_8 = predict_8state(records, pssm_csv, job_dir, verbose=verbose, **kwargs8) if run_8state else {}

    # ── 3. Cartoon plots ───────────────────────────────────────────────
    if make_plot and results_3:
        try:
            from s8kpred.plotting.cartoon import plot_secondary_structure_cartoon
            from s8kpred.config import MAP_8TO3

            for idx, (seq_id, ss3) in enumerate(results_3.items(), start=1):
                plot_path = job_dir / f"Seq_{idx}_cartoon.png"
                # Convert H/E/L → H/E for cartoon (L shown as connector)
                ss_cartoon = "".join(
                    "H" if c == "H" else ("E" if c == "E" else "L") for c in ss3
                )
                plot_secondary_structure_cartoon(ss_cartoon, plot_path)
                if verbose:
                    print(f"  [Plot] Saved cartoon: {plot_path.name}")
        except ImportError:
            if verbose:
                print("  [Plot] Skipping cartoon (biotite not installed). "
                      "Run: pip install s8kpred[plot]")

    if verbose:
        print(f"\n✓ Done.  Results saved to: {job_dir}")

    return PredictionResult(results_3, results_8, job_dir)


# ── Public convenience wrappers ────────────────────────────────────────────

def predict(
    sequence: str,
    seq_id: str = "query",
    output_dir: str | Path = "s8kpred_jobs",
    job_id: Optional[str] = None,
    blastdb: str = BLASTDB_DEFAULT,
    psiblast_bin: str = PSIBLAST_DEFAULT,
    num_iterations: int = PSIBLAST_ITERATIONS,
    num_threads: Optional[int] = None,
    model_3state: Optional[Path] = None,
    model_8state: Optional[Path] = None,
    run_3state: bool = True,
    run_8state: bool = True,
    make_plot: bool = True,
    verbose: bool = True,
) -> PredictionResult:
    """
    Predict secondary structure for a single amino-acid sequence string.

    Parameters
    ----------
    sequence      : raw amino acid string (single-letter codes)
    seq_id        : identifier for this sequence (used in output files)
    output_dir    : parent directory for job output folders
    job_id        : explicit job folder name; auto-generated if None
    blastdb       : path prefix to the PSI-BLAST database (required)
    psiblast_bin  : path to the psiblast executable
    num_iterations: PSI-BLAST iterations (default 3)
    num_threads   : CPU threads for PSI-BLAST (default: all available)
    model_3state  : override path for the 3-state XGBoost model
    model_8state  : override path for the 8-state XGBoost model
    run_3state    : run 3-state predictor (default True)
    run_8state    : run 8-state predictor (default True)
    make_plot     : generate cartoon PNG (requires biotite)
    verbose       : print progress

    Returns
    -------
    PredictionResult
    """
    records = [SeqRecord(id=seq_id, description=seq_id, seq=sequence.upper())]

    for w in validate_sequence(sequence):
        if verbose:
            print(f"  [Warning] {w}")

    return _run_pipeline(
        records=records,
        output_dir=output_dir,
        job_id=job_id,
        blastdb=blastdb,
        psiblast_bin=psiblast_bin,
        num_iterations=num_iterations,
        num_threads=num_threads,
        model_3state=model_3state,
        model_8state=model_8state,
        run_3state=run_3state,
        run_8state=run_8state,
        make_plot=make_plot,
        verbose=verbose,
    )


def predict_file(
    fasta_file: str | Path,
    output_dir: str | Path = "s8kpred_jobs",
    job_id: Optional[str] = None,
    blastdb: str = BLASTDB_DEFAULT,
    psiblast_bin: str = PSIBLAST_DEFAULT,
    num_iterations: int = PSIBLAST_ITERATIONS,
    num_threads: Optional[int] = None,
    model_3state: Optional[Path] = None,
    model_8state: Optional[Path] = None,
    run_3state: bool = True,
    run_8state: bool = True,
    make_plot: bool = True,
    verbose: bool = True,
) -> PredictionResult:
    """
    Predict secondary structure for all sequences in a FASTA file.

    Parameters
    ----------
    fasta_file : path to FASTA file (single or multi-sequence)
    (all other parameters same as :func:`predict`)

    Returns
    -------
    PredictionResult
    """
    fasta_file = Path(fasta_file)
    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    records = list(parse_fasta(fasta_file))
    if not records:
        raise ValueError(f"No sequences found in {fasta_file}")

    if verbose:
        print(f"Loaded {len(records)} sequence(s) from {fasta_file.name}")

    for rec in records:
        for w in validate_sequence(rec.seq):
            if verbose:
                print(f"  [Warning] {rec.id}: {w}")

    return _run_pipeline(
        records=records,
        output_dir=output_dir,
        job_id=job_id,
        blastdb=blastdb,
        psiblast_bin=psiblast_bin,
        num_iterations=num_iterations,
        num_threads=num_threads,
        model_3state=model_3state,
        model_8state=model_8state,
        run_3state=run_3state,
        run_8state=run_8state,
        make_plot=make_plot,
        verbose=verbose,
    )


def predict_fasta_string(
    fasta_string: str,
    output_dir: str | Path = "s8kpred_jobs",
    job_id: Optional[str] = None,
    **kwargs,
) -> PredictionResult:
    """
    Predict secondary structure from a raw FASTA-formatted string.

    Parameters
    ----------
    fasta_string : FASTA-formatted text (one or more sequences)
    (all other keyword arguments passed to :func:`predict_file`)

    Returns
    -------
    PredictionResult
    """
    records = list(parse_fasta(fasta_string))
    if not records:
        raise ValueError("No sequences found in the provided FASTA string.")

    verbose = kwargs.pop("verbose", True)
    if verbose:
        print(f"Parsed {len(records)} sequence(s) from string.")

    return _run_pipeline(
        records=records,
        output_dir=output_dir,
        job_id=job_id,
        verbose=verbose,
        **{k: kwargs.get(k, v) for k, v in dict(
            blastdb=BLASTDB_DEFAULT,
            psiblast_bin=PSIBLAST_DEFAULT,
            num_iterations=PSIBLAST_ITERATIONS,
            num_threads=None,
            model_3state=None,
            model_8state=None,
            run_3state=True,
            run_8state=True,
            make_plot=True,
        ).items()},
    )
