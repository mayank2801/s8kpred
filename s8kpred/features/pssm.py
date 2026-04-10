"""
s8kpred/features/pssm.py
-------------------------
Generate PSSM files with PSI-BLAST and extract sliding-window features.

This is a clean, standalone rewrite of PSSMBasedFeaturesGeneration.py that:
  - accepts explicit paths for psiblast and the BLAST database
  - writes only to the job directory provided
  - exposes a simple function-level API
  - avoids any global state
"""
from __future__ import annotations

import csv
import multiprocessing
import os
import re
import subprocess
import time
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from s8kpred.config import WINDOW_SIZE, PADDING_LEN, PADDING_RESIDUE
from s8kpred.utils.fasta import parse_fasta
from s8kpred.utils.io import ensure_csv

# ── Column headers for the 17-residue window PSSM feature CSV ─────────────
_AA_ORDER = list("ARNDCQEGHILKMFPSTWYV")

def _build_csv_header() -> List[str]:
    cols = ["ID", "sequence"]
    for r in range(1, WINDOW_SIZE + 1):
        for aa in _AA_ORDER:
            cols.append(f"R{r}_P_{aa}")
    return cols


# ── Low-level PSSM reader ──────────────────────────────────────────────────

def read_pssm_matrix(pssm_path: Path) -> Tuple[str, np.ndarray]:
    """
    Parse a PSI-BLAST ASCII PSSM file.

    Returns
    -------
    amino_acid_sequence : str
    pssm_matrix         : np.ndarray, shape (L, 20)
    """
    numeric_pattern = re.compile(r"-?\d+")
    amino_acids: List[str] = []
    matrix: List[List[int]] = []
    header_skipped = False

    with open(pssm_path) as fh:
        for line in fh:
            if not header_skipped:
                if line.startswith("Last position-specific scoring matrix computed"):
                    header_skipped = True
                continue
            cols = line.split()
            if not cols or not cols[0].isdigit():
                continue
            amino_acids.append(cols[1])
            if len(cols) >= 22:
                matrix.append([int(x) for x in cols[2:22]])

    if not matrix:
        raise ValueError(f"No PSSM data found in {pssm_path}")

    return "".join(amino_acids), np.array(matrix, dtype=float)


# ── PSI-BLAST runner ───────────────────────────────────────────────────────

def run_psiblast(
    fasta_path: Path,
    pssm_out: Path,
    database: str,
    psiblast_bin: str = "psiblast",
    num_iterations: int = 3,
    num_threads: Optional[int] = None,
    timeout: int = 3600,
) -> bool:
    """
    Run PSI-BLAST for a single sequence.

    Returns True on success, False on failure.
    """
    if num_threads is None:
        num_threads = multiprocessing.cpu_count()

    cmd = [
        psiblast_bin,
        "-query",          str(fasta_path),
        "-db",             database,
        "-num_iterations", str(num_iterations),
        "-num_threads",    str(num_threads),
        "-out_ascii_pssm", str(pssm_out),
    ]

    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=timeout,
        )
        return pssm_out.exists() and pssm_out.stat().st_size > 0
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        return False


# ── Main entry points ──────────────────────────────────────────────────────

def generate_pssm_files(
    job_dir: Path,
    fasta_path: Path,
    database: str,
    psiblast_bin: str = "psiblast",
    num_iterations: int = 3,
    num_threads: Optional[int] = None,
    verbose: bool = True,
) -> List[Path]:
    """
    Run PSI-BLAST for every sequence in *fasta_path*.

    Parameters
    ----------
    job_dir       : job output directory (must already exist)
    fasta_path    : input FASTA file
    database      : path to BLAST database (prefix, no extension)
    psiblast_bin  : path to the psiblast executable
    num_iterations: PSI-BLAST iterations (default 3)
    num_threads   : CPU threads (default: all available)
    verbose       : print progress to stdout

    Returns
    -------
    List of successfully generated .pssm file paths.
    """
    pssm_dir = job_dir / "pssm_outputs"
    pssm_dir.mkdir(exist_ok=True)
    log_path = job_dir / "log.dat"

    generated: List[Path] = []

    with open(log_path, "a") as log:
        log.write("=== PSSM Calculation ===\n")
        log.write("SeqID\tTimeSec\tStatus\n")

        for idx, record in enumerate(parse_fasta(fasta_path), start=1):
            seq_id  = f"Seq_{idx}"
            tmp_fa  = pssm_dir / f"{seq_id}.fasta"
            pssm_out = pssm_dir / f"{seq_id}.pssm"

            # Write single-sequence FASTA
            tmp_fa.write_text(f">{record.id}\n{record.seq}\n")

            if verbose:
                print(f"  [PSSM] {record.id} ({len(record.seq)} aa) ...", end=" ", flush=True)

            t0 = time.time()
            ok = run_psiblast(
                fasta_path=tmp_fa,
                pssm_out=pssm_out,
                database=database,
                psiblast_bin=psiblast_bin,
                num_iterations=num_iterations,
                num_threads=num_threads,
            )
            elapsed = round(time.time() - t0, 2)

            # Remove temp single-seq FASTA
            tmp_fa.unlink(missing_ok=True)

            if ok:
                generated.append(pssm_out)
                status = "SUCCESS"
                if verbose:
                    print(f"done ({elapsed}s)")
            else:
                status = "FAILED"
                if verbose:
                    print("FAILED")

            log.write(f"{seq_id}\t{elapsed}\t{status}\n")

    return generated


def extract_pssm_features(
    job_dir: Path,
    verbose: bool = True,
) -> Path:
    """
    Read all .pssm files in *job_dir/pssm_outputs/* and write a sliding-window
    feature CSV to *job_dir/PSSM_Features_ML_17W.csv*.

    Returns the path to the CSV file.
    """
    pssm_dir   = job_dir / "pssm_outputs"
    output_csv = job_dir / "PSSM_Features_ML_17W.csv"
    header     = _build_csv_header()

    ensure_csv(output_csv, header)

    pssm_files = sorted(pssm_dir.glob("*.pssm"))
    if not pssm_files:
        raise FileNotFoundError(f"No .pssm files found in {pssm_dir}")

    with open(output_csv, "a", newline="") as fh:
        writer = csv.writer(fh)

        for pssm_path in pssm_files:
            seq_id = pssm_path.stem  # e.g. "Seq_1"
            if verbose:
                print(f"  [Features] Extracting from {pssm_path.name} ...", end=" ", flush=True)

            try:
                aa_seq, pssm = read_pssm_matrix(pssm_path)
            except Exception as exc:
                print(f"ERROR: {exc}")
                continue

            # Sigmoid-scale the first 20 PSSM columns
            scaled = 1.0 / (1.0 + np.exp(-pssm[:, :20]))

            # Padding: 8 random-low rows on each side (mimics original behaviour)
            pad_rows = np.round(np.random.uniform(0, 0.3, (PADDING_LEN, 20)), 4)
            padded   = np.vstack([pad_rows, scaled, pad_rows])  # (L+16) × 20
            flat     = padded.flatten()                          # row-major

            padded_seq = PADDING_RESIDUE * PADDING_LEN + aa_seq + PADDING_RESIDUE * PADDING_LEN

            n_written = 0
            for j in range(len(padded_seq) - WINDOW_SIZE + 1):
                window = padded_seq[j: j + WINDOW_SIZE]
                if len(window) != WINDOW_SIZE:
                    continue
                start   = j * 20
                features = flat[start: start + WINDOW_SIZE * 20]
                row = [seq_id, window] + [round(x, 4) for x in features.tolist()]
                writer.writerow(row)
                n_written += 1

            if verbose:
                print(f"{n_written} windows")

    return output_csv
