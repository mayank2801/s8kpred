"""
s8kpred/utils/io.py
-------------------
Helpers for job-directory management and output file writing.
"""
from __future__ import annotations

import csv
import os
import uuid
from datetime import datetime
from pathlib import Path
from typing import List


def make_job_dir(base: str | Path = "s8kpred_jobs", job_id: str | None = None) -> Path:
    """
    Create (and return) a job output directory.

    Parameters
    ----------
    base    : parent directory for all jobs
    job_id  : explicit ID; if None a timestamped UUID is generated
    """
    if job_id is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        job_id = f"{ts}_{uuid.uuid4().hex[:6]}"
    job_dir = Path(base) / job_id
    (job_dir / "pssm_outputs").mkdir(parents=True, exist_ok=True)
    (job_dir / "FASTA").mkdir(parents=True, exist_ok=True)
    return job_dir


def write_ss2(path: Path, data, ss, probabilities, header: str, n_classes: int):
    """Write a PSIPRED-style .ss2 file."""
    with open(path, "a") as fh:
        if n_classes == 3:
            fh.write(f"# S8kPred VFORMAT {header}\n")
            for i, (aa, pred) in enumerate(zip(data, ss)):
                fh.write("%4d %c %c  %6.3f %6.3f %6.3f\n" % (
                    i + 1, aa, pred,
                    probabilities[i, 0], probabilities[i, 1], probabilities[i, 2]
                ))
        else:
            fh.write(f"# S8kPred VFORMAT \n{header}\n # AA  SS  B  E  G  H  I  L  S  T \n")
            for i, (aa, pred) in enumerate(zip(data, ss)):
                fh.write("%4d %c %c  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n" % (
                    i + 1, aa, pred,
                    probabilities[i, 0], probabilities[i, 1], probabilities[i, 2],
                    probabilities[i, 3], probabilities[i, 4], probabilities[i, 5],
                    probabilities[i, 6], probabilities[i, 7]
                ))


def write_horiz(path: Path, data, ss, probabilities, header: str):
    """Write a PSIPRED-style .horiz file."""
    import numpy as np

    def chunkstring(s, n):
        return [s[i:i + n] for i in range(0, len(s), n)]

    with open(path, "a") as fh:
        fh.write(f"# S8kPred HFORMAT  {header}\n")
        sub_seqs = chunkstring("".join(str(a) for a in data), 60)
        sub_ss   = chunkstring("".join(str(s) for s in ss), 60)

        num_len = int(np.floor(len(data) / 10))
        num_seq = "".join(f"{str((i + 1) * 10):>10}" for i in range(num_len + 1))
        num_seq_chunks = chunkstring(num_seq, 60)

        conf_idxs = probabilities.argmax(-1)
        confs     = probabilities[np.arange(len(conf_idxs)), conf_idxs]
        conf_str  = "".join(str(x) for x in np.floor(confs * 10).astype(int))
        conf_chunks = chunkstring(conf_str, 60)

        for idx, subsq in enumerate(sub_seqs):
            fh.write(f"\nConf: {conf_chunks[idx]}\n")
            fh.write(f"Pred: {sub_ss[idx]}\n")
            fh.write(f"  AA: {subsq}\n")
            fh.write(f"      {num_seq_chunks[idx]}\n\n")


def write_fasta_result(path: Path, data, ss, header: str):
    """Write a pseudo-FASTA result file."""
    with open(path, "a") as fh:
        fh.write(f"> S8kPred Fasta {header}\n")
        fh.write("".join(str(a) for a in data) + "\n")
        fh.write("".join(str(s) for s in ss) + "\n")


def ensure_csv(path: Path, header: List[str]):
    """Create a CSV with header if it does not exist."""
    if not path.exists():
        with open(path, "w", newline="") as fh:
            csv.writer(fh).writerow(header)
