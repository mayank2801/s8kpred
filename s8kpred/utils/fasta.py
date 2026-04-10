"""
s8kpred/utils/fasta.py
----------------------
Lightweight FASTA parser that does NOT require Biopython so the package
remains installable in minimal environments.  Biopython is still an optional
dependency (used for PSSM generation) but not required just to read FASTA.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, List


VALID_AA = re.compile(r"[^ACDEFGHIKLMNPQRSTVWYBXZUO*-]", re.IGNORECASE)


@dataclass
class SeqRecord:
    id: str
    description: str
    seq: str

    def __len__(self) -> int:
        return len(self.seq)

    def __repr__(self) -> str:
        return f"SeqRecord(id={self.id!r}, len={len(self.seq)})"


def _parse_header(header: str):
    """Split '>id desc' → (id, description)."""
    header = header.lstrip(">").strip()
    parts = header.split(None, 1)
    seq_id = parts[0] if parts else "unknown"
    desc = parts[1] if len(parts) > 1 else seq_id
    return seq_id, desc


def parse_fasta(source) -> Iterator[SeqRecord]:
    """
    Yield SeqRecord objects from a FASTA file or string.

    Parameters
    ----------
    source : str | Path | file-like
        A file path, pathlib.Path, raw FASTA string, or any object with a
        ``read()`` method.
    """
    if isinstance(source, (str, Path)):
        p = Path(source)
        if p.exists():
            text = p.read_text()
        else:
            # Treat as raw FASTA string
            text = str(source)
    elif hasattr(source, "read"):
        text = source.read()
    else:
        raise TypeError(f"Unsupported source type: {type(source)}")

    current_id = current_desc = None
    seq_parts: List[str] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                yield SeqRecord(
                    id=current_id,
                    description=current_desc,
                    seq="".join(seq_parts).upper(),
                )
            current_id, current_desc = _parse_header(line)
            seq_parts = []
        else:
            seq_parts.append(line)

    if current_id is not None and seq_parts:
        yield SeqRecord(
            id=current_id,
            description=current_desc,
            seq="".join(seq_parts).upper(),
        )


def validate_sequence(seq: str) -> List[str]:
    """Return a list of warning strings for unusual characters."""
    warnings = []
    bad = set(VALID_AA.findall(seq))
    if bad:
        warnings.append(f"Non-standard characters found (will be replaced with G): {bad}")
    if len(seq) < 10:
        warnings.append("Sequence is very short (< 10 residues); predictions may be unreliable.")
    return warnings
