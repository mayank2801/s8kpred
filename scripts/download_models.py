#!/usr/bin/env python3
"""
scripts/download_models.py
--------------------------
Download S8kPred model data files from the GitHub release and place them
in the correct location inside the installed package.

Usage
-----
    python scripts/download_models.py
    python scripts/download_models.py --release v0.1.0
    python scripts/download_models.py --dest /custom/data/dir
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path

# ── Files to download ──────────────────────────────────────────────────────
# Each entry: (filename, sha256_hex_or_None)
DATA_FILES = [
    ("TriPeptidePropensityThreeStateSecStructure2AND.csv", None),
    ("TriPeptidePropensityEightStateSecStructure.csv",     None),
    ("TripeptideBinaryTable_60.csv",                       None),
    ("model_3state.json",                                  None),
    ("model_8state.ubj",                                   None),
]

GITHUB_RELEASES_BASE = (
    "https://github.com/mayank2801/s8kpred/releases/download/{release}/{filename}"
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def download_models(release: str = "latest", dest: Path | None = None):
    if dest is None:
        # Default: package data directory
        import s8kpred
        dest = Path(s8kpred.__file__).parent / "data"

    dest.mkdir(parents=True, exist_ok=True)
    print(f"Downloading model files to: {dest}")

    if release == "latest":
        # Resolve latest release tag via GitHub API
        import json
        api_url = "https://api.github.com/repos/mayank2801/s8kpred/releases/latest"
        try:
            with urllib.request.urlopen(api_url, timeout=10) as resp:
                info = json.loads(resp.read())
                release = info["tag_name"]
        except Exception as exc:
            print(f"Could not resolve latest release: {exc}", file=sys.stderr)
            print("Please specify --release explicitly (e.g. --release v0.1.0)")
            sys.exit(1)

    print(f"Release: {release}\n")

    for filename, expected_sha256 in DATA_FILES:
        out_path = dest / filename
        if out_path.exists():
            print(f"  ✓ Already present: {filename}")
            continue

        url = GITHUB_RELEASES_BASE.format(release=release, filename=filename)
        print(f"  ↓ {filename} ... ", end="", flush=True)
        try:
            urllib.request.urlretrieve(url, out_path)
            print("done")
        except Exception as exc:
            print(f"FAILED ({exc})", file=sys.stderr)
            out_path.unlink(missing_ok=True)
            continue

        if expected_sha256:
            actual = _sha256(out_path)
            if actual != expected_sha256:
                print(f"    WARNING: checksum mismatch for {filename}!")
                print(f"      expected: {expected_sha256}")
                print(f"      got:      {actual}")

    print("\nAll done.")


def main():
    parser = argparse.ArgumentParser(description="Download S8kPred model data files.")
    parser.add_argument(
        "--release",
        default="latest",
        help="GitHub release tag (e.g. v0.1.0). Default: latest.",
    )
    parser.add_argument(
        "--dest",
        default=None,
        type=Path,
        help="Destination directory. Default: s8kpred/data/ inside the package.",
    )
    args = parser.parse_args()
    download_models(release=args.release, dest=args.dest)


if __name__ == "__main__":
    main()
