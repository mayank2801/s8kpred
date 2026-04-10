"""
s8kpred/cli.py
--------------
Command-line interface.

Usage examples
--------------
# Single FASTA file
s8kpred predict -i protein.fasta --blastdb /data/uniref50

# Multiple FASTA files
s8kpred predict -i seq1.fasta seq2.fasta --blastdb /data/uniref50

# Inline sequence string
s8kpred predict --sequence MKTAYIAKQRQ... --blastdb /data/uniref50

# Custom output directory and job name
s8kpred predict -i input.fasta --blastdb /data/uniref50 -o ./results --job myrun

# Skip 8-state prediction
s8kpred predict -i input.fasta --blastdb /data/uniref50 --no-8state

# Check installed version
s8kpred --version
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from s8kpred import __version__
from s8kpred.config import PSIBLAST_DEFAULT, BLASTDB_DEFAULT, PSIBLAST_ITERATIONS


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="s8kpred",
        description=(
            "S8kPred — Protein Secondary Structure Prediction\n"
            "Predicts 3-state (H/E/L) and 8-state (H/G/I/E/B/T/S/L) "
            "secondary structure from amino acid sequences."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"s8kpred {__version__}")

    sub = parser.add_subparsers(dest="command", required=True)

    # ── predict sub-command ───────────────────────────────────────────
    pred = sub.add_parser(
        "predict",
        help="Run secondary structure prediction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Predict secondary structure.\n\n"
            "Supply either --input (FASTA files) or --sequence (raw string)."
        ),
    )

    # Input
    input_group = pred.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-i", "--input",
        nargs="+",
        metavar="FASTA",
        help="One or more FASTA files (single- or multi-sequence).",
    )
    input_group.add_argument(
        "-s", "--sequence",
        metavar="SEQ",
        help="Raw amino-acid sequence string (single-letter codes).",
    )

    # BLAST
    pred.add_argument(
        "--blastdb",
        default=BLASTDB_DEFAULT,
        metavar="DB",
        help=(
            "Path prefix to PSI-BLAST database (e.g. /data/uniref50/uniref50).\n"
            "Can also be set via the S8KPRED_BLASTDB environment variable."
        ),
    )
    pred.add_argument(
        "--psiblast",
        default=PSIBLAST_DEFAULT,
        metavar="BIN",
        help="Path to psiblast executable (default: 'psiblast' on PATH).",
    )
    pred.add_argument(
        "--iterations",
        type=int,
        default=PSIBLAST_ITERATIONS,
        metavar="N",
        help=f"PSI-BLAST iterations (default: {PSIBLAST_ITERATIONS}).",
    )
    pred.add_argument(
        "--threads",
        type=int,
        default=None,
        metavar="N",
        help="CPU threads for PSI-BLAST (default: all available).",
    )

    # Output
    pred.add_argument(
        "-o", "--output-dir",
        default="s8kpred_jobs",
        metavar="DIR",
        help="Parent directory for job output (default: ./s8kpred_jobs).",
    )
    pred.add_argument(
        "--job",
        default=None,
        metavar="ID",
        help="Job identifier / sub-folder name (auto-generated if omitted).",
    )

    # Model overrides
    pred.add_argument(
        "--model-3state",
        default=None,
        metavar="PATH",
        help="Override path to 3-state XGBoost model (.json).",
    )
    pred.add_argument(
        "--model-8state",
        default=None,
        metavar="PATH",
        help="Override path to 8-state XGBoost model (.ubj).",
    )

    # Prediction mode flags
    pred.add_argument(
        "--no-3state",
        action="store_true",
        help="Skip 3-state prediction.",
    )
    pred.add_argument(
        "--no-8state",
        action="store_true",
        help="Skip 8-state prediction.",
    )
    pred.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip cartoon plot generation.",
    )

    # Verbosity
    pred.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress progress output.",
    )

    # Sequence ID (used only with --sequence)
    pred.add_argument(
        "--id",
        default="query",
        metavar="ID",
        help="Sequence ID when using --sequence (default: 'query').",
    )

    return parser


def main(argv=None):
    parser = _build_parser()
    args   = parser.parse_args(argv)

    if args.command == "predict":
        _cmd_predict(args)


def _cmd_predict(args):
    from s8kpred.api import predict, predict_file

    verbose  = not args.quiet
    run_3    = not args.no_3state
    run_8    = not args.no_8state
    make_plot = not args.no_plot

    common = dict(
        output_dir     = args.output_dir,
        job_id         = args.job,
        blastdb        = args.blastdb,
        psiblast_bin   = args.psiblast,
        num_iterations = args.iterations,
        num_threads    = args.threads,
        model_3state   = Path(args.model_3state) if args.model_3state else None,
        model_8state   = Path(args.model_8state) if args.model_8state else None,
        run_3state     = run_3,
        run_8state     = run_8,
        make_plot      = make_plot,
        verbose        = verbose,
    )

    try:
        if args.sequence:
            result = predict(sequence=args.sequence, seq_id=args.id, **common)
            if verbose:
                print(result.summary())

        else:
            # One or more FASTA files
            for fasta_path in args.input:
                if verbose:
                    print(f"\n{'='*60}")
                    print(f"Processing: {fasta_path}")
                    print('='*60)
                result = predict_file(fasta_file=fasta_path, **common)
                if verbose:
                    print(result.summary())

    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        sys.exit(130)


if __name__ == "__main__":
    main()
