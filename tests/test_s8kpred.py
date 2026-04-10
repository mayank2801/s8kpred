"""
tests/test_s8kpred.py
---------------------
Unit tests for S8kPred.  These tests do NOT require PSI-BLAST or a BLAST
database — all external calls are mocked so the test suite runs in any CI
environment.
"""
from __future__ import annotations

import csv
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

# ── helpers ────────────────────────────────────────────────────────────────

SAMPLE_SEQ   = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"
SAMPLE_FASTA = f">seq1 test protein\n{SAMPLE_SEQ}\n"

MULTI_FASTA = textwrap.dedent("""\
    >prot_A
    ACDEFGHIKLMNPQRSTVWY
    >prot_B
    MKTAYIAKQRQISFVKSHFS
    >prot_C
    ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY
""")


# ══════════════════════════════════════════════════════════════════════════
# utils/fasta.py
# ══════════════════════════════════════════════════════════════════════════

class TestFastaParser:
    def test_single_record(self):
        from s8kpred.utils.fasta import parse_fasta
        records = list(parse_fasta(SAMPLE_FASTA))
        assert len(records) == 1
        assert records[0].id  == "seq1"
        assert records[0].seq == SAMPLE_SEQ.upper()

    def test_multi_record(self):
        from s8kpred.utils.fasta import parse_fasta
        records = list(parse_fasta(MULTI_FASTA))
        assert len(records) == 3
        assert records[0].id  == "prot_A"
        assert records[1].id  == "prot_B"
        assert records[2].id  == "prot_C"

    def test_sequence_uppercased(self):
        from s8kpred.utils.fasta import parse_fasta
        records = list(parse_fasta(">p\nacdefg\n"))
        assert records[0].seq == "ACDEFG"

    def test_empty_string_raises(self):
        from s8kpred.utils.fasta import parse_fasta
        records = list(parse_fasta(""))
        assert records == []

    def test_file_path(self, tmp_path):
        from s8kpred.utils.fasta import parse_fasta
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(SAMPLE_FASTA)
        records = list(parse_fasta(fasta_file))
        assert len(records) == 1
        assert records[0].seq == SAMPLE_SEQ.upper()

    def test_len(self):
        from s8kpred.utils.fasta import parse_fasta
        records = list(parse_fasta(SAMPLE_FASTA))
        assert len(records[0]) == len(SAMPLE_SEQ)

    def test_multiline_sequence(self):
        from s8kpred.utils.fasta import parse_fasta
        fa = ">p\nACDEF\nGHIKL\nMNPQR\n"
        records = list(parse_fasta(fa))
        assert records[0].seq == "ACDEFGHIKLMNPQR"


class TestValidateSequence:
    def test_valid_sequence(self):
        from s8kpred.utils.fasta import validate_sequence
        warnings = validate_sequence("MKTAYIAKQR")
        assert warnings == []

    def test_short_sequence_warning(self):
        from s8kpred.utils.fasta import validate_sequence
        warnings = validate_sequence("MKTAY")
        assert any("short" in w.lower() for w in warnings)

    def test_non_standard_characters(self):
        from s8kpred.utils.fasta import validate_sequence
        warnings = validate_sequence("MKTAY123")
        assert any("Non-standard" in w for w in warnings)

    def test_x_is_handled(self):
        from s8kpred.utils.fasta import validate_sequence
        # X is a valid ambiguity code — should NOT trigger non-standard warning
        warnings = validate_sequence("MKTAYXAKQR")
        assert not any("Non-standard" in w for w in warnings)


# ══════════════════════════════════════════════════════════════════════════
# utils/io.py
# ══════════════════════════════════════════════════════════════════════════

class TestMakeJobDir:
    def test_creates_subdirs(self, tmp_path):
        from s8kpred.utils.io import make_job_dir
        job_dir = make_job_dir(base=tmp_path, job_id="testjob")
        assert job_dir.exists()
        assert (job_dir / "pssm_outputs").is_dir()
        assert (job_dir / "FASTA").is_dir()

    def test_auto_job_id(self, tmp_path):
        from s8kpred.utils.io import make_job_dir
        job_dir = make_job_dir(base=tmp_path)
        assert job_dir.exists()
        # Auto ID contains a timestamp prefix
        assert "_" in job_dir.name

    def test_explicit_job_id(self, tmp_path):
        from s8kpred.utils.io import make_job_dir
        job_dir = make_job_dir(base=tmp_path, job_id="myjob")
        assert job_dir.name == "myjob"


class TestEnsureCsv:
    def test_creates_file_with_header(self, tmp_path):
        from s8kpred.utils.io import ensure_csv
        path = tmp_path / "out.csv"
        ensure_csv(path, ["col1", "col2", "col3"])
        assert path.exists()
        with open(path) as f:
            reader = csv.reader(f)
            header = next(reader)
        assert header == ["col1", "col2", "col3"]

    def test_does_not_overwrite_existing(self, tmp_path):
        from s8kpred.utils.io import ensure_csv
        path = tmp_path / "out.csv"
        path.write_text("existing,content\n1,2\n")
        ensure_csv(path, ["col1", "col2"])
        content = path.read_text()
        assert content.startswith("existing,content")


class TestWriteOutputFiles:
    def _make_probs(self, n, k):
        p = np.random.dirichlet(np.ones(k), size=n)
        return p

    def test_write_ss2_3state(self, tmp_path):
        from s8kpred.utils.io import write_ss2
        path    = tmp_path / "result.ss2"
        residues = list("MKTAYIAKQR")
        ss       = list("HHHEEELLLH")
        probs    = self._make_probs(10, 3)
        write_ss2(path, residues, ss, probs, header="test_seq", n_classes=3)
        content = path.read_text()
        assert "S8kPred VFORMAT" in content
        assert "test_seq" in content
        assert len(content.strip().splitlines()) == 11  # header + 10 residues

    def test_write_ss2_8state(self, tmp_path):
        from s8kpred.utils.io import write_ss2
        path     = tmp_path / "result8.ss2"
        residues = list("MKTAYIAKQR")
        ss       = list("HHHEEELLLH")
        probs    = self._make_probs(10, 8)
        write_ss2(path, residues, ss, probs, header="seq8", n_classes=8)
        content = path.read_text()
        assert "B  E  G  H  I  L  S  T" in content

    def test_write_fasta_result(self, tmp_path):
        from s8kpred.utils.io import write_fasta_result
        path     = tmp_path / "result.fas"
        residues = list("MKTAYIAKQR")
        ss       = list("HHHEEELLLH")
        write_fasta_result(path, residues, ss, header="myseq")
        lines = path.read_text().strip().splitlines()
        assert lines[0].startswith("> S8kPred Fasta myseq")
        assert lines[1] == "MKTAYIAKQR"
        assert lines[2] == "HHHEEELLLH"


# ══════════════════════════════════════════════════════════════════════════
# features/pssm.py
# ══════════════════════════════════════════════════════════════════════════

class TestReadPssmMatrix:
    def _write_fake_pssm(self, tmp_path: Path) -> Path:
        """Write a minimal valid PSI-BLAST PSSM file."""
        pssm_file = tmp_path / "Seq_1.pssm"
        lines = [
            "Some header line\n",
            "Last position-specific scoring matrix computed, nah blah\n",
            "           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V\n",
        ]
        aa_seq = "MKTAY"
        for i, aa in enumerate(aa_seq, 1):
            vals = " ".join(str(i % 5 - 2) for _ in range(20))
            lines.append(f"{i:4d} {aa}  {vals}  0.00  0.00\n")
        lines.append("\n\nsome trailing info\n")
        pssm_file.write_text("".join(lines))
        return pssm_file

    def test_reads_sequence_and_matrix(self, tmp_path):
        from s8kpred.features.pssm import read_pssm_matrix
        pssm_path = self._write_fake_pssm(tmp_path)
        seq, mat = read_pssm_matrix(pssm_path)
        assert seq == "MKTAY"
        assert mat.shape == (5, 20)

    def test_raises_on_empty_pssm(self, tmp_path):
        from s8kpred.features.pssm import read_pssm_matrix
        empty = tmp_path / "empty.pssm"
        empty.write_text("no data here\n")
        with pytest.raises(ValueError, match="No PSSM data"):
            read_pssm_matrix(empty)


class TestExtractPssmFeatures:
    def _write_fake_pssm(self, pssm_dir: Path, name: str, seq: str):
        pssm_file = pssm_dir / f"{name}.pssm"
        lines = ["Last position-specific scoring matrix computed\n"]
        for i, aa in enumerate(seq, 1):
            vals = " ".join("1" for _ in range(20))
            lines.append(f"{i:4d} {aa}  {vals}  0.00  0.00\n")
        pssm_file.write_text("".join(lines))

    def test_csv_written(self, tmp_path):
        from s8kpred.features.pssm import extract_pssm_features
        (tmp_path / "pssm_outputs").mkdir()
        seq = "A" * 20
        self._write_fake_pssm(tmp_path / "pssm_outputs", "Seq_1", seq)
        csv_path = extract_pssm_features(tmp_path, verbose=False)
        assert csv_path.exists()
        df = pd.read_csv(csv_path)
        assert "ID" in df.columns
        assert "sequence" in df.columns
        assert len(df) == len(seq)   # 1 window per residue (after padding)

    def test_no_pssm_files_raises(self, tmp_path):
        from s8kpred.features.pssm import extract_pssm_features
        (tmp_path / "pssm_outputs").mkdir()
        with pytest.raises(FileNotFoundError):
            extract_pssm_features(tmp_path, verbose=False)


# ══════════════════════════════════════════════════════════════════════════
# features/propensity.py
# ══════════════════════════════════════════════════════════════════════════

def _make_fake_propensity_3state(tmp_path: Path) -> Path:
    """Create a minimal propensity CSV for 3-state tests."""
    import itertools
    AAs = list("ACDEFGHIKLMNPQRSTVWY") + ["G"]
    rows = []
    for a, b, c in itertools.product("ACDEFG", repeat=3):
        rows.append({"TriPeptide": a + b + c, "dummy": 0, "H": 0.3, "E": 0.3, "L": 0.4})
    df = pd.DataFrame(rows)
    path = tmp_path / "prop3.csv"
    df.to_csv(path, index=False)
    return path


def _make_fake_binary_table(tmp_path: Path) -> Path:
    """Create a minimal binary lookup table."""
    import itertools
    rows = []
    for a, b, c in itertools.product("ACDEFG", repeat=3):
        row = {"Tripeptide": a + b + c}
        row.update({f"f{i}": 0 for i in range(60)})
        rows.append(row)
    df = pd.DataFrame(rows)
    path = tmp_path / "binary.csv"
    df.to_csv(path, index=False)
    return path


class TestBuildFeatures:
    def test_3state_returns_dataframe(self, tmp_path):
        from s8kpred.features.propensity import build_3state_features
        from s8kpred.utils.fasta import SeqRecord
        prop_csv   = _make_fake_propensity_3state(tmp_path)
        binary_csv = _make_fake_binary_table(tmp_path)
        records    = [SeqRecord(id="p1", description="p1", seq="ACDEFG" * 5)]
        df, bdf    = build_3state_features(records, prop_csv, binary_csv)
        assert isinstance(df, pd.DataFrame)
        assert "ID" in df.columns
        assert "sequence" in df.columns
        assert "Residue9th" in df.columns
        assert len(df) == len(records[0].seq)

    def test_3state_window_size(self, tmp_path):
        from s8kpred.features.propensity import build_3state_features
        from s8kpred.utils.fasta import SeqRecord
        prop_csv   = _make_fake_propensity_3state(tmp_path)
        binary_csv = _make_fake_binary_table(tmp_path)
        records    = [SeqRecord(id="p1", description="p1", seq="ACDEFG" * 3)]
        df, _      = build_3state_features(records, prop_csv, binary_csv)
        # Every window entry must be exactly 17 characters
        assert (df["sequence"].str.len() == 17).all()


# ══════════════════════════════════════════════════════════════════════════
# config.py
# ══════════════════════════════════════════════════════════════════════════

class TestConfig:
    def test_window_size(self):
        from s8kpred.config import WINDOW_SIZE
        assert WINDOW_SIZE == 17

    def test_padding_len(self):
        from s8kpred.config import PADDING_LEN
        assert PADDING_LEN == 8

    def test_3state_map_covers_all(self):
        from s8kpred.config import IDX2CHAR_3STATE
        assert set(IDX2CHAR_3STATE.values()) == {"E", "H", "L"}

    def test_8state_map_covers_all(self):
        from s8kpred.config import IDX2CHAR_8STATE
        assert set(IDX2CHAR_8STATE.values()) == set("BGHIELST")


# ══════════════════════════════════════════════════════════════════════════
# cli.py
# ══════════════════════════════════════════════════════════════════════════

class TestCLI:
    def test_version(self, capsys):
        from s8kpred.cli import main
        with pytest.raises(SystemExit) as exc:
            main(["--version"])
        assert exc.value.code == 0

    def test_predict_requires_blastdb(self, tmp_path, capsys):
        """Without --blastdb the pipeline should raise ValueError."""
        from s8kpred.cli import main
        import os
        # Make sure env var is not set
        os.environ.pop("S8KPRED_BLASTDB", None)

        fasta = tmp_path / "test.fasta"
        fasta.write_text(SAMPLE_FASTA)

        with pytest.raises(SystemExit) as exc:
            main(["predict", "-i", str(fasta)])
        assert exc.value.code == 1

    def test_predict_sequence_flag_missing_blastdb(self, capsys):
        from s8kpred.cli import main
        import os
        os.environ.pop("S8KPRED_BLASTDB", None)
        with pytest.raises(SystemExit) as exc:
            main(["predict", "--sequence", SAMPLE_SEQ, "--id", "myq"])
        assert exc.value.code == 1

    def test_no_input_errors(self, capsys):
        from s8kpred.cli import main
        with pytest.raises(SystemExit):
            main(["predict"])


# ══════════════════════════════════════════════════════════════════════════
# api.py  (mocked pipeline — no real BLAST needed)
# ══════════════════════════════════════════════════════════════════════════

class TestApiMocked:
    """Test the public API with all external calls mocked."""

    def _make_fake_pssm_csv(self, job_dir: Path, seq_id="Seq_1", n_windows=30):
        """Write a fake PSSM feature CSV that the predictors will read."""
        from s8kpred.features.pssm import _build_csv_header
        header = _build_csv_header()
        csv_path = job_dir / "PSSM_Features_ML_17W.csv"
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            for i in range(n_windows):
                window = ("ACDEFGHIKLMNPQRST"[:17])
                row = [seq_id, window] + [round(np.random.rand(), 4) for _ in range(17 * 20)]
                w.writerow(row)
        return csv_path

    @patch("s8kpred.features.pssm.run_psiblast", return_value=True)
    @patch("s8kpred.features.pssm.read_pssm_matrix")
    def test_predict_file_runs(self, mock_read_pssm, mock_run_psiblast, tmp_path):
        """predict_file() should complete without errors when BLAST is mocked."""
        # Stub read_pssm_matrix to return a fake matrix
        mock_read_pssm.return_value = (
            SAMPLE_SEQ,
            np.ones((len(SAMPLE_SEQ), 20), dtype=float),
        )

        fasta = tmp_path / "input.fasta"
        fasta.write_text(SAMPLE_FASTA)

        # We also need to mock the XGBoost models
        mock_model = MagicMock()
        mock_model.predict.return_value = np.zeros(len(SAMPLE_SEQ), dtype=int)
        mock_model.predict_proba.return_value = np.ones(
            (len(SAMPLE_SEQ), 3)
        ) / 3.0

        with patch("xgboost.XGBClassifier") as MockXGB:
            MockXGB.return_value = mock_model
            from s8kpred.api import predict_file
            result = predict_file(
                fasta_file=fasta,
                output_dir=tmp_path / "jobs",
                blastdb="/fake/db",
                run_3state=True,
                run_8state=False,
                make_plot=False,
                verbose=False,
            )

        assert isinstance(result.results_3state, dict)
        assert result.job_dir.exists()

    def test_predict_file_not_found(self, tmp_path):
        from s8kpred.api import predict_file
        with pytest.raises(FileNotFoundError):
            predict_file(tmp_path / "nonexistent.fasta", blastdb="/db")

    def test_predict_empty_fasta(self, tmp_path):
        from s8kpred.api import predict_file
        fasta = tmp_path / "empty.fasta"
        fasta.write_text("")
        with pytest.raises(ValueError, match="No sequences"):
            predict_file(fasta, blastdb="/db")

    def test_prediction_result_repr(self, tmp_path):
        from s8kpred.api import PredictionResult
        r = PredictionResult(
            results_3state={"seq1": "HHHEEEL"},
            results_8state={"seq1": "HHHEEELT"},
            job_dir=tmp_path,
        )
        assert "1 sequence" in repr(r)

    def test_prediction_result_summary(self, tmp_path):
        from s8kpred.api import PredictionResult
        r = PredictionResult(
            results_3state={"seq1": "HHHEEEL"},
            results_8state={"seq1": "HHHEEELT"},
            job_dir=tmp_path,
        )
        summary = r.summary()
        assert "seq1" in summary
        assert "HHHEEEL" in summary
