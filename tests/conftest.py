"""
conftest.py
-----------
Shared pytest fixtures available to all tests.
"""
import pytest


@pytest.fixture(scope="session")
def sample_seq():
    return "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"


@pytest.fixture(scope="session")
def sample_fasta_str(sample_seq):
    return f">test_protein Test protein sequence\n{sample_seq}\n"


@pytest.fixture
def tmp_job_dir(tmp_path):
    """A temporary job directory with the expected sub-structure."""
    (tmp_path / "pssm_outputs").mkdir()
    (tmp_path / "FASTA").mkdir()
    return tmp_path
