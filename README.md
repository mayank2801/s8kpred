# S8kPred — Protein Secondary Structure Prediction

[![PyPI](https://img.shields.io/pypi/v/s8kpred)](https://pypi.org/project/s8kpred/)
[![Python](https://img.shields.io/pypi/pyversions/s8kpred)](https://pypi.org/project/s8kpred/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**S8kPred** predicts protein secondary structure directly from amino acid sequence using XGBoost models trained on PSSM (Position-Specific Scoring Matrix) and tripeptide propensity features.

- **3-state**: Helix (H), Beta-strand (E), Coil/Loop (L)
- **8-state**: H, G, I, E, B, T, S, L (full DSSP alphabet)

---

## Requirements

| Dependency | Purpose |
| :-- | :-- |
| Python ≥ 3.9 | Runtime |
| `numpy`, `pandas`, `xgboost`, `scikit-learn` | Core ML pipeline |
| `biopython` | FASTA I/O for PSI-BLAST |
| **NCBI PSI-BLAST** | PSSM generation (external binary) |
| **UniRef50 (or similar) BLAST database** | PSI-BLAST database |
| `biotite`, `matplotlib` *(optional)* | Cartoon structure plots |


---

## Installation

### From PyPI (recommended)

```bash
pip install s8kpred
```


### With cartoon plot support

```bash
pip install s8kpred[plot]
```


### From GitHub (latest development version)

```bash
pip install git+https://github.com/mayank2801/s8kpred.git
```


### From source

```bash
git clone [https://github.com/mayank2801/s8kpred.git](https://github.com/mayank2801/s8kpred.git)
cd s8kpred
pip install -e .          # editable install
pip install -e .[plot]    # with plotting extras
```


---

## Setting up PSI-BLAST

S8kPred requires NCBI PSI-BLAST to generate evolutionary features. You have two options:

### Option A — System install

```bash
# Ubuntu / Debian
sudo apt install ncbi-blast+

# macOS (Homebrew)
brew install blast

# Conda
conda install -c bioconda -c conda-forge "blast>=2.14"
```


### Option B — Manual download

Download the NCBI BLAST+ toolkit from:
[https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

Then either add the `bin/` folder to your `PATH` or pass the full path via `--psiblast`.

---

## Setting up a BLAST database

S8kPred works best with **UniRef50**. Download and format it:

```bash
# Download
wget [https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz](https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz)
gunzip uniref50.fasta.gz

# Build BLAST database
mkdir -p ~/blast_dbs/uniref50
makeblastdb -in uniref50.fasta \\
            -dbtype prot \\
            -out ~/blast_dbs/uniref50/uniref50 \\
            -title "UniRef50"
```

Then point s8kpred at it:

```bash
export S8KPRED_BLASTDB=~/blast_dbs/uniref50/uniref50
```

or pass `--blastdb ~/blast_dbs/uniref50/uniref50` on every invocation.

---

## Model data files

The trained XGBoost models and lookup tables are **not** bundled in the PyPI wheel because of their size. Download them from the [Releases page](https://github.com/mayank2801/s8kpred/releases) and place them in the `s8kpred/data/` directory inside your Python environment:

```
s8kpred/data/
  TriPeptidePropensityThreeStateSecStructure2AND.csv
  TriPeptidePropensityEightStateSecStructure.csv
  TripeptideBinaryTable_60.csv
  model_3state.json
  model_8state.ubj
```

Or override paths at runtime:

```bash
s8kpred predict -i input.fasta \\
  --blastdb ~/blast_dbs/uniref50/uniref50 \\
  --model-3state /path/to/model_3state.json \\
  --model-8state /path/to/model_8state.ubj
```


---

## Quick start

### Command line

```bash
# Single FASTA file
s8kpred predict -i protein.fasta --blastdb ~/blast_dbs/uniref50/uniref50

# Multi-sequence FASTA
s8kpred predict -i multi_seq.fasta --blastdb ~/blast_dbs/uniref50/uniref50

# Multiple separate FASTA files in one run
s8kpred predict -i seq1.fasta seq2.fasta seq3.fasta \\
                --blastdb ~/blast_dbs/uniref50/uniref50

# Inline sequence (no file needed)
s8kpred predict \\
  --sequence MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD \\
  --id my_protein \\
  --blastdb ~/blast_dbs/uniref50/uniref50

# Custom output folder and job name
s8kpred predict -i input.fasta \\
  --blastdb ~/blast_dbs/uniref50/uniref50 \\
  --output-dir ./results \\
  --job experiment_01

# Skip 8-state prediction
s8kpred predict -i input.fasta --blastdb ... --no-8state

# Skip cartoon plots
s8kpred predict -i input.fasta --blastdb ... --no-plot

# Quiet mode (suppress progress output)
s8kpred predict -i input.fasta --blastdb ... --quiet

# Use more PSI-BLAST threads
s8kpred predict -i input.fasta --blastdb ... --threads 16

# Override PSI-BLAST location
s8kpred predict -i input.fasta \\
  --psiblast /opt/ncbi-blast/bin/psiblast \\
  --blastdb ~/blast_dbs/uniref50/uniref50
```


### Python API

```python
from s8kpred import predict, predict_file

# ── Single sequence ──────────────────────────────────────────────────
result = predict(
    sequence="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD",
    seq_id="my_protein",
    blastdb="/data/blast/uniref50/uniref50",
)

print(result.results_3state["my_protein"])   # e.g. "CCCHHHHHHCCEEEEC..."
print(result.results_8state["my_protein"])   # e.g. "LLLHHHHHHLLEEEELL..."
print(result.job_dir)                         # Path to all output files

# ── Single FASTA file ────────────────────────────────────────────────
result = predict_file(
    fasta_file="proteins.fasta",
    blastdb="/data/blast/uniref50/uniref50",
    output_dir="./results",
)
print(result.summary())

# ── Multi-sequence FASTA ─────────────────────────────────────────────
result = predict_file("multi_seq.fasta", blastdb="...")
for seq_id, ss in result.results_3state.items():
    print(f"{seq_id}: {ss}")

# ── Custom model paths ────────────────────────────────────────────────
from pathlib import Path
result = predict_file(
    "proteins.fasta",
    blastdb="...",
    model_3state=Path("/models/model_3state.json"),
    model_8state=Path("/models/model_8state.ubj"),
)

# ── Skip 8-state to save time ─────────────────────────────────────────
result = predict("MKTAYI...", blastdb="...", run_8state=False)
```


---

## Output files

All outputs are written to a timestamped job directory under `--output-dir`:

```
s8kpred_jobs/
└── 20250210_153042_a1b2c3/
    ├── FASTA/
    │   └── input_sequence.fasta        # combined input
    ├── pssm_outputs/
    │   ├── Seq_1.pssm                  # raw PSI-BLAST PSSM
    │   └── ...
    ├── PSSM_Features_ML_17W.csv        # sliding-window PSSM features
    ├── ResultThreeState.ss2            # PSIPRED-style vertical format
    ├── ResultThreeState.horiz           ├── ResultThreeState.csv            # per-residue probabilities   # PSIPRED-style horizontal format
    ├── ResultThreeState.fas            # pseudo-FASTA format
    ├── ResultThreeState.csv            # per-residue probabilities
    ├── ResultEightState.ss2
    ├── ResultEightState.horiz
    ├── ResultEightState.fas
    ├── ResultEightState.csv
    ├── Seq_1_cartoon.png               # helix/sheet cartoon (requires biotite)
    └── log.dat                         # timing and status log
```


### Secondary structure codes

| Code | State |
| :-- | :-- |
| **3-state** |  |
| H | α-Helix |
| E | β-Strand |
| L | Loop / Coil |
| **8-state** |  |
| H | α-Helix |
| G | 3₁₀-Helix |
| I | π-Helix |
| E | β-Strand |
| B | β-Bridge |
| T | Turn |
| S | Bend |
| L | Loop / Coil |


---

## Environment variables

| Variable | Default | Description |
| :-- | :-- | :-- |
| `S8KPRED_BLASTDB` | *(empty)* | BLAST database path prefix |
| `S8KPRED_PSIBLAST` | `psiblast` | PSI-BLAST binary path |
| `S8KPRED_ITERATIONS` | `3` | PSI-BLAST iterations |


---

## CLI reference

```
s8kpred predict --help
```

```
usage: s8kpred predict [-h] (-i FASTA [FASTA ...] | -s SEQ)
                       [--blastdb DB] [--psiblast BIN] [--iterations N]
                       [--threads N] [-o DIR] [--job ID]
                       [--model-3state PATH] [--model-8state PATH]
                       [--no-3state] [--no-8state] [--no-plot] [-q]
                       [--id ID]
```


---

## Citation

If you use S8kPred in your research, please cite:

> [Your citation here]

---

## License

MIT — see [LICENSE](LICENSE) for details.

Here is a **combined, cleaned, publication‑ready README** saved as a block of plain text that you can copy into a file named `README.txt` (or `README.md`) for your GitHub repo.

You can paste this directly into `README.txt` and commit it.

***

```
# S8kPred — Protein Secondary Structure Prediction

[![PyPI](https://img.shields.io/pypi/v/s8kpred)](https://pypi.org/project/s8kpred/)
[![Python](https://img.shields.io/pypi/pyversions/s8kpred)](https://pypi.org/project/s8kpred/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**S8kPred** predicts protein secondary structure directly from amino acid sequence using XGBoost models trained on PSSM (Position‑Specific Scoring Matrix) and tripeptide propensity features. It supports both fast coarse‑grained 3‑state labels and fine‑grained 8‑state DSSP‑style annotations.

- **3-state**: Helix (H), Beta‑strand (E), Coil/Loop (L)  
- **8-state**: H, G, I, E, B, T, S, L (full DSSP alphabet)

Protein secondary structure prediction is a fundamental problem in structural bioinformatics, providing insights into folding, function, and interaction landscapes. S8kPred combines evolutionary information from PSI‑BLAST with local sequence‑motif features to deliver robust predictions from sequence alone.

---

## 🌐 Web Server

A freely accessible web server is available at:  
👉 **http://s8kpred.in**

- No installation required  
- Accepts FASTA input  
- Returns interactive and downloadable results

---

## 📦 Model & Data Availability

Pretrained XGBoost models and lookup tables are hosted on:

- GitHub **Releases**: https://github.com/mayank2801/s8kpred/releases  
- Zenodo (search for `"S8kPred"`)

These files are excluded from the PyPI wheel due to size and must be downloaded separately.

---

## Requirements

| Dependency | Purpose |
|-----------|---------|
| Python ≥ 3.9 | Runtime environment |
| `numpy`, `pandas`, `xgboost`, `scikit‑learn` | Core machine learning pipeline |
| `biopython` | FASTA input handling |
| **NCBI PSI‑BLAST** | PSSM generation (external binary) |
| **UniRef50 (or similar) BLAST database** | Sequence database for PSI‑BLAST |
| `biotite`, `matplotlib` *(optional)* | Structure cartoons and plotting |

---

## Installation

### From PyPI (recommended)

```bash
pip install s8kpred
```


### With cartoon plot support

```bash
pip install s8kpred[plot]
```


### From GitHub (latest development version)

```bash
pip install git+https://github.com/mayank2801/s8kpred.git
```


### From source (editable)

```bash
git clone https://github.com/mayank2801/s8kpred.git
cd s8kpred
pip install -e .          # editable install
pip install -e .[plot]    # with plotting extras
```


---

## Setting up PSI‑BLAST

S8kPred requires NCBI PSI‑BLAST to generate PSSM‑based evolutionary features.

### Option A — System install

```bash
# Ubuntu / Debian
sudo apt install ncbi‑blast+

# macOS (Homebrew)
brew install blast

# Conda
conda install -c bioconda -c conda‑forge "blast>=2.14"
```


### Option B — Manual download

Download from the NCBI BLAST+ site:
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Then either add the `bin/` folder to your `PATH`, or pass the full path via the `--psiblast` option.

---

## Setting up a BLAST database

For best results, use **UniRef50**.

```bash
# Download
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
gunzip uniref50.fasta.gz

# Build BLAST database
mkdir -p ~/blast_dbs/uniref50
makeblastdb -in uniref50.fasta \
            -dbtype prot \
            -out ~/blast_dbs/uniref50/uniref50 \
            -title "UniRef50"
```

Point S8kPred at the database:

```bash
export S8KPRED_BLASTDB=~/blast_dbs/uniref50/uniref50
```

or pass it at runtime:

```bash
s8kpred predict -i input.fasta --blastdb ~/blast_dbs/uniref50/uniref50
```


---

## Model data files

Download the model artifacts from the **Releases** page and place them in:

```
s8kpred/data/
```

Required files:

- `TriPeptidePropensityThreeStateSecStructure2AND.csv`
- `TriPeptidePropensityEightStateSecStructure.csv`
- `TripeptideBinaryTable_60.csv`
- `model_3state.json`
- `model_8state.ubj`

Alternatively, override paths at runtime:

```bash
s8kpred predict -i input.fasta \
  --blastdb ~/blast_dbs/uniref50/uniref50 \
  --model-3state /path/to/model_3state.json \
  --model-8state /path/to/model_8state.ubj
```


---

## Quick start

### Command line

Single FASTA file:

```bash
s8kpred predict -i protein.fasta --blastdb ~/blast_dbs/uniref50/uniref50
```

Multiple sequences in one file:

```bash
s8kpred predict -i multi_seq.fasta --blastdb ~/blast_dbs/uniref50/uniref50
```

Multiple separate files:

```bash
s8kpred predict -i seq1.fasta seq2.fasta seq3.fasta \
                --blastdb ~/blast_dbs/uniref50/uniref50
```

Inline sequence (no file):

```bash
s8kpred predict \
  --sequence MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD \
  --id my_protein \
  --blastdb ~/blast_dbs/uniref50/uniref50
```

Custom output folder and job name:

```bash
s8kpred predict -i input.fasta \
  --blastdb ~/blast_dbs/uniref50/uniref50 \
  --output-dir ./results \
  --job experiment_01
```

Skip 8‑state prediction:

```bash
s8kpred predict -i input.fasta --blastdb ... --no-8state
```

Skip cartoon plots:

```bash
s8kpred predict -i input.fasta --blastdb ... --no-plot
```

Quiet mode:

```bash
s8kpred predict -i input.fasta --blastdb ... --quiet
```

Use more threads:

```bash
s8kpred predict -i input.fasta --blastdb ... --threads 16
```

Override PSI‑BLAST binary:

```bash
s8kpred predict -i input.fasta \
  --psiblast /opt/ncbi‑blast/bin/psiblast \
  --blastdb ~/blast_dbs/uniref50/uniref50
```


### Python API

Single sequence:

```python
from s8kpred import predict

result = predict(
    sequence="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD",
    seq_id="my_protein",
    blastdb="/data/blast/uniref50/uniref50",
)

print(result.results_3state["my_protein"])   # e.g. "CCCHHHHHHCCEEEEC..."
print(result.results_8state["my_protein"])   # e.g. "LLLHHHHHHLLEEEELL..."
print(result.job_dir)                         # Path to all output files
```

Single FASTA file:

```python
from s8kpred import predict_file

result = predict_file(
    fasta_file="proteins.fasta",
    blastdb="/data/blast/uniref50/uniref50",
    output_dir="./results",
)
print(result.summary())
```

Multi‑sequence FASTA:

```python
result = predict_file("multi_seq.fasta", blastdb="...")
for seq_id, ss in result.results_3state.items():
    print(f"{seq_id}: {ss}")
```

Custom model paths:

```python
from pathlib import Path
result = predict_file(
    "proteins.fasta",
    blastdb="...",
    model_3state=Path("/models/model_3state.json"),
    model_8state=Path("/models/model_8state.ubj"),
)
```

Skip 8‑state:

```python
result = predict("MKTAYI...", blastdb="...", run_8state=False)
```


---

## Output files

All outputs are written to a timestamped job directory under `--output-dir`:

```
s8kpred_jobs/
└── 20250210_153042_a1b2c3/
    ├── FASTA/
    │   └── input_sequence.fasta        # combined input
    ├── pssm_outputs/
    │   ├── Seq_1.pssm                  # raw PSI‑BLAST PSSM
    │   └── ...
    ├── PSSM_Features_ML_17W.csv        # sliding‑window PSSM features
    ├── ResultThreeState.ss2            # 3‑state, PSIPRED‑style vertical
    ├── ResultThreeState.horiz          # 3‑state, PSIPRED‑style horizontal
    ├── ResultThreeState.fas            # pseudo‑FASTA format
    ├── ResultThreeState.csv            # per‑residue 3‑state probabilities
    ├── ResultEightState.ss2            # 8‑state, vertical
    ├── ResultEightState.horiz          # 8‑state, horizontal
    ├── ResultEightState.fas            # pseudo‑FASTA format
    ├── ResultEightState.csv            # per‑residue 8‑state probabilities
    ├── Seq_1_cartoon.png               # helix/sheet cartoon (requires biotite)
    └── log.dat                         # timing and status log
```


### Secondary structure codes

**3‑state:**


| Code | State |
| :-- | :-- |
| H | α‑Helix |
| E | β‑Strand |
| L | Loop / Coil |

**8‑state (DSSP):**


| Code | State |
| :-- | :-- |
| H | α‑Helix |
| G | 3₁₀‑Helix |
| I | π‑Helix |
| E | β‑Strand |
| B | β‑Bridge |
| T | Turn |
| S | Bend |
| L | Loop / Coil |


---

## Environment variables

| Variable | Default | Description |
| :-- | :-- | :-- |
| `S8KPRED_BLASTDB` | *(empty)* | BLAST database path prefix |
| `S8KPRED_PSIBLAST` | `psiblast` | PSI‑BLAST binary path |
| `S8KPRED_ITERATIONS` | `3` | PSI‑BLAST iterations |


---

## CLI reference

```bash
s8kpred predict --help
```

Typical usage:

```
usage: s8kpred predict [-h] (-i FASTA [FASTA ...] | -s SEQ)
                       [--blastdb DB] [--psiblast BIN] [--iterations N]
                       [--threads N] [-o DIR] [--job ID]
                       [--model-3state PATH] [--model-8state PATH]
                       [--no-3state] [--no-8state] [--no-plot] [-q]
                       [--id ID]
```


---

## 📖 Citation

If you use **S8kPred** in your research, please cite:

> Kumar, M., & Rathore, R.S. (2026). S8kPred : A Novel Approach for Protein Secondary Structure Prediction Using 8000 Tripeptide Propensities. Peptide Science.
> DOI: **10.1002/pep2.70029**

(You can expand this line with bibtex or full journal details as needed.)

---

## 🔬 Reproducibility

- **Code**: https://github.com/mayank2801/s8kpred
- **Models**: GitHub Releases \& Zenodo (search `"S8kPred"`)
- **Web server**: http://s8kpred.in

---

## 📄 License

This project is licensed under the **MIT License** – see the `LICENSE` file for details.

```

