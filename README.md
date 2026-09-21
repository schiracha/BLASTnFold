# AF3 + USalign Structural Modelling Pipeline

**Project:** Evolutionary origins of caveolin proteins beyond Metazoa  
**PI:** Anne K. Kenworthy — University of Virginia, Dept. of Molecular Physiology & Biological Physics  
**Cluster:** UVA Rivanna/Afton HPC (SLURM; account `cavorigins`)  
**AF3 version:** AlphaFold3 v3.0.0  

---

## Overview

This pipeline takes sequence hits from upstream DELTA-BLAST searches, predicts
homo-oligomeric assemblies with AlphaFold3, optionally embeds the protein in a
lipid environment for improved membrane protein modelling, and evaluates
structural similarity to caveolin reference structures via USalign TM-scores.

```
Input (FASTA / region list / accession list / SummaryHits CSV)
            ↓
  make_json_unified.sbatch    → 1-mer JSON per sequence  +  manifest.tsv
            ↓
     AF3_data.sbatch          → MSA/template search; fan-out to N-mer data JSONs
                                token budget check; lipid injection
            ↓
   AF3_inference.sbatch       → GPU structure prediction
                                2 seeds × 5 diffusion samples = 10 models per input
            ↓
      US_align.sbatch         → TM-score vs 6 reference structures
                                skipped CIFs logged to Align_Out/skipped.log
            ↓
       summarize.sh           → per-reference CSVs + SummaryHits.csv
                                skipped CIF report at end of run
```

`submit_af3.sh` handles array sizing, SLURM dependency chaining, and all
tunable parameters from a single command.

---

## Directory Layout

All scripts expect to run from a single **working directory**:

```
<working_dir>/
├── AF3_Inputs/               1-mer JSON files + manifest.tsv
├── AF3_Data/                 *_data.json fan-out files (flat)
├── AF3_Inference/            predicted CIF files
│   └── <basename>/
│       ├── <basename>_0.cif .. _4.cif        (single-seed: sample index only)
│       └── <basename>_seed1_sample0.cif ...  (multi-seed: seed + sample)
├── Align_Out/                USalign stats files + skipped.log
├── cluster_out/              CD-HIT output (from upstream ClusterHits.sbatch)
├── blast_out/                BLAST output (from upstream DBLAST pipeline)
├── logs/                     SLURM stdout/stderr
├── refs/                     download sentinel (.downloads_complete)
├── 7SC0.cif                  reference structures (auto-downloaded)
├── 9DN1.cif
├── 7D60.cif
├── 7RSL.cif
├── AF-Q03135.cif
└── AF-Complex.cif            ← place manually before running US_align
```

---

## Scripts

### `make_json_unified.sbatch`

Unified input preprocessor. Fetches sequences (or reads from FASTA) and writes
one 1-mer AF3 JSON per entry, plus an append-only TSV provenance manifest.

**Input formats (auto-detected, or override with `--mode`):**

| Mode | Format | Typical source |
|------|--------|----------------|
| `fasta` | Standard FASTA (`.fasta`, `.fa`, `.afa`) | `ClusterHits.sbatch` output |
| `region` | `ACCESSION:START-END` or `ACC\tSTART\tEND` | BLAST `*_cd.out` files |
| `accession` | One accession per line | Plain accession lists |
| `hits` | `SummaryHits.csv` from `summarize.sh` | Iterative re-modelling |

**Usage:**
```bash
sbatch make_json_unified.sbatch INPUT_FILE [OPTIONS]

Options:
  --mode  fasta|region|accession|hits   override auto-detection
  --outdir DIR                          output directory (default: ./AF3_Inputs)
  --min-tm-avg  FLOAT                   hits mode: min TMscore_avgL (default 0.2)
  --min-tm-refl FLOAT                   hits mode: min TMscore_refL (default 0.0)
  --source-tag  TAG                     provenance label in manifest
  --skip-existing                       skip JSONs that already exist (resume)
```

**Key behaviours:**
- FASTA mode reads sequences directly — no BLAST DB fetch needed
- Gap characters (`-`) in aligned FASTA are stripped automatically
- Version suffixes (`.1`, `.2`) stripped from accessions for consistent naming
- `region` mode parser splits on `:`, `\t`, and `-`; NCBI accessions (WP_, XP_,
  NP_) are safe; UniProt accessions with hyphens should use `--mode accession`
- Parallel workers write to per-PID temp files, merged after completion to avoid
  manifest corruption
- `modelSeeds: [1, 2]` — two seeds per model by default

**Manifest columns:**
`json_basename`, `accession`, `region`, `seq_len`, `source_tag`, `input_file`, `timestamp`

---

### `AF3_data.sbatch`

Runs the AF3 data pipeline (MSA + template search) for one 1-mer JSON per array
task, then fans out the resulting `*_data.json` into one copy per oligomer size.
Injects a lipid environment into each fan-out JSON.

**Configuration block (edit at top of script):**
```bash
CHAIN_MIN=5          # minimum oligomer size to model
CHAIN_MAX=15         # maximum oligomer size to model
LIPID_TOKENS=258     # lipid token budget (0 = disabled)
```

**Lipid environment:**
Three lipids from the PDB CCD database are used in equal proportion:

| CCD | Name | Tokens/molecule |
|-----|------|----------------|
| `PLW` | POPC | 52 |
| `EPE` | POPE | 49 |
| `CLR` | Cholesterol | 28 |

One set = 1 POPC + 1 POPE + 1 CLR = **129 tokens**.
`LIPID_TOKENS=258` → 2 sets (258 tokens). Set to `0` to disable entirely.

**Token budget check:**
The fan-out Python computes `seq_len × n_chains + lipid_tokens` and caps
`chain_max` at the largest oligomer that fits within AF3's 5000-token limit.
If even a monomer + lipids exceeds the limit, lipids are disabled for that
sequence with a warning. If a bare monomer exceeds 5000 tokens, the task exits
with an error.

**Resources:** 16 CPUs, 32G RAM per task (standard partition).

**Output:** `AF3_Data/<acc>[_START-END]_<N>mer_data.json` for N in [CHAIN_MIN, CHAIN_MAX]

---

### `AF3_inference.sbatch`

GPU structure prediction for one `*_data.json` per array task.

**Key settings:**
- Partition: `gpu`, requires A6000 or A100 (`-C "a6000|a100"`)
- Default: **2 seeds × 5 diffusion samples = 10 models per input**
- Seed count detected automatically from `modelSeeds` array in the JSON

**Memory tuning (automatic):**

| Oligomer size | `XLA_CLIENT_MEM_FRACTION` |
|---------------|--------------------------|
| ≤ 8-mer | 2.5 |
| 9–11-mer | 3.2 |
| ≥ 12-mer | 3.5 |

**Output naming:**

| Condition | Convention | Example |
|-----------|-----------|---------|
| Single seed | `<basename>_<sample>.cif` | `wp_123_11mer_0.cif` |
| Multi-seed | `<basename>_seed<S>_sample<M>.cif` | `wp_123_11mer_seed1_sample0.cif` |

**AF3 version notes:**
Scripts target **AF3 3.0.0**. `PYTHONNOUSERSITE=1` is set in both data and
inference scripts to prevent user-installed AF3 packages in `~/.local` from
shadowing the module — critical when switching between 3.0.0 and 3.0.1.

Flags that exist in **3.0.1 only** (do not add these for 3.0.0):

| Flag | Added in |
|------|----------|
| `--num_diffusion_samples` | 3.0.1 |
| `--diffusion_num_samples` | 3.0.1 |
| `--num_recycles` | 3.0.1 |
| `--num_seeds` | 3.0.1 |
| `--max_template_date` | 3.0.1 |
| `--conformer_max_iterations` | 3.0.1 |

---

### `US_align.sbatch`

Structural alignment of one model CIF per array task against six reference
structures. Uses `flock`-based sentinel for reference downloads (no race
conditions across parallel tasks).

**References:**

| Key | File | Role | `-mm 1` |
|-----|------|------|---------|
| `7SC0` | `7SC0.cif` | Positive control | ✓ |
| `9DN1` | `9DN1.cif` | Positive control | ✓ |
| `7D60` | `7D60.cif` | Negative control | ✓ |
| `7RSL` | `7RSL.cif` | Negative control | ✓ |
| `AF-Q03135` | `AF-Q03135.cif` | Neutral — AF2 human Cav-1 monomer | ✗ |
| `AFComp` | `AF-Complex.cif` | Neutral — AF-predicted complex | ✓ |

7SC0, 9DN1, 7D60, 7RSL, and AF-Q03135 are auto-downloaded by task 1.
**`AF-Complex.cif` must be placed manually** in the working directory.

**Skipped CIF logging:**
If the input CIF is missing (inference task failed), the task logs a line to
`Align_Out/skipped.log` and exits cleanly (code 0). The align array uses
`--dependency=afterany` on the inference array so it runs even if some
inference tasks failed.

`skipped.log` format: `TIMESTAMP\tMISSING_CIF\tTASK_ID\tCIF_PATH`

---

### `summarize.sh`

Parses all USalign `_stats.out` files and produces structured CSV summaries.

**Usage:**
```bash
bash summarize.sh [Align_Out_dir]
# Default: ./Align_Out relative to CWD
```

**Per-reference CSVs** (`<REF>_<PARDIR>.csv`):

| Column | Description |
|--------|-------------|
| `basename` | Full CIF basename |
| `accession` | Parsed from basename |
| `region` | Coordinates (`START-END`) or `full` |
| `oligomer_size` | Oligomeric state |
| `sample` | Diffusion sample index |
| `aligned_len` | USalign aligned residue count |
| `RMSD` | Root-mean-square deviation (Å) |
| `Seq_ID` | Sequence identity of aligned region |
| `TMscore_refL` | TM-score normalised to reference length |
| `TMscore_L1` | TM-score normalised to query length |
| `TMscore_L2` | TM-score normalised to reference length |
| `TMscore_avgL` | TM-score normalised to average length |

**SummaryHits CSV** (`<PARDIR>_SummaryHits.csv`):
Filter threshold `TMscore_avgL > 0.2`. Columns: `accession`, `region`,
`oligomer_size`, `sample`, `reference`, `TMscore_refL`, `TMscore_avgL`.
This file is the direct input for the next pipeline iteration (hits mode).

**TM-score thresholds:**

| Threshold | Interpretation |
|-----------|---------------|
| > 0.5 | Definitive structural similarity (same fold) |
| > 0.4 | Probable structural similarity |
| > 0.2 | Possible — screening filter in SummaryHits |

At the end of the run, any entries in `Align_Out/skipped.log` are reported
with instructions for re-running failed inference tasks.

---

### `patch_add_lipids.sbatch`

Adds (or removes) lipid environment entries to all existing `*_data.json` files
in `AF3_Data/`. Use to retrofit lipids onto a data pipeline run completed
without them.

**Usage:**
```bash
sbatch patch_add_lipids.sbatch --dry-run          # preview
sbatch patch_add_lipids.sbatch                    # apply default (258 tokens)
sbatch patch_add_lipids.sbatch --lipid-tokens 387 # custom budget
sbatch patch_add_lipids.sbatch --lipid-tokens 0   # remove lipids
sbatch patch_add_lipids.sbatch --data-dir /path/to/AF3_Data
```

Idempotent: files that already contain the correct lipid entries are skipped.
Files with wrong/old ligand entries are replaced.

---

### `submit_af3.sh`

Master orchestration script. Run directly on the login node (`bash submit_af3.sh`).
Do **not** `sbatch` this script.

**Full option reference:**

| Option | Default | Description |
|--------|---------|-------------|
| `--from STAGE` | `make_json` | Resume from: `make_json\|data\|inference\|align\|summarize` |
| `--mode MODE` | auto | Input flavour: `fasta\|region\|accession\|hits` |
| `--source-tag TAG` | input filename | Provenance label written to manifest |
| `--chain-min N` | `5` | Minimum oligomer size |
| `--chain-max N` | `15` | Maximum oligomer size |
| `--min-tm-avg F` | `0.2` | Hits mode: min TMscore_avgL filter |
| `--lipid-tokens N` | `258` | Lipid token budget (0 = disabled) |
| `--throttle-data N` | `16` | Max concurrent data array tasks |
| `--throttle-infer N` | `4` | Max concurrent inference array tasks |
| `--throttle-align N` | `32` | Max concurrent align array tasks |
| `--chain-estimate` | off | Size arrays from input file line count |
| `--skip-existing` | off | Pass through to make_json |
| `--dry-run` | off | Print sbatch commands without submitting |

**Array sizing rules:**

| `--from` | data sized from | inference sized from | align sized from |
|---|---|---|---|
| `make_json` | stops; prints re-run command | — | — |
| `data` | `AF3_Inputs/` disk | derived: N × oligomers | derived: N × 10 |
| `inference` | `AF3_Inputs/` disk | `AF3_Data/` disk | derived: N × 10 |
| `align` | `AF3_Inputs/` disk | `AF3_Data/` disk | `AF3_Inference/` disk |

**Dependency chain:**
```
make_json  ──afterok──►  data array  ──afterok──►  inference array  ──afterany──►  align array  ──afterok──►  summarize
```
`afterany` on the align stage means the align array starts even if some
inference tasks failed; missing CIFs are logged and skipped rather than
blocking the entire alignment run.

---

## Typical Workflows

### Workflow A — From clustered FASTA (recommended starting point)
```bash

sbatch make_json_unified.sbatch clustered_sequences.fasta --source-tag Run01

# After make_json finishes:
bash submit_af3.sh --from data --chain-min 11 --chain-max 11 --source-tag Run01
```

### Workflow B — From BLAST region list (two-step)
```bash
bash submit_af3.sh blast_out/ACC000_nonmet_consolidated_cd.tsv \
    --mode region --chain-min 11 --chain-max 11 --source-tag ACC000
# Script submits make_json then prints the --from data command to run next.
```

### Workflow C — One-shot with chain estimate
```bash
bash submit_af3.sh accessions.txt --mode accession \
    --chain-min 11 --chain-max 11 --chain-estimate --source-tag Run01 \
    --throttle-data 16 --throttle-infer 4
```

### Workflow D — Second-round modelling of top structural hits
```bash
bash submit_af3.sh Run01_SummaryHits.csv --mode hits \
    --min-tm-avg 0.3 --chain-min 5 --chain-max 15 --source-tag Run02_hits
```

### Workflow E — Add lipids to existing AF3_Data/ before inference
```bash
sbatch patch_add_lipids.sbatch --lipid-tokens 258
bash submit_af3.sh --from inference --chain-min 11 --chain-max 11
```

### Workflow F — Resume after partial inference failure
```bash
cat Align_Out/skipped.log          # see what failed
bash submit_af3.sh --from align --chain-min 11 --chain-max 11
# Completed alignments are skipped (idempotent)
```
