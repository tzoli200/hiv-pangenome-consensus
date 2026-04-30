# How to run the pipeline

This walks through running the pipeline on the bundled example dataset (HIV-1 SIM_AE replicate 10). The same commands work on your own data, just swap the filenames.

## Setup

Place all input files in your working directory. For the bundled example:

```bash
cp example/inputs/* .
```

## Map reads to the pangenome graph

```bash
docker run --rm -v $(pwd):/data quay.io/biocontainers/graphaligner:1.0.17b--h21ec9f0_2 \
    GraphAligner -g /data/pggb_optimized_16ref.gfa \
    -f /data/HIV1_CON_AE_simref_QSsim_R1_paired.fastq \
    -f /data/HIV1_CON_AE_simref_QSsim_R2_paired.fastq \
    -a /data/aligned.gaf \
    -x vg
```

## Build correspondence table

```bash
python3 scripts/step1_correspondence.py \
    --msa HIV1_CON_2021_optimized_16_pansn_corrected_gapped.fasta \
    --gfa pggb_optimized_16ref.gfa \
    --output correspondence.tsv
```

## Strip `#0#genome` suffix (required, dont skip it)

Step 1 emits PanSN-tagged sequence names; step 2 expects them stripped. Don't skip this, without it, step 2 will report `Mapped X nodes with 0 positions` and produce an empty consensus.

```bash
cp correspondence.tsv correspondence_original.tsv
sed 's/#0#genome//g' correspondence_original.tsv > correspondence.tsv
```

## Generate consensus

```bash
python3 scripts/step2_consensus.py \
    --gfa pggb_optimized_16ref.gfa \
    --gaf aligned.gaf \
    --fastq1 HIV1_CON_AE_simref_QSsim_R1_paired.fastq \
    --fastq2 HIV1_CON_AE_simref_QSsim_R2_paired.fastq \
    --correspondence correspondence.tsv \
    --output-consensus consensus.fasta \
    --output-name SIM_AE_rep10
```

## Mismatch analysis

Optional, used for benchmarking against a known ground truth and a Shiver consensus.

```bash
python3 scripts/mismatch_analysis.py \
    --v5 consensus.fasta \
    --shiver shiver_genome.fasta \
    --gt HIV1_CON_AE_simref_QSsim_aligned_consensus_fixed.fasta \
    --label SIM_AE_rep10 \
    --output-dir .
```

## Running on your own data

Replace the example filenames with yours. The required inputs are:

- A PGGB-built variation graph (`.gfa`) with PanSN-tagged paths (sequences ending in `#0#genome`)
- The gapped MSA used to build the graph, with matching PanSN-tagged headers
- Paired-end FASTQ reads for one sample
- (Optional) A ground truth FASTA and a Shiver consensus, if you want to run Stage 3
