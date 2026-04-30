# HIV Pangenome Consensus Pipeline

A pipeline for generating sample-specific HIV-1 consensus sequences by mapping short reads against a **pangenome variation graph** (PGGB) instead of a single reference. Aiming to avoid the reference bias that affects single-reference workflows on a virus as diverse as HIV-1.

## Pipeline overview

1. **Map** paired-end short reads to a PGGB pangenome graph with GraphAligner
2. **Build** a coordinate translation table between graph positions and MSA columns
3. **Vote** at each MSA column to call a consensus base, with insertion clustering for indel-sliding artefacts
4. **Compare** (optional) the consensus against ground truth and a Shiver consensus

## Scope

This repository covers everything **from a pre-built pangenome graph to a sample-specific consensus** (and an optional benchmarking step). Graph construction itself is upstream of this pipeline and was done with [Panalyze](https://github.com/downingtim/Panalyze); the graph and its accompanying MSA are provided as inputs in `example/inputs/`.

## Repository layout
```
scripts/             Pipeline scripts 
example/inputs/      Demo dataset (HIV-1 SIM_AE simulated reads + pre-built graph)
HOW_TO_RUN.md        Step-by-step run instructions
```

## Requirements

- Python 3 
- Docker (for GraphAligner and MAFFT, pulled automatically on first use)

## Quick start

See **`HOW_TO_RUN.md`** for the full command sequence. The bundled example in `example/inputs/` lets you run the pipeline from read mapping through consensus generation on a known-good dataset.

## Citation

Built with [Panalyze](https://github.com/downingtim/Panalyze) (Downing et al.) for graph construction.
