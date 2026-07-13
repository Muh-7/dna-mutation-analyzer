# DNA Mutation Analyzer — Project Plan

## Project Goal

Build an AI-powered platform that compares a reference DNA sequence
with a mutated DNA sequence and analyzes the predicted mutation effects.

## Main Inputs

- Reference DNA sequence
- Mutated DNA sequence
- Genomic location when available
- Gene name
- Tissue or cell type

## Main Outputs

- Mutation position and type
- Reference and mutated RNA comparison
- Predicted gene-expression change
- Predicted splicing change
- Transcription-factor binding changes
- Scientific visualizations and final report

## Main Technologies

- AlphaGenome for expression, RNA, splicing, and contextual TF analysis
- JASPAR and FIMO for interpretable TF motif analysis
- Carbon as a later sequence-likelihood extension
- DNABERT-2 as an experimental embedding baseline
- Gradio for the user interface

## MVP Scope

- Human genome
- hg38
- One SNV mutation
- One gene per analysis
- One selected tissue or cell type
- Research and educational use only

## Implementation Order

1. Sequence validation
2. Mutation detection
3. Mutation-window extraction
4. Genomic annotation
5. DNA-to-RNA transcription
6. AlphaGenome expression analysis
7. Splicing analysis
8. JASPAR/FIMO motif analysis
9. Result integration
10. Visualizations
11. Gradio interface
12. Validation
13. Carbon extension
14. DNABERT comparison
15. Deployment

## Current Phase

Version 0.1:

- clean_sequence()
- validate_sequence()
- detect_snv()
- extract_mutation_window()
