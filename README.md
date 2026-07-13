# DNA Mutation Analyzer

An AI-powered research project for analyzing the potential effects of
DNA mutations on RNA, gene expression, splicing, and transcription-factor
binding.

## Project Idea

The system compares two DNA sequences:

- A reference DNA sequence
- A mutated DNA sequence

It identifies the mutation and analyzes its predicted biological effects.

## Planned Analyses

- DNA sequence validation
- Mutation detection
- RNA sequence comparison
- Gene-expression prediction
- Splicing-effect prediction
- Transcription-factor binding analysis
- Scientific visualization
- Final mutation report

## Planned Models and Tools

- AlphaGenome
- JASPAR
- FIMO
- Carbon
- DNABERT-2
- Gradio

## Current Development Phase

The current phase focuses on:

- Cleaning DNA sequences
- Validating DNA sequences
- Detecting a single nucleotide variant
- Extracting a sequence window around the mutation

## Project Structure

```text
src/preprocessing/   DNA validation and mutation detection
notebooks/           Research and experimentation notebooks
tests/               Automated tests
data/examples/       Non-sensitive example sequences
assets/              Project images and figures
```


## Installation
```git
git clone https://github.com/Muh-7/dna-mutation-analyzer.git
cd dna-mutation-analyzer

python -m venv .venv
source .venv/bin/activate
```
pip install -r requirements.txt
