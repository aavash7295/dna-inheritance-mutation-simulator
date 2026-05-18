# DNA Inheritance and Mutation Simulator

## Overview

This bioinformatics project simulates DNA inheritance between male and female DNA sequences to generate a simulated child DNA sequence.

The project performs:
- DNA inheritance simulation
- ATCG nucleotide analysis
- GC and AT content analysis
- Mutation simulation
- DNA → RNA transcription
- RNA → Protein translation
- Protein comparison
- Mutation visualization
- Mutation heatmap generation
- FASTA file support

---

## Pipeline

```text
Male DNA + Female DNA
↓
Simulated Child DNA
↓
ATCG Count + GC/AT Content
↓
Mutation Simulation
↓
DNA → RNA
↓
RNA → Protein
↓
Mutation Visualization
```

---

## Technologies Used

- Python
- Biopython
- NumPy
- Matplotlib
- Seaborn

---

## Installation

```bash
pip install -r requirements.txt
```

---

## Running the Project

```bash
python main.py
```

---

## Input Files

```text
male_dna.txt
female_dna.txt
```

---

## Output Files

```text
child_dna_analysis.txt
child_dna.fasta
mutated_child_dna.fasta
```

---

## Future Improvements

- SNP analysis
- Sequence alignment
- Real genome datasets
- Protein similarity scoring
- Machine learning mutation prediction
