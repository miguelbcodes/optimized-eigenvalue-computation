# Optimized Eigenvalue Computation for Square Matrices

This repository contains a practical project for the **Ordinary Differential Equations** course at the **Federal University of Goiás (UFG)**.

The goal is to implement and analyze numerical methods for computing eigenvalues of square matrices, including scalability tests for execution time and memory usage.

For the full project description, see the PDF: [Project Statement](./project-statement.pdf).

## Build

```bash
make            # builds the main program (compute-eigenvalues)
make benchmark  # builds the benchmark program
```

## Running the Benchmark

The benchmark outputs CSV data to stdout. To collect results across multiple runs:

```bash
# First run — creates the file with the CSV header
./benchmark > benchmark_results.csv

# Subsequent runs — appends only the data rows (skips the header)
./benchmark | tail -n +2 >> benchmark_results.csv
```

`tail -n +2` starts output from line 2 onwards, skipping the header so it doesn't appear duplicated in the middle of the CSV.

## Generating Scalability Plots

```bash
pip install -r requirements.txt
python3 scalability_analysis.py benchmark_results.csv
```

This produces `scalability_time.png` and `scalability_memory.png` with scatter plots, power-law fits, and error bars (when multiple samples per N are available).