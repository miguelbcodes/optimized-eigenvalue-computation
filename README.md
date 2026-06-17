# Optimized Eigenvalue Computation for Square Matrices

This repository contains a practical project for the **Ordinary Differential Equations** course at the **Federal University of Goiás (UFG)**.

The goal is to implement and analyze numerical methods for computing eigenvalues of square matrices, including scalability tests for execution time and memory usage.

For the full project description, see the PDF: [Project Statement](./project-statement.pdf).

## Build

```bash
make            # builds the main program (build/compute-eigenvalues)
make benchmark  # builds the benchmark program (build/benchmark)
```

## Running the Main Program

```bash
./build/compute-eigenvalues             # demo with default tolerance
./build/compute-eigenvalues --epsilon 1e-12
```

## Running the Benchmark

The benchmark outputs CSV data to stdout. To collect results across multiple runs:

```bash
# First run — creates the file with the CSV header
./build/benchmark > benchmark_results.csv

# Subsequent runs — appends only the data rows (skips the header)
./build/benchmark | tail -n +2 >> benchmark_results.csv
```

`tail -n +2` starts output from line 2 onwards, skipping the header so it doesn't appear duplicated in the middle of the CSV.

## Generating Scalability Plots

```bash
pip install -r bench/requirements.txt
python3 bench/scalability_analysis.py benchmark_results.csv
```

This produces `scalability_time.png` and `scalability_memory.png` with scatter plots, power-law fits, and error bars (when multiple samples per N are available). When the CSV contains both algorithms, the two series are overlaid on the same axes and a per-N **percentage performance gain** table is printed to the console.

## Structural Optimization (Phase 3)

For **symmetric** matrices (`A = Aᵀ`) the project ships a dedicated, faster path:

- `Matrix::fill_random_symmetric` generates random symmetric inputs.
- `reduce_to_tridiagonal` applies Householder reflectors with a symmetry-exploiting rank-2 update; on a symmetric matrix the Hessenberg form collapses to a **symmetric tridiagonal** matrix, stored as two `O(N)` bands instead of a dense `O(N²)` matrix.
- `qr_iterate_tridiagonal` runs an implicit-shift QL iteration directly on those bands in `O(N²)`, versus `O(N³)` for the general dense Hessenberg QR. Because the matrix is symmetric, all eigenvalues are real.

Running `./build/compute-eigenvalues` prints a demonstration (N=5) that the general Householder reduction applied to a symmetric matrix collapses to symmetric tridiagonal form, then computes the eigenvalues through the optimized path.

The benchmark runs **both** the general and the optimized algorithm on the same symmetric matrix for every `N`, emitting `hessenberg_qr` and `tridiagonal_qr` rows so the plots and gain table can compare them directly.

## Practical Application — Reactor Network (Phase 4)

Implements §5 of the [Project Statement](./project-statement.pdf): eigenvalue analysis of a CSTR reactor network for stability and stiffness classification.

```bash
make reactor                            # builds build/reactor-network
./build/reactor-network                 # default: N=100, tolerance=1e-10
./build/reactor-network --n 50 --epsilon 1e-12
```
