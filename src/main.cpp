#include "../include/matrix.h"
#include "../include/hessenberg.h"
#include "../include/qr.h"
#include "../include/tridiagonal.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <vector>

// Largest magnitude among the entries that lie outside the tridiagonal band.
static double max_off_tridiagonal(const Matrix& M) {
    double worst = 0.0;
    std::size_t n = M.size();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            long band = static_cast<long>(i) - static_cast<long>(j);
            if (band < -1 || band > 1)
                worst = std::max(worst, std::abs(M(i, j)));
        }
    return worst;
}

// Largest |M(i,j) - M(j,i)| — a measure of how far the matrix is from symmetric.
static double max_asymmetry(const Matrix& M) {
    double worst = 0.0;
    std::size_t n = M.size();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j)
            worst = std::max(worst, std::abs(M(i, j) - M(j, i)));
    return worst;
}

int main(int argc, char* argv[]) {
    double tolerance = 1e-10;
    for (int i = 1; i < argc - 1; ++i) {
        if (std::strcmp(argv[i], "--epsilon") == 0)
            tolerance = std::atof(argv[i + 1]);
    }

    // === Phase 1: general (non-symmetric) matrix ===
    // Phase 1.1: Redução à Forma de Hessenberg (Householder)
    // Phase 1.2: Iteração QR até convergência da subdiagonal
    std::cout << "=== Fase 1: Matriz generica (N=5) ===\n\n";
    Matrix A(5);
    Matrix::fill_random(A);
    std::cout << "Matriz A:\n";
    A.print();

    std::cout << "\nForma de Hessenberg:\n";
    Matrix H = reduce_to_hessenberg(A);
    H.print();

    std::cout << "\nAutovalores (diagonal apos QR):\n";
    Matrix schur = qr_iterate(H, tolerance);
    schur.print();

    // === Phase 3.2: a symmetric matrix collapses to symmetric tridiagonal ===
    std::cout << "\n=== Fase 3.2: Matriz simetrica (N=5) ===\n\n";
    Matrix S(5);
    Matrix::fill_random_symmetric(S);
    std::cout << "Matriz simetrica S (S = S^T):\n";
    S.print();

    // Apply the *same general* Householder Hessenberg reduction: on a symmetric
    // input it collapses past Hessenberg form into a symmetric tridiagonal.
    std::cout << "\nApos as reflexoes de Householder (mesmo algoritmo geral):\n";
    Matrix Hs = reduce_to_hessenberg(S);
    Hs.print();

    std::cout << std::scientific << std::setprecision(2);
    std::cout << "\nMaior entrada fora da banda tridiagonal: " << max_off_tridiagonal(Hs) << '\n';
    std::cout << "Maior assimetria |H(i,j) - H(j,i)|:      " << max_asymmetry(Hs) << '\n';
    std::cout << std::defaultfloat;
    std::cout << "=> A matriz colapsou para uma Matriz Tridiagonal Simetrica.\n";

    // Eigenvalues via the optimized O(N^2) tridiagonal path.
    std::cout << "\nAutovalores via algoritmo otimizado (tridiagonal, O(N^2)):\n";
    SymmetricTridiagonal T = reduce_to_tridiagonal(S);
    std::vector<double> eigenvalues = qr_iterate_tridiagonal(T, tolerance);
    std::cout << std::fixed << std::setprecision(4);
    for (double lambda : eigenvalues)
        std::cout << std::setw(12) << lambda << '\n';

    // Phase 2 & 3.3/3.4: Estudo de escalabilidade — ver bench/benchmark.cpp

    return 0;
}
