#include "../include/tridiagonal.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
    int n = 100;
    double tolerance = 1e-10;

    for (int i = 1; i < argc - 1; ++i) {
        if (std::strcmp(argv[i], "--n") == 0)
            n = std::atoi(argv[i + 1]);
        else if (std::strcmp(argv[i], "--epsilon") == 0)
            tolerance = std::atof(argv[i + 1]);
    }

    // §4.1 — Build the CSTR tridiagonal matrix directly as bands.
    // Each reactor exchanges mass only with its immediate neighbours:
    //   diagonal = -2.5 (degradation + outflow), off-diagonal = 1.0 (inter-reactor flow).
    SymmetricTridiagonal A;
    A.diag.assign(n, -2.5);
    A.offdiag.assign(n - 1, 1.0);

    std::cout << "=== Fase 4: Aplicacao Pratica — Rede de Reatores CSTR ===\n\n";
    std::cout << "Numero de reatores (N): " << n << "\n";
    std::cout << "Matriz tridiagonal simetrica: diagonal = -2.5, sub/super-diagonal = 1.0\n";
    std::cout << "Tolerancia: " << std::scientific << std::setprecision(1) << tolerance
              << std::defaultfloat << "\n\n";

    // §4.2 — Eigenvalues via the Phase 3 optimized path (O(N²), O(N) memory).
    std::vector<double> eigenvalues = qr_iterate_tridiagonal(A, tolerance);
    std::sort(eigenvalues.begin(), eigenvalues.end());

    std::cout << "--- Autovalores (" << n << " valores, ordem crescente) ---\n";
    std::cout << std::fixed << std::setprecision(6);
    for (int i = 0; i < static_cast<int>(eigenvalues.size()); ++i) {
        std::cout << std::setw(12) << eigenvalues[i];
        if ((i + 1) % 5 == 0)
            std::cout << '\n';
    }
    if (eigenvalues.size() % 5 != 0)
        std::cout << '\n';

    // §4.2 (validation) — Closed-form spectrum: λ_k = -2.5 + 2·cos(kπ/(N+1)).
    std::vector<double> analytical(n);
    for (int k = 1; k <= n; ++k)
        analytical[k - 1] = -2.5 + 2.0 * std::cos(k * M_PI / (n + 1));
    std::sort(analytical.begin(), analytical.end());

    double max_diff = 0.0;
    for (int i = 0; i < n; ++i)
        max_diff = std::max(max_diff, std::abs(eigenvalues[i] - analytical[i]));

    std::cout << std::scientific << std::setprecision(2);
    std::cout << "\nValidacao (formula analitica): max |lambda_num - lambda_analitico| = "
              << max_diff << "\n";

    // §4.3 — Stability: all eigenvalues real; stable iff max(λ) < 0.
    double lambda_max = eigenvalues.back();
    double lambda_min = eigenvalues.front();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n--- Estabilidade ---\n";
    std::cout << "Maior autovalor (mais proximo de 0): " << lambda_max << "\n";
    std::cout << "Menor autovalor (mais negativo):     " << lambda_min << "\n";
    std::cout << "Veredicto: " << (lambda_max < 0 ? "ESTAVEL" : "INSTAVEL")
              << " (todos os autovalores " << (lambda_max < 0 ? "negativos" : "nao sao todos negativos")
              << ")\n";

    // §4.4 — Stiffness ratio: S = max|λ| / min|λ|.
    double max_abs = std::max(std::abs(lambda_max), std::abs(lambda_min));
    double min_abs = std::min(std::abs(lambda_max), std::abs(lambda_min));

    std::cout << "\n--- Razao de Rigidez ---\n";
    if (min_abs > 0) {
        double S = max_abs / min_abs;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "S = max|lambda| / min|lambda| = " << max_abs << " / " << min_abs
                  << " = " << S << "\n";
        bool stiff = S > 1e3;
        std::cout << "Classificacao: " << (stiff ? "RIGIDO (S > 1000)" : "NAO-RIGIDO (S < 1000)")
                  << "\n";
        std::cout << "Recomendacao: solver " << (stiff ? "IMPLICITO" : "EXPLICITO (ex: RK4)")
                  << "\n";
    } else {
        std::cout << "min|lambda| = 0 — razao de rigidez indefinida.\n";
    }

    return 0;
}
