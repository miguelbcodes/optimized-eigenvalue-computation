#include "../include/matrix.h"
#include "../include/hessenberg.h"
#include "../include/qr.h"
#include <cstdlib>
#include <cstring>

int main(int argc, char* argv[]) {
    double tolerance = 1e-10;
    for (int i = 1; i < argc - 1; ++i) {
        if (std::strcmp(argv[i], "--epsilon") == 0)
            tolerance = std::atof(argv[i + 1]);
    }

    // Phase 1.1: Redução à Forma de Hessenberg (Householder)
    // Phase 1.2: Iteração QR até convergência da subdiagonal
    Matrix A(5);
    Matrix::fill_random(A);
    A.print();

    std::cout << "\nForma de Hessenberg:\n";
    Matrix H = reduce_to_hessenberg(A);
    H.print();

    std::cout << "\nAutovalores (diagonal apos QR):\n";
    Matrix schur = qr_iterate(H, tolerance);
    schur.print();

    // Phase 2: Estudo de Escalabilidade — ver bench/benchmark.cpp

    return 0;
}
