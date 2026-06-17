#ifndef TRIDIAGONAL_H
#define TRIDIAGONAL_H

#include "matrix.h"
#include <vector>

// Compact representation of a symmetric tridiagonal matrix:
//   diag[i]    = T(i, i)               for i in [0, n)
//   offdiag[i] = T(i+1, i) = T(i, i+1) for i in [0, n-1)
// Stores only the two non-trivial bands, so memory is O(N) instead of O(N^2).
struct SymmetricTridiagonal {
    std::vector<double> diag;
    std::vector<double> offdiag;
};

// Reduce a symmetric matrix to symmetric tridiagonal form via Householder
// reflectors. Exploits symmetry with a rank-2 update (B <- B - w*v^T - v*w^T),
// touching only the trailing block at each step, and returns just the two bands.
SymmetricTridiagonal reduce_to_tridiagonal(Matrix A);

// Compute all eigenvalues of a symmetric tridiagonal matrix using the
// implicit-shift QL algorithm with Wilkinson shift and deflation. Costs O(N^2),
// versus O(N^3) for the general dense Hessenberg QR. Because the matrix is
// symmetric, every eigenvalue is real. The returned values are unordered.
std::vector<double> qr_iterate_tridiagonal(SymmetricTridiagonal T,
                                           double tolerance,
                                           int max_iterations = 10000);

#endif // TRIDIAGONAL_H
