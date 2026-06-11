#ifndef QR_H
#define QR_H

#include "matrix.h"
#include "vector.h"

struct QRResult {
    Matrix Q;
    Matrix R;
};

// QR factorization of an upper Hessenberg matrix using Givens rotations.
QRResult qr_factorize_hessenberg(const Matrix& H);

// Iterative QR algorithm with Wilkinson shift and deflation.
// Repeatedly applies H = QR, H <- RQ until the subdiagonal converges to zero.
// Returns the real Schur form: real eigenvalues on the diagonal,
// complex conjugate pairs as 2x2 blocks.
Matrix qr_iterate(Matrix H, double tolerance, int max_iterations = 10000);

#endif // QR_H
