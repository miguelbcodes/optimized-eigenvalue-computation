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

// Iterative QR algorithm: repeatedly factorize H = QR and form H <- RQ
// until subdiagonal elements fall below epsilon.
// Returns the matrix with eigenvalues on the diagonal.
Matrix qr_iterate(Matrix H, double epsilon, int max_iter = 10000);

#endif // QR_H
