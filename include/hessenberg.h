#ifndef HESSENBERG_H
#define HESSENBERG_H

#include "matrix.h"

// Reduce a dense matrix to upper Hessenberg form via similarity
// transformation H = Q^T * A * Q using Householder reflectors.
Matrix reduce_to_hessenberg(Matrix A);

#endif // HESSENBERG_H
