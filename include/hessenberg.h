#ifndef HESSENBERG_H
#define HESSENBERG_H

#include "matrix.h"

// Reduce a dense matrix A to upper Hessenberg form H via
// similarity transformation H = Q^T A Q using Householder reflectors.
Matrix reduce_to_hessenberg(Matrix A);

#endif // HESSENBERG_H
