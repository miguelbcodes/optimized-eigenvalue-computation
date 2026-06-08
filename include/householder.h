#ifndef HOUSEHOLDER_H
#define HOUSEHOLDER_H

#include "vector.h"
#include "matrix.h"

// Compute the Householder vector v (unit norm) such that
// (I - 2*v*v^T) * x = -sign(x0)*||x||*e1
Vector householder_vector(const Vector& x);

// Apply reflector (I - 2*v*v^T) from the left to rows [row_start, row_start+v.size())
// and columns [col_start, col_end) of A, in-place.
void apply_householder_left(Matrix& A, const Vector& v,
                            std::size_t row_start, std::size_t col_start, std::size_t col_end);

// Apply reflector (I - 2*v*v^T) from the right to rows [row_start, row_end)
// and columns [col_start, col_start+v.size()) of A, in-place.
void apply_householder_right(Matrix& A, const Vector& v,
                             std::size_t row_start, std::size_t row_end, std::size_t col_start);

#endif // HOUSEHOLDER_H
