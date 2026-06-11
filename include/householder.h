#ifndef HOUSEHOLDER_H
#define HOUSEHOLDER_H

#include "vector.h"
#include "matrix.h"

// Compute a unit Householder reflector v such that (I - 2*v*v^T)*x = -sign(x0)*||x||*e1.
Vector householder_vector(const Vector& x);

// Left-multiply a subregion of A by the Householder reflector (I - 2*v*v^T).
// Acts on rows [row_start, row_start + v.size()) and columns [col_start, col_end).
void apply_householder_left(Matrix& A, const Vector& v,
                            std::size_t row_start, std::size_t col_start, std::size_t col_end);

// Right-multiply a subregion of A by the Householder reflector (I - 2*v*v^T).
// Acts on rows [row_start, row_end) and columns [col_start, col_start + v.size()).
void apply_householder_right(Matrix& A, const Vector& v,
                             std::size_t row_start, std::size_t row_end, std::size_t col_start);

#endif // HOUSEHOLDER_H
