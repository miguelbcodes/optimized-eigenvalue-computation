#include "include/householder.h"
#include <cmath>

Vector householder_vector(const Vector& x) {
    std::size_t len = x.size();
    Vector v(len);

    double x_norm = x.norm();
    if (x_norm == 0.0) {
        v[0] = 1.0;
        return v;
    }

    // Choose sign to avoid cancellation: v = x + sign(x0)*||x||*e1
    double sign = (x[0] >= 0.0) ? 1.0 : -1.0;
    v[0] = x[0] + sign * x_norm;
    for (std::size_t i = 1; i < len; ++i)
        v[i] = x[i];

    // Normalize to unit length
    double v_norm = v.norm();
    for (std::size_t i = 0; i < len; ++i)
        v[i] /= v_norm;

    return v;
}

void apply_householder_left(Matrix& A, const Vector& v,
                            std::size_t row_start, std::size_t col_start, std::size_t col_end) {
    std::size_t len = v.size();

    // For each column j: compute w_j = v^T * A(:,j), then A(:,j) -= 2*v*w_j
    for (std::size_t col = col_start; col < col_end; ++col) {
        double projection = 0.0;
        for (std::size_t i = 0; i < len; ++i)
            projection += v[i] * A(row_start + i, col);
        for (std::size_t i = 0; i < len; ++i)
            A(row_start + i, col) -= 2.0 * v[i] * projection;
    }
}

void apply_householder_right(Matrix& A, const Vector& v,
                             std::size_t row_start, std::size_t row_end, std::size_t col_start) {
    std::size_t len = v.size();

    // For each row i: compute w_i = A(i,:) * v, then A(i,:) -= 2*w_i*v^T
    for (std::size_t row = row_start; row < row_end; ++row) {
        double projection = 0.0;
        for (std::size_t j = 0; j < len; ++j)
            projection += A(row, col_start + j) * v[j];
        for (std::size_t j = 0; j < len; ++j)
            A(row, col_start + j) -= 2.0 * projection * v[j];
    }
}
