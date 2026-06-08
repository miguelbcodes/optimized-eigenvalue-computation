#include "include/householder.h"
#include <cmath>

Vector householder_vector(const Vector& x) {
    std::size_t m = x.size();
    Vector v(m);

    double x_norm = x.norm();
    if (x_norm == 0.0) {
        v[0] = 1.0;
        return v;
    }

    double sign = (x[0] >= 0.0) ? 1.0 : -1.0;
    v[0] = x[0] + sign * x_norm;
    for (std::size_t i = 1; i < m; ++i)
        v[i] = x[i];

    double v_norm = v.norm();
    for (std::size_t i = 0; i < m; ++i)
        v[i] /= v_norm;

    return v;
}

void apply_householder_left(Matrix& A, const Vector& v,
                            std::size_t row_start, std::size_t col_start, std::size_t col_end) {
    // A(row_start:row_start+m, col_start:col_end) -= 2 * v * (v^T * A_sub)
    std::size_t m = v.size();

    for (std::size_t j = col_start; j < col_end; ++j) {
        double dot = 0.0;
        for (std::size_t i = 0; i < m; ++i)
            dot += v[i] * A(row_start + i, j);
        for (std::size_t i = 0; i < m; ++i)
            A(row_start + i, j) -= 2.0 * v[i] * dot;
    }
}

void apply_householder_right(Matrix& A, const Vector& v,
                             std::size_t row_start, std::size_t row_end, std::size_t col_start) {
    // A(row_start:row_end, col_start:col_start+m) -= 2 * (A_sub * v) * v^T
    std::size_t m = v.size();

    for (std::size_t i = row_start; i < row_end; ++i) {
        double dot = 0.0;
        for (std::size_t j = 0; j < m; ++j)
            dot += A(i, col_start + j) * v[j];
        for (std::size_t j = 0; j < m; ++j)
            A(i, col_start + j) -= 2.0 * dot * v[j];
    }
}
