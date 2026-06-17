#include "../include/tridiagonal.h"
#include "../include/householder.h"
#include <cmath>

SymmetricTridiagonal reduce_to_tridiagonal(Matrix A) {
    std::size_t n = A.size();
    SymmetricTridiagonal T;
    T.diag.assign(n, 0.0);
    T.offdiag.assign(n > 0 ? n - 1 : 0, 0.0);

    // n - 2 Householder steps drive the matrix to tridiagonal form. Each step
    // zeroes the entries below the subdiagonal in one column.
    for (std::size_t col = 0; col + 2 < n; ++col) {
        std::size_t m = n - col - 1;  // length of subcolumn and trailing block

        // Subcolumn to be reduced: A(col+1 : n-1, col).
        Vector x(m);
        for (std::size_t i = 0; i < m; ++i)
            x[i] = A(col + 1 + i, col);

        Vector v = householder_vector(x);  // unit reflector, P = I - 2*v*v^T

        // Symmetric similarity update of the trailing block
        // B = A(col+1:, col+1:): with p = 2*B*v and w = p - (v^T*p)*v,
        //   P*B*P = B - w*v^T - v*w^T.
        Vector p(m);
        for (std::size_t i = 0; i < m; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < m; ++j)
                sum += A(col + 1 + i, col + 1 + j) * v[j];
            p[i] = 2.0 * sum;
        }

        double vtp = 0.0;
        for (std::size_t i = 0; i < m; ++i)
            vtp += v[i] * p[i];

        Vector w(m);
        for (std::size_t i = 0; i < m; ++i)
            w[i] = p[i] - vtp * v[i];

        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j)
                A(col + 1 + i, col + 1 + j) -= w[i] * v[j] + v[i] * w[j];

        // Apply the reflector to the subcolumn: P*x leaves only the first entry
        // nonzero. Store it as the new sub-/super-diagonal coupling and zero the
        // rest to keep the structure explicit.
        double vtx = 0.0;
        for (std::size_t i = 0; i < m; ++i)
            vtx += v[i] * x[i];
        double new_offdiag = x[0] - 2.0 * v[0] * vtx;

        A(col + 1, col) = new_offdiag;
        A(col, col + 1) = new_offdiag;
        for (std::size_t i = 1; i < m; ++i) {
            A(col + 1 + i, col) = 0.0;
            A(col, col + 1 + i) = 0.0;
        }
    }

    for (std::size_t i = 0; i < n; ++i)
        T.diag[i] = A(i, i);
    for (std::size_t i = 0; i + 1 < n; ++i)
        T.offdiag[i] = A(i + 1, i);

    return T;
}

// Implicit-shift QL with Wilkinson shift and deflation (cf. tqli, Numerical
// Recipes). Operates in place on the diagonal d and off-diagonal e bands; e[i]
// couples d[i] and d[i+1], and e[n-1] is unused workspace. On return d holds the
// eigenvalues (unordered).
static void tql_implicit_shift(std::vector<double>& d, std::vector<double>& e,
                               int n, double tolerance, int max_iterations) {
    for (int l = 0; l < n; ++l) {
        int iter = 0;
        int m;
        do {
            // Look for a negligible off-diagonal entry that splits the matrix.
            for (m = l; m < n - 1; ++m) {
                double scale = std::abs(d[m]) + std::abs(d[m + 1]);
                if (std::abs(e[m]) <= tolerance * scale)
                    break;
            }
            if (m != l) {
                if (iter++ == max_iterations)
                    return;  // give up gracefully rather than spin forever

                // Wilkinson shift from the leading 2x2 of the active block.
                double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                double r = std::hypot(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + std::copysign(r, g));

                double s = 1.0, c = 1.0, p = 0.0;
                int i = m - 1;
                for (; i >= l; --i) {
                    double f = s * e[i];
                    double b = c * e[i];
                    r = std::hypot(f, g);
                    e[i + 1] = r;
                    if (r == 0.0) {  // recover from underflow
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;
                }
                if (r == 0.0 && i >= l)
                    continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }
}

std::vector<double> qr_iterate_tridiagonal(SymmetricTridiagonal T,
                                           double tolerance,
                                           int max_iterations) {
    int n = static_cast<int>(T.diag.size());
    std::vector<double> d = T.diag;
    if (n <= 1)
        return d;

    // Pad the off-diagonal to length n; the last slot is scratch used by the sweep.
    std::vector<double> e(n, 0.0);
    for (int i = 0; i + 1 < n; ++i)
        e[i] = T.offdiag[i];

    tql_implicit_shift(d, e, n, tolerance, max_iterations);
    return d;
}
