#include "include/qr.h"
#include <cmath>
#include <vector>

// Compute Givens rotation coefficients (c, s) that zero out `lower`:
// [ c  -s ] [ upper ] = [ r ]
// [ s   c ] [ lower ]   [ 0 ]
static void givens_rotation(double upper, double lower, double& c, double& s) {
    if (lower == 0.0) {
        c = 1.0;
        s = 0.0;
    } else if (std::abs(lower) > std::abs(upper)) {
        double tau = -upper / lower;
        s = 1.0 / std::sqrt(1.0 + tau * tau);
        c = s * tau;
    } else {
        double tau = -lower / upper;
        c = 1.0 / std::sqrt(1.0 + tau * tau);
        s = c * tau;
    }
}

QRResult qr_factorize_hessenberg(const Matrix& H) {
    std::size_t n = H.size();
    Matrix R = H;
    Matrix Q = Matrix::identity(n);

    for (std::size_t row = 0; row < n - 1; ++row) {
        double c, s;
        givens_rotation(R(row, row), R(row + 1, row), c, s);

        // Apply G^T from the left to rows (row, row+1)
        for (std::size_t col = row; col < n; ++col) {
            double val_curr = R(row, col);
            double val_next = R(row + 1, col);
            R(row, col)     =  c * val_curr - s * val_next;
            R(row + 1, col) =  s * val_curr + c * val_next;
        }

        // Accumulate Q = Q * G (columns row, row+1)
        for (std::size_t i = 0; i < n; ++i) {
            double val_curr = Q(i, row);
            double val_next = Q(i, row + 1);
            Q(i, row)     =  c * val_curr - s * val_next;
            Q(i, row + 1) =  s * val_curr + c * val_next;
        }
    }

    return {Q, R};
}

// Check if subdiagonal entry H(row, row-1) is negligible relative to its neighbors.
static bool is_negligible(const Matrix& H, int row, double tolerance) {
    double scale = std::abs(H(row - 1, row - 1)) + std::abs(H(row, row));
    if (scale == 0.0) scale = 1.0;
    return std::abs(H(row, row - 1)) < tolerance * scale;
}

// Check if a 2x2 block has complex eigenvalues (negative discriminant).
static bool has_complex_eigenvalues(double h11, double h12, double h21, double h22) {
    double discriminant = (h11 - h22) * (h11 - h22) + 4.0 * h12 * h21;
    return discriminant < 0.0;
}

// Compute Wilkinson shift: eigenvalue of the trailing 2x2 block closest to h22.
static double wilkinson_shift(double h11, double h12, double h21, double h22) {
    double trace = h11 + h22;
    double det = h11 * h22 - h12 * h21;
    double discriminant = trace * trace - 4.0 * det;

    if (discriminant < 0.0)
        return h22;

    double sqrt_disc = std::sqrt(discriminant);
    double eigenvalue_1 = (trace + sqrt_disc) / 2.0;
    double eigenvalue_2 = (trace - sqrt_disc) / 2.0;

    return (std::abs(eigenvalue_1 - h22) < std::abs(eigenvalue_2 - h22))
        ? eigenvalue_1
        : eigenvalue_2;
}

// Apply one shifted QR step via Givens rotations on H[block_start..block_end].
static void qr_step(Matrix& H, int block_start, int block_end, double shift) {
    for (int i = block_start; i <= block_end; ++i)
        H(i, i) -= shift;

    int block_size = block_end - block_start;
    std::vector<double> cosines(block_size), sines(block_size);

    // QR factorization via Givens: G_{n-1}^T ... G_1^T * (H - shift*I) = R
    for (int row = block_start; row < block_end; ++row) {
        double c, s;
        givens_rotation(H(row, row), H(row + 1, row), c, s);
        cosines[row - block_start] = c;
        sines[row - block_start] = s;

        for (int col = row; col <= block_end; ++col) {
            double val_curr = H(row, col);
            double val_next = H(row + 1, col);
            H(row, col)     =  c * val_curr - s * val_next;
            H(row + 1, col) =  s * val_curr + c * val_next;
        }
    }

    // Form R * Q by applying Givens rotations from the right
    for (int i = block_start; i < block_end; ++i) {
        double c = cosines[i - block_start];
        double s = sines[i - block_start];
        for (int row = block_start; row <= i + 1; ++row) {
            double val_curr = H(row, i);
            double val_next = H(row, i + 1);
            H(row, i)     =  c * val_curr - s * val_next;
            H(row, i + 1) =  s * val_curr + c * val_next;
        }
    }

    for (int i = block_start; i <= block_end; ++i)
        H(i, i) += shift;
}

Matrix qr_iterate(Matrix H, double tolerance, int max_iterations) {
    int n = static_cast<int>(H.size());
    if (n <= 1)
        return H;

    int active_end = n - 1;

    for (int iter = 0; iter < max_iterations && active_end > 0; ++iter) {
        // Deflate converged real eigenvalues from the bottom
        while (active_end > 0 && is_negligible(H, active_end, tolerance)) {
            H(active_end, active_end - 1) = 0.0;
            --active_end;
        }
        if (active_end == 0) break;

        // Deflate converged 2x2 complex block at the bottom
        if (active_end >= 2 && is_negligible(H, active_end - 1, tolerance)) {
            double h11 = H(active_end - 1, active_end - 1);
            double h12 = H(active_end - 1, active_end);
            double h21 = H(active_end, active_end - 1);
            double h22 = H(active_end, active_end);
            if (has_complex_eigenvalues(h11, h12, h21, h22)) {
                H(active_end - 1, active_end - 2) = 0.0;
                active_end -= 2;
                continue;
            }
        }

        // Irreducible 2x2 complex block — nothing left to do
        if (active_end == 1) {
            double h11 = H(0, 0), h12 = H(0, 1);
            double h21 = H(1, 0), h22 = H(1, 1);
            if (has_complex_eigenvalues(h11, h12, h21, h22)) break;
        }

        // Find the top of the active unreduced block
        int active_start = active_end - 1;
        while (active_start > 0 && !is_negligible(H, active_start, tolerance))
            --active_start;
        if (active_start > 0) H(active_start, active_start - 1) = 0.0;

        double shift = wilkinson_shift(
            H(active_end - 1, active_end - 1), H(active_end - 1, active_end),
            H(active_end, active_end - 1),     H(active_end, active_end));

        qr_step(H, active_start, active_end, shift);
    }

    return H;
}
