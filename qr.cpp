#include "include/qr.h"
#include <cmath>
#include <algorithm>

static void givens_rotation(double a, double b, double& c, double& s) {
    if (b == 0.0) {
        c = 1.0;
        s = 0.0;
    } else if (std::abs(b) > std::abs(a)) {
        double t = -a / b;
        s = 1.0 / std::sqrt(1.0 + t * t);
        c = s * t;
    } else {
        double t = -b / a;
        c = 1.0 / std::sqrt(1.0 + t * t);
        s = c * t;
    }
}

QRResult qr_factorize_hessenberg(const Matrix& H) {
    std::size_t n = H.size();
    Matrix R = H;
    Matrix Q = Matrix::identity(n);

    for (std::size_t i = 0; i < n - 1; ++i) {
        double c, s;
        givens_rotation(R(i, i), R(i + 1, i), c, s);

        for (std::size_t j = i; j < n; ++j) {
            double r0 = R(i, j);
            double r1 = R(i + 1, j);
            R(i, j)     =  c * r0 - s * r1;
            R(i + 1, j) =  s * r0 + c * r1;
        }

        for (std::size_t j = 0; j < n; ++j) {
            double q0 = Q(j, i);
            double q1 = Q(j, i + 1);
            Q(j, i)     =  c * q0 - s * q1;
            Q(j, i + 1) =  s * q0 + c * q1;
        }
    }

    return {Q, R};
}

static bool is_negligible(const Matrix& H, int i, double epsilon) {
    double scale = std::abs(H(i - 1, i - 1)) + std::abs(H(i, i));
    if (scale == 0.0) scale = 1.0;
    return std::abs(H(i, i - 1)) < epsilon * scale;
}

// Apply one implicit shifted QR step on the active submatrix H[l..p, l..p]
static void qr_step(Matrix& H, int l, int p, double shift) {
    for (int i = l; i <= p; ++i)
        H(i, i) -= shift;

    std::vector<double> cs(p - l), ss(p - l);

    for (int i = l; i < p; ++i) {
        double c, s;
        givens_rotation(H(i, i), H(i + 1, i), c, s);
        cs[i - l] = c;
        ss[i - l] = s;

        for (int j = i; j <= p; ++j) {
            double r0 = H(i, j);
            double r1 = H(i + 1, j);
            H(i, j)     =  c * r0 - s * r1;
            H(i + 1, j) =  s * r0 + c * r1;
        }
    }

    for (int i = l; i < p; ++i) {
        double c = cs[i - l], s = ss[i - l];
        for (int j = l; j <= i + 1; ++j) {
            double r0 = H(j, i);
            double r1 = H(j, i + 1);
            H(j, i)     =  c * r0 - s * r1;
            H(j, i + 1) =  s * r0 + c * r1;
        }
    }

    for (int i = l; i <= p; ++i)
        H(i, i) += shift;
}

Matrix qr_iterate(Matrix H, double epsilon, int max_iter) {
    int n = static_cast<int>(H.size());
    if (n <= 1)
        return H;

    int p = n - 1;

    for (int iter = 0; iter < max_iter && p > 0; ++iter) {
        // Deflate converged entries from the bottom
        while (p > 0 && is_negligible(H, p, epsilon)) {
            H(p, p - 1) = 0.0;
            --p;
        }
        if (p == 0) break;

        // Check for converged 2x2 complex block at bottom
        if (p >= 2 && is_negligible(H, p - 1, epsilon)) {
            double a = H(p - 1, p - 1), b = H(p - 1, p);
            double c = H(p, p - 1),     d = H(p, p);
            double disc = (a - d) * (a - d) + 4.0 * b * c;
            if (disc < 0.0) {
                H(p - 1, p - 2) = 0.0;
                p -= 2;
                continue;
            }
        }
        if (p == 1) {
            double a = H(0, 0), b = H(0, 1);
            double c = H(1, 0), d = H(1, 1);
            double disc = (a - d) * (a - d) + 4.0 * b * c;
            if (disc < 0.0) break;
        }

        // Find the start of the active unreduced block
        int l = p - 1;
        while (l > 0 && !is_negligible(H, l, epsilon)) --l;
        if (l > 0) H(l, l - 1) = 0.0;

        // Wilkinson shift from trailing 2x2 of active block
        double a = H(p - 1, p - 1);
        double b = H(p - 1, p);
        double c = H(p, p - 1);
        double d = H(p, p);
        double trace = a + d;
        double det = a * d - b * c;
        double disc = trace * trace - 4.0 * det;
        double mu;
        if (disc < 0.0) {
            mu = d;
        } else {
            double sqrt_disc = std::sqrt(disc);
            double lam1 = (trace + sqrt_disc) / 2.0;
            double lam2 = (trace - sqrt_disc) / 2.0;
            mu = (std::abs(lam1 - d) < std::abs(lam2 - d)) ? lam1 : lam2;
        }

        qr_step(H, l, p, mu);
    }

    return H;
}
