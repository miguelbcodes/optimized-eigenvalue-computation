#include "include/hessenberg.h"
#include "include/householder.h"

Matrix reduce_to_hessenberg(Matrix A) {
    std::size_t n = A.size();
    if (n <= 2)
        return A;

    for (std::size_t k = 0; k < n - 2; ++k) {
        std::size_t m = n - k - 1;
        Vector x(m);
        for (std::size_t i = 0; i < m; ++i)
            x[i] = A(k + 1 + i, k);

        Vector v = householder_vector(x);

        apply_householder_left(A, v, k + 1, k, n);
        apply_householder_right(A, v, 0, n, k + 1);
    }

    return A;
}
