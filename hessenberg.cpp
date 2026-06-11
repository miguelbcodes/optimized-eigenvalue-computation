#include "include/hessenberg.h"
#include "include/householder.h"

Matrix reduce_to_hessenberg(Matrix A) {
    std::size_t n = A.size();
    if (n <= 2)
        return A;

    for (std::size_t col = 0; col < n - 2; ++col) {
        // Extract the subcolumn A(col+1 : n-1, col) to be zeroed out
        std::size_t subvec_len = n - col - 1;
        Vector subcolumn(subvec_len);
        for (std::size_t i = 0; i < subvec_len; ++i)
            subcolumn[i] = A(col + 1 + i, col);

        Vector reflector = householder_vector(subcolumn);

        // Similarity transformation: A <- Q_k^T * A * Q_k
        apply_householder_left(A, reflector, col + 1, col, n);
        apply_householder_right(A, reflector, 0, n, col + 1);
    }

    return A;
}
