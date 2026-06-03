#include "include/matrix.h"
#include <iomanip>

Matrix::Matrix() : n_(0) {}

Matrix::Matrix(std::size_t n) : n_(n), data_(n * n, 0.0) {}

std::size_t Matrix::size() const { return n_; }

double& Matrix::operator()(std::size_t i, std::size_t j) {
    return data_[i * n_ + j];
}

double Matrix::operator()(std::size_t i, std::size_t j) const {
    return data_[i * n_ + j];
}

Matrix Matrix::operator*(const Matrix& other) const {
    Matrix result(n_);
    for (std::size_t i = 0; i < n_; ++i)
        for (std::size_t k = 0; k < n_; ++k)
            for (std::size_t j = 0; j < n_; ++j)
                result.data_[i * n_ + j] += data_[i * n_ + k] * other.data_[k * n_ + j];
    return result;
}

Matrix Matrix::operator+(const Matrix& other) const {
    Matrix result(n_);
    for (std::size_t i = 0; i < n_ * n_; ++i)
        result.data_[i] = data_[i] + other.data_[i];
    return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
    Matrix result(n_);
    for (std::size_t i = 0; i < n_ * n_; ++i)
        result.data_[i] = data_[i] - other.data_[i];
    return result;
}

Matrix Matrix::transpose() const {
    Matrix result(n_);
    for (std::size_t i = 0; i < n_; ++i)
        for (std::size_t j = 0; j < n_; ++j)
            result.data_[j * n_ + i] = data_[i * n_ + j];
    return result;
}

Matrix Matrix::identity(std::size_t n) {
    Matrix m(n);
    for (std::size_t i = 0; i < n; ++i)
        m.data_[i * n + i] = 1.0;
    return m;
}

void Matrix::print(std::ostream& out) const {
    out << std::fixed << std::setprecision(4);
    for (std::size_t i = 0; i < n_; ++i) {
        for (std::size_t j = 0; j < n_; ++j) {
            if (j > 0) out << "  ";
            out << std::setw(10) << data_[i * n_ + j];
        }
        out << '\n';
    }
}
