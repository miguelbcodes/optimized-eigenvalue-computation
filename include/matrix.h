#ifndef MATRIX_H
#define MATRIX_H

#include <cstddef>
#include <iostream>
#include <vector>

class Matrix {
public:
    Matrix();
    explicit Matrix(std::size_t n);

    std::size_t size() const;

    double& operator()(std::size_t i, std::size_t j);
    double operator()(std::size_t i, std::size_t j) const;

    Matrix operator*(const Matrix& other) const;
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;

    Matrix transpose() const;

    static Matrix identity(std::size_t n);

    void print(std::ostream& out = std::cout) const;
    static void fill_random(Matrix& m);

private:
    std::size_t n_;
    std::vector<double> data_;
};

#endif // MATRIX_H
