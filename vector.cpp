#include "include/vector.h"
#include <cmath>
#include <iomanip>

Vector::Vector() : n_(0) {}

Vector::Vector(std::size_t n) : n_(n), data_(n, 0.0) {}

std::size_t Vector::size() const { return n_; }

double& Vector::operator[](std::size_t i) { return data_[i]; }

double Vector::operator[](std::size_t i) const { return data_[i]; }

Vector Vector::operator*(double scalar) const {
    Vector result(n_);
    for (std::size_t i = 0; i < n_; ++i)
        result.data_[i] = data_[i] * scalar;
    return result;
}

Vector Vector::operator+(const Vector& other) const {
    Vector result(n_);
    for (std::size_t i = 0; i < n_; ++i)
        result.data_[i] = data_[i] + other.data_[i];
    return result;
}

Vector Vector::operator-(const Vector& other) const {
    Vector result(n_);
    for (std::size_t i = 0; i < n_; ++i)
        result.data_[i] = data_[i] - other.data_[i];
    return result;
}

double Vector::dot(const Vector& other) const {
    double sum = 0.0;
    for (std::size_t i = 0; i < n_; ++i)
        sum += data_[i] * other.data_[i];
    return sum;
}

double Vector::norm() const {
    return std::sqrt(dot(*this));
}

void Vector::print(std::ostream& out) const {
    out << std::fixed << std::setprecision(4);
    for (std::size_t i = 0; i < n_; ++i)
        out << std::setw(10) << data_[i] << '\n';
}

Vector operator*(double scalar, const Vector& v) {
    return v * scalar;
}
