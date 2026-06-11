#ifndef VECTOR_H
#define VECTOR_H

#include <cstddef>
#include <iostream>
#include <vector>

// Dense vector of doubles with basic linear algebra operations.
class Vector {
public:
    Vector();
    explicit Vector(std::size_t length);

    std::size_t size() const;

    double& operator[](std::size_t index);
    double operator[](std::size_t index) const;

    Vector operator*(double scalar) const;
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;

    double dot(const Vector& other) const;
    // Euclidean (L2) norm.
    double norm() const;

    void print(std::ostream& out = std::cout) const;

private:
    std::size_t n_;
    std::vector<double> data_;
};

Vector operator*(double scalar, const Vector& v);

#endif // VECTOR_H
