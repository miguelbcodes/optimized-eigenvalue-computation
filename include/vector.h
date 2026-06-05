#ifndef VECTOR_H
#define VECTOR_H

#include <cstddef>
#include <iostream>
#include <vector>

class Vector {
public:
    Vector();
    explicit Vector(std::size_t n);

    std::size_t size() const;

    double& operator[](std::size_t i);
    double operator[](std::size_t i) const;

    Vector operator*(double scalar) const;
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;

    double dot(const Vector& other) const;
    double norm() const;

    void print(std::ostream& out = std::cout) const;

private:
    std::size_t n_;
    std::vector<double> data_;
};

Vector operator*(double scalar, const Vector& v);

#endif // VECTOR_H
