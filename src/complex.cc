//
// Created by brian on 11/20/18.
//

#include "complex.h"

#include <cmath>

const float PI = 3.14159265358979f;

Complex::Complex() : real(0.0f), imag(0.0f) {}

Complex::Complex(float r) : real(r), imag(0.0f) {}

Complex::Complex(float r, float i) : real(r), imag(i) {}

Complex Complex::operator+(const Complex &b) const
{
    Complex a;
    a.real = real + b.real;
    a.imag = imag + b.imag;
    return a;
}

Complex Complex::operator-(const Complex &b) const
{
    Complex a;
    a.real = real - b.real;
    a.imag = imag - b.imag;
    return a;
}

Complex Complex::operator*(const Complex &b) const
{
    Complex a;
    a.real = real * b.real - imag * b.imag;
    a.imag = imag * b.real + real * b.imag;
    return a;
}

float Complex::mag() const
{
    float a;
    a = pow(real * real + imag * imag, 0.5);
    return a;
}

float Complex::angle() const
{
    float a;
    a = atan(imag / real);
    return a;
}

Complex Complex::conj() const
{
    Complex a;
    a.real = real;
    a.imag = -imag;
    return a;
}

std::ostream &operator<<(std::ostream &os, const Complex &rhs)
{
    Complex c(rhs);
    if (fabsf(rhs.imag) < 1e-10)
        c.imag = 0.0f;
    if (fabsf(rhs.real) < 1e-10)
        c.real = 0.0f;

    if (c.imag == 0)
        os << c.real;
    else
        os << "(" << c.real << "," << c.imag << ")";
    return os;
}