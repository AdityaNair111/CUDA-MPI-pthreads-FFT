//
// Created by brian on 11/20/18.
//

#pragma once

#include <iostream>

class Complex
{
public:
  Complex();
  Complex(float r, float i);
  Complex(float r);
  Complex operator+(const Complex &b) const;
  Complex operator-(const Complex &b) const;
  Complex operator*(const Complex &b) const;

  float mag() const;
  float angle() const;
  Complex conj() const;

  float real;
  float imag;
};

std::ostream &operator<<(std::ostream &os, const Complex &rhs);
