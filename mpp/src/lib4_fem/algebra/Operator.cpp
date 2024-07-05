#include "Operator.hpp"
#include "Algebra.hpp"

void Operator::apply(Vector &u, Vector &v) {
  u = 0;
  apply_plus(u, v);
}

void Operator::apply_minus(Vector &u, Vector &v) {
  Vector tmp(u);
  apply(tmp, v);
  u -= tmp;
}

void Operator::apply_transpose(Vector &u, Vector &v) {
  u = 0;
  apply_transpose_plus(u, v);
}

void Operator::apply_transpose_minus(Vector &u, Vector &v) {
  Vector tmp(u);
  apply_transpose(tmp, v);
  u -= tmp;
}

void Operator::multiply(Vector &u, const Vector &v) const {
  u = 0;
  multiply_plus(u, v);
}

void Operator::multiply_minus(Vector &u, const Vector &v) const {
  Vector tmp(u);
  multiply(tmp, v);
  u -= tmp;
}

void Operator::multiply_transpose(Vector &u, const Vector &v) const {
  u = 0;
  multiply_transpose_plus(u, v);
}

void Operator::multiply_transpose_minus(Vector &u, const Vector &v) const {
  Vector tmp(u);
  multiply_transpose(tmp, v);
  u -= tmp;
}

//**** --> with Vectors
void Operator::apply(Vectors &u, Vectors &v) {
  for (int n = 0; n < u.size(); ++n)
    apply(u[n], v[n]);
}

void Operator::apply_minus(Vectors &u, Vectors &v) {
  for (int n = 0; n < u.size(); ++n)
    apply_minus(u[n], v[n]);
}

void Operator::apply_transpose(Vectors &u, Vectors &v) {
  for (int n = 0; n < u.size(); ++n)
    apply_transpose(u[n], v[n]);
}

void Operator::apply_transpose_minus(Vectors &u, Vectors &v) {
  for (int n = 0; n < u.size(); ++n)
    apply_transpose_minus(u[n], v[n]);
}

void Operator::multiply(Vectors &u, const Vectors &v) const {
  //     u = 0; multiply_plus(u,v);
  //     return;
  for (int n = 0; n < u.size(); ++n)
    multiply(u[n], v[n]);
}

void Operator::multiply_minus(Vectors &u, const Vectors &v) const {
  for (int n = 0; n < u.size(); ++n)
    multiply_minus(u[n], v[n]);
}

void Operator::multiply_transpose(Vectors &u, const Vectors &v) const {
  for (int n = 0; n < u.size(); ++n)
    multiply_transpose(u[n], v[n]);
}

void Operator::multiply_transpose_minus(Vectors &u, const Vectors &v) const {
  for (int n = 0; n < u.size(); ++n)
    multiply_transpose_minus(u[n], v[n]);
}

void Operator::apply_plus(Vectors &u, Vectors &v) {
  for (int n = 0; n < u.size(); ++n)
    apply_plus(u[n], v[n]);
}

void Operator::apply_transpose_plus(Vectors &u, Vectors &v) {
  for (int n = 0; n < u.size(); ++n)
    apply_transpose_plus(u[n], v[n]);
}

void Operator::multiply_plus(Vectors &u, const Vectors &v) const {
  for (int n = 0; n < u.size(); ++n)
    multiply_plus(u[n], v[n]);
}

void Operator::multiply_transpose_plus(Vectors &u, const Vectors &v) const {
  for (int n = 0; n < u.size(); ++n)
    multiply_transpose_plus(u[n], v[n]);
}
