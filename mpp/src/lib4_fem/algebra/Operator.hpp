#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include "AlgebraFwd.hpp"
#include "Assertion.hpp"

class Operator {
public:
  virtual void apply_plus(Vector &, Vector &) { THROW("apply_plus not implemented") }

  virtual void apply(Vector &, Vector &);

  virtual void apply_minus(Vector &, Vector &);

  virtual void apply_transpose_plus(Vector &, Vector &) {
    THROW("apply_transpose_plus not implemented")
  }

  virtual void apply_transpose(Vector &u, Vector &v);

  virtual void apply_transpose_minus(Vector &u, Vector &v);

  virtual void multiply_plus(Vector &, const Vector &) const {
    THROW("multiply_plus not implemented")
  }

  virtual void multiply(Vector &u, const Vector &v) const;

  virtual void multiply_minus(Vector &u, const Vector &v) const;

  virtual void multiply_transpose_plus(Vector &, const Vector &) const {
    THROW("multiply_transpose_plus not implemented")
  }

  virtual void multiply_transpose(Vector &u, const Vector &v) const;

  virtual void multiply_transpose_minus(Vector &u, const Vector &v) const;

  virtual ~Operator() = default;

  virtual void apply_plus(Vectors &u, Vectors &v);

  virtual void apply(Vectors &, Vectors &);

  virtual void apply_minus(Vectors &, Vectors &);

  virtual void apply_transpose_plus(Vectors &u, Vectors &v);

  virtual void apply_transpose(Vectors &u, Vectors &v);

  virtual void apply_transpose_minus(Vectors &u, Vectors &v);

  virtual void multiply_plus(Vectors &u, const Vectors &v) const;

  virtual void multiply(Vectors &u, const Vectors &v) const;

  virtual void multiply_minus(Vectors &u, const Vectors &v) const;

  virtual void multiply_transpose_plus(Vectors &u, const Vectors &v) const;

  virtual void multiply_transpose(Vectors &u, const Vectors &v) const;

  virtual void multiply_transpose_minus(Vectors &u, const Vectors &v) const;
};

#endif // OPERATOR_HPP
