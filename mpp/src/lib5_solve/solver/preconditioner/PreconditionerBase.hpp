#ifndef PRECONDITIONERBASE_HPP
#define PRECONDITIONERBASE_HPP

#include "Matrix.hpp"

class IAssemble;

class LinearSolver;

class Transfer;

class Preconditioner : public Operator {
protected:
  int verbose;
public:
  virtual void Construct(const Matrix &) = 0;

  virtual void Construct(const Matrix &A, const IAssemble &) { Construct(A); }

  virtual void Destruct() = 0;

  Preconditioner() : verbose(0) { Config::Get("PreconditionerVerbose", verbose); }

  virtual void multiply_transpose(Vector &u, const Vector &b) const {
    multiply(u, b);
    return;
  }

  virtual ~Preconditioner() {}

  // transpose is used for the adjoint equation. For example by a multigrid preconditioner.
  virtual void transpose() const {}

  // reset for preconditioners with state
  virtual void reset() {}

  virtual string Name() const = 0;

  friend std::ostream &operator<<(std::ostream &s, const Preconditioner &PC) {
    return s << PC.Name();
  }

  void SetVerbose(int v) { verbose = v; }
};

inline constAB<Operator, Vector> operator*(const Preconditioner &PC, const Vector &v) {
  return constAB<Operator, Vector>(PC, v);
}

inline constAB<Vector, Operator> operator*(const Vector &v, const Preconditioner &PC) {
  return constAB<Vector, Operator>(v, PC);
}

inline constAB<Operator, Vectors> operator*(const Preconditioner &PC, const Vectors &v) {
  return constAB<Operator, Vectors>(PC, v);
}

inline constAB<Vectors, Operator> operator*(const Vectors &v, const Preconditioner &PC) {
  return constAB<Vectors, Operator>(v, PC);
}

#endif // PRECONDITIONERBASE_HPP
