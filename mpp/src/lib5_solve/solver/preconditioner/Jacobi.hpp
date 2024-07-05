#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "PreconditionerBase.hpp"
#include "SparseRMatrix.hpp"

class Jacobi : public Preconditioner {
protected:
  Vector *D;
  double theta;
  bool shift_special;
public:
  Jacobi();

  void Construct(const Matrix &A) override;

  void Destruct() override;

  virtual ~Jacobi() { Destruct(); }

  virtual void multiply(Vector &u, const Vector &b) const override;

  virtual string Name() const { return "Jacobi"; }

  friend std::ostream &operator<<(std::ostream &s, const Jacobi &Jac) { return s << *(Jac.D); }
};

class DampedJacobi : public Jacobi {
  double damp;
public:
  DampedJacobi() {
    damp = 1.0;
    Config::Get("PreconditionerDamp", damp);
  }

  virtual ~DampedJacobi() { Destruct(); }

  void multiply(Vector &u, const Vector &b) const override;

  string Name() const { return "DampedJacobi"; }
};

class BlockJacobi : public Preconditioner {
protected:
  Scalar *dd;
  std::unordered_map<int, int> dd_entry;
  int blocksize;
public:
  BlockJacobi();

  void Construct(const Matrix &A) override;

  void Construct_Matrix_free(const Matrix &A);

  void Destruct() override;

  virtual ~BlockJacobi() { Destruct(); }

  virtual void multiply(Vector &u, const Vector &b) const override;

  virtual void multiply_transpose(Vector &u, const Vector &b) const override;

  virtual string Name() const { return "BlockJacobi"; }
};

class PointBlockJacobi : public Preconditioner {
protected:
  Scalar *dd;
  int *dd_entry;
  double theta;
  bool LAPACKFlag = false;
public:
  PointBlockJacobi();

  void apply(int n, Scalar *u, const Scalar *a, const Scalar *b) const;

  void Construct(const Matrix &A) override;

  void Destruct() override;

  virtual ~PointBlockJacobi() { Destruct(); }

  virtual void multiply(Vector &u, const Vector &b) const override;

  virtual void multiply_transpose(Vector &u, const Vector &b) const override;

  virtual string Name() const { return "PointBlockJacobi"; }
};

class PointBlockJacobi_dG : public Preconditioner {
protected:
  Matrix *D;
  double theta;
public:
  PointBlockJacobi_dG();

  void Construct(const Matrix &A) override;

  void Destruct() override;

  virtual ~PointBlockJacobi_dG() { Destruct(); }

  virtual void multiply(Vector &u, const Vector &b) const override;

  virtual string Name() const { return "PointBlockJacobi_dG"; }

  friend std::ostream &operator<<(std::ostream &s, const PointBlockJacobi_dG &Jac) {
    return s << *(Jac.D);
  }
};

class CellWiseJacobi : public Preconditioner {
protected:
  std::unordered_map<Point, SparseRMatrix> InvMap;
  double theta;
public:
  CellWiseJacobi();

  void Construct(const Matrix &A) override;

  void Destruct() override;

  virtual ~CellWiseJacobi() { Destruct(); }

  virtual void multiply(Vector &u, const Vector &b) const override;

  virtual string Name() const { return "CellWiseJacobi"; }
};
#endif // JACOBI_HPP
