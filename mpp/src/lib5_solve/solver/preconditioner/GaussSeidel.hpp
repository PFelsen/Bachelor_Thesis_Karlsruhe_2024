#ifndef GAUSSSEIDEL_HPP
#define GAUSSSEIDEL_HPP

#include "Jacobi.hpp"

class GaussSeidel : public Preconditioner {
  Matrix *A;
  SparseMatrix *Sp;
  bool shift_special;
public:
  GaussSeidel();

  void Destruct() override;

  void Construct(const Matrix &_A) override;

  virtual ~GaussSeidel() { Destruct(); }

  void multiply(Vector &u, const Vector &b) const override;

  void multiply_transpose(Vector &u, const Vector &b) const override;

  virtual string Name() const { return "GaussSeidel"; }

  friend std::ostream &operator<<(std::ostream &s, const GaussSeidel &GS);
};

class BlockGaussSeidel : public BlockJacobi {
protected:
  const Matrix *A;
private:
  int *transposed_entries;
  int *transposed_j;
  int *transposed_nj;
public:
  BlockGaussSeidel() : A(0), transposed_entries(0), transposed_j(0), transposed_nj(0) {}

  void Construct(const Matrix &_A) override;

  void Destruct() override;

  virtual ~BlockGaussSeidel() { Destruct(); }

  virtual void multiply(Vector &u, const Vector &b) const override;

  virtual string Name() const { return "BlockGaussSeidel"; }
};

class PointBlockGaussSeidel : public PointBlockJacobi {
  mutable bool forward = true;
protected:
  const Matrix *A = 0;
public:
  PointBlockGaussSeidel(bool forward = true);

  void Construct(const Matrix &_A) override;

  void Destruct() override;

  ~PointBlockGaussSeidel() { Destruct(); }

  void multiply(Vector &u, const Vector &b) const override;

  void multiply_transpose(Vector &u, const Vector &b) const override;

  string Name() const { return "PointBlockGaussSeidel"; }

  void transpose() const override { forward = !forward; }
};

#endif // GAUSSSEIDEL_HPP
