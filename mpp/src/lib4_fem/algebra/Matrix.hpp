#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "AlgebraFwd.hpp"
#include "Operator.hpp"
#include "RVector.hpp"
#include "VectorMatrixBase.hpp"

class Matrix : public VectorMatrixBase, public Operator {
  BasicVector data;
  const Vector &Vec;
  SparseMatrix *sm;
  bool collectAfterProd = true;
  bool addBC = true;
public:
  mutable int applications = 0;

  void setCollectAfterProd(bool);

  double Max() const { return data.MaxAccumulated(CommSplit()); }

  double Min() const { return data.MinAccumulated(CommSplit()); }

  double norm() const { return data.normAccumulated(CommSplit()); }

  Matrix(const Vector &U, bool correctionAddBCArgyris = true);

  Matrix(const Matrix &A);

  Matrix(Matrix &&A);

  virtual ~Matrix();

  void PrintInfo() const { graph.MatrixMemoryInfo().PrintInfo(); }

  const SparseMatrix *GetSparse() const { return sm; }

  void CreateSparse();

  Matrix &operator=(const Matrix &A);

  Matrix &operator=(Matrix &&A);

  Matrix &operator=(Scalar b) {
    data = b;
    return *this;
  }

  const Vector &GetVector() const { return Vec; }

  int size() const { return VectorMatrixBase::size(); }

  const BasicVector &GetData() const { return data; }

  BasicVector &GetData() { return data; }

  const Scalar *operator()(int d) const { return data() + Entry(d); }

  Scalar *operator()(int d) { return data() + Entry(d); }

  const Scalar *operator()(const row &r0, const row &r1) const { return data() + GetEntry(r0, r1); }

  Scalar *operator()(const row &r0, const row &r1) { return data() + GetEntry(r0, r1); }

  Matrix &operator+=(const constAB<Scalar, Matrix> &aA);

  Matrix &operator+=(const Matrix &A) { return *this += constAB<Scalar, Matrix>(1.0, A); };

  Matrix &operator*=(const Scalar &a);

  void multiply(Vector &u, const Vector &v) const override;

  void multiply(Vectors &u, const Vectors &v) const override;

  void multiply_plus(Vector &, const Vector &) const override;

  void multiply_transpose_plus(Vector &, const Vector &) const override;

  void multiply_plus(Vectors &, const Vectors &) const override;

  void multiply_transpose_plus(Vectors &, const Vectors &) const override;

  void multiply_minus(Vector &b, const Vector &u) const override;

  void copy(Scalar *, int *, int *) const;

  void Transpose();

  void Symmetric();

  bool IsSymmetric();

  void print() const;

  int rowsize() const { return VectorMatrixBase::rowsize(); }

  void ClearDirichletValues();

  void EliminateDirichlet();

  void Accumulate();
private:
  void CommunicateMatrix(ExchangeBuffer &exBuffer);

  void AccumulateIdentify();

  void AccumulateParallel();
protected:
  void plusMatTVec(Vector &b, const Vector &u) const;

  void plusMatTVecs(Vectors &b, const Vectors &u) const;

  void plusMatVec(Vector &b, const Vector &u) const;

  void matVecs(Vectors &u, const Vectors &v) const;

  void plusMatVecs(Vectors &b, const Vectors &u) const;

  void matVec(Vector &u, const Vector &v) const;

  void minusMatVec(Vector &b, const Vector &u) const;
};

std::ostream &operator<<(std::ostream &s, const Matrix &u);

constAB<Operator, Vector> operator*(const Matrix &A, const Vector &v);

constAB<Operator, Vectors> operator*(const Matrix &A, const Vectors &v);

constAB<Scalar, Matrix> operator*(const Scalar &, const Matrix &);

#endif // MATRIX_HPP
