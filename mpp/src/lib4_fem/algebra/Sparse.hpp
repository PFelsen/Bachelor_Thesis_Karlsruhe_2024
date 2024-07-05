#ifndef _SPARSE_H_
#define _SPARSE_H_

#include "AlgebraFwd.hpp"

#include "BasicSparseSolver.hpp"
#include "Operator.hpp"

// Needed?
#include <cmath>

template<class D>
bool isZero(D z) {
  return abs(z) < nearToZero;
}

template<class D>
bool isnotZero(D z) {
  return abs(z) >= nearToZero;
}

class SparseMatrix : public SparseRMatrix, public Operator {
  int size_(const SparseMatrix &M, const vector<bool> &mask) const;

  int Size_(const SparseMatrix &M, const vector<bool> &mask) const;

  int Size_(const SparseMatrix &M, const vector<int> &rindex) const;

  int Size_(const SparseMatrix &M, const vector<int> &indexrow, const vector<int> &indexcol) const;

  int int_in_vector(int x, const vector<int> &y) const;

  int _Size(const SparseMatrix &M, const vector<int> &Indices) const;

  typedef std::map<int, Scalar> SparseVecT;
  typedef SparseVecT::iterator SparVecIt;
  typedef SparseVecT::const_iterator conSparVecIt;
  typedef vector<SparseVecT> SparseMatT;
  typedef vector<int> DofDofId;
  DofDofId ParDofId, ParDofRes;

  void convert_sm(SparseMatT &) const;

  void print_sm(const SparseMatT &) const;

  void print_sm_dense(const SparseMatT &) const;

  void build_ident(const Matrix &);

  void apply_PC_mat(SparseMatT &, const SparseMatT &) const;
public:
  void convert_sm_back(const SparseMatT &);

  SparseMatrix(const Matrix &M);

  SparseMatrix(const SparseMatrix &);

  SparseMatrix(const SparseMatrix &M, const vector<bool> &mask);

  SparseMatrix(const SparseMatrix &M, const vector<int> &Indices);

  SparseMatrix(const SparseMatrix &M, const vector<int> &index, const vector<int> &rindex);

  SparseMatrix(const SparseMatrix &M, const vector<int> &indexrow, const vector<int> &indexcol,
               bool dummy);

  void multiply_plus(Vector &b, const Vector &u) const;

  void pc_mat_convert(const Matrix &);

  void ShrinkIdentify(Vector &, const Vector &) const;

  void ExpandIdentify(Vector &) const;

  Scalar *ref(int j) const;

  void plusMatVec(Vector &b, const Vector &u) const;

  void minusMatVec(Vector &b, const Vector &u) const;

  void multiply_minus(Vector &b, const Vector &u) const;
};

constAB<Operator, Vector> operator*(const SparseMatrix &S, const Vector &v);

constAB<Operator, Vectors> operator*(const SparseMatrix &S, const Vectors &v);

#endif // of #ifndef _SPARSE_H_
