#ifndef BASICOPERATOR_HPP
#define BASICOPERATOR_HPP

#include "CMatrix.hpp"

template<template<typename> class MATRIX, typename T>
using RHS_DATATYPE = typename std::conditional<
    std::is_same_v<MATRIX<T>, RMatrixT<T>> || std::is_same_v<MATRIX<T>, SymRMatrixT<T>>
        || std::is_same_v<MATRIX<T>, AntisymRMatrixT<T>>,
    T,
    typename std::conditional<std::is_same_v<MATRIX<T>, CMatrixT<T>>
                                  || std::is_same_v<MATRIX<T>, HermCMatrixT<T>>,
                              COMPLEX_TYPE<T>, void>::type>::type;

template<template<typename> class MATRIX, typename T>
using RHS_VECTORTYPE = typename std::conditional<
    std::is_same_v<MATRIX<T>, RMatrixT<T>> || std::is_same_v<MATRIX<T>, SymRMatrixT<T>>
        || std::is_same_v<MATRIX<T>, AntisymRMatrixT<T>>,
    RVectorT<T>,
    typename std::conditional<std::is_same_v<MATRIX<T>, CMatrixT<T>>
                                  || std::is_same_v<MATRIX<T>, HermCMatrixT<T>>,
                              CVectorT<T>, void>::type>::type;

template<template<typename> class MATRIX, typename T>
using RHS_MATRIXTYPE = typename std::conditional<
    std::is_same_v<MATRIX<T>, RMatrixT<T>> || std::is_same_v<MATRIX<T>, SymRMatrixT<T>>
        || std::is_same_v<MATRIX<T>, AntisymRMatrixT<T>>,
    RMatrixT<T>,
    typename std::conditional<std::is_same_v<MATRIX<T>, CMatrixT<T>>
                                  || std::is_same_v<MATRIX<T>, HermCMatrixT<T>>,
                              CMatrixT<T>, void>::type>::type;

template<template<typename> class MATRIX, typename T>
class BasicLinearSolverBaseT {
protected:
  const MATRIX<T> &A;
public:
  BasicLinearSolverBaseT(const MATRIX<T> &A) : A(A) {}

  /// Solves the system A * x = rhs
  virtual RHS_VECTORTYPE<MATRIX, T> Solve(const RHS_VECTORTYPE<MATRIX, T> &rhs) = 0;

  /// Solves the system A * x = RHS, i.e., A * x = r is solved for each column of RHS.
  /// The results are returned in the columns of a matrix of the same size
  virtual RHS_MATRIXTYPE<MATRIX, T> Solve(const RHS_MATRIXTYPE<MATRIX, T> &RHS) = 0;
};

template<template<typename> class MATRIX, typename T>
class DirectLinearSolverT : public BasicLinearSolverBaseT<MATRIX, T> {
  MATRIX<T> Ainv;
public:
  DirectLinearSolverT(const MATRIX<T> &A) : BasicLinearSolverBaseT<MATRIX, T>(A), Ainv(A) {
    Ainv.Invert();
  }

  RHS_VECTORTYPE<MATRIX, T> Solve(const RHS_VECTORTYPE<MATRIX, T> &rhs) override {
    return Ainv * rhs;
  }

  RHS_MATRIXTYPE<MATRIX, T> Solve(const RHS_MATRIXTYPE<MATRIX, T> &RHS) override {
    RHS_MATRIXTYPE<MATRIX, T> sol(RHS.rows(), RHS.cols());
    for (int j = 0; j < RHS.cols(); ++j) {
      RHS_VECTORTYPE<MATRIX, T> sol_j = Solve(RHS.col(j));
      for (int i = 0; i < RHS.rows(); ++i)
        sol(i, j) = sol_j[i];
    }
    return sol;
  }
};

void LAPACKSolve(const RMatrix &A, std::vector<double> &b, int Nb);

void LAPACKSolve(const SymRMatrix &A, std::vector<double> &b, int Nb);

void LAPACKSolve(const AntisymRMatrix &A, std::vector<double> &b, int Nb);

void LAPACKSolve(const CMatrix &A, std::vector<std::complex<double>> &b, int Nb);

void LAPACKSolve(const HermCMatrix &A, std::vector<std::complex<double>> &b, int Nb);

template<template<typename> class MATRIX, typename T>
class LAPACKLinearSolverT : public BasicLinearSolverBaseT<MATRIX, T> {
public:
  LAPACKLinearSolverT(const MATRIX<T> &A) : BasicLinearSolverBaseT<MATRIX, T>(A) {}

  RHS_VECTORTYPE<MATRIX, T> Solve(const RHS_VECTORTYPE<MATRIX, T> &rhs) override {
    RHS_VECTORTYPE<MATRIX, T> sol(rhs);
    LAPACKSolve(this->A, sol.asVector(), 1);
    return sol;
  }

  RHS_MATRIXTYPE<MATRIX, T> Solve(const RHS_MATRIXTYPE<MATRIX, T> &RHS) override {
    std::vector<RHS_DATATYPE<MATRIX, T>> b(RHS.rows() * RHS.cols());
    for (int i = 0; i < RHS.rows(); ++i)
      for (int j = 0; j < RHS.cols(); ++j)
        b[j * RHS.rows() + i] = RHS(i, j);
    LAPACKSolve(this->A, b, RHS.cols());
    RHS_MATRIXTYPE<MATRIX, T> SOL(RHS.rows(), RHS.cols());
    for (int i = 0; i < RHS.rows(); ++i)
      for (int j = 0; j < RHS.cols(); ++j)
        SOL(i, j) = b[j * RHS.rows() + i];
    return SOL;
  }
};

template<template<typename> class MATRIX, typename T>
class BasicLinearSolverT {
protected:
  BasicLinearSolverBaseT<MATRIX, T> *solver = nullptr;
public:
  BasicLinearSolverT(const MATRIX<T> &A, const std::string &name = "Direct") {
    if (name == "Direct") solver = new DirectLinearSolverT<MATRIX, T>(A);
    else if (name == "LAPACK" || name == "Lapack") solver = new LAPACKLinearSolverT<MATRIX, T>(A);
    else THROW("BasicLinearSolver " + name + " not implemented")
  }

  RHS_VECTORTYPE<MATRIX, T> Solve(const RHS_VECTORTYPE<MATRIX, T> &rhs) {
    return solver->Solve(rhs);
  }

  RHS_MATRIXTYPE<MATRIX, T> Solve(const RHS_MATRIXTYPE<MATRIX, T> &rhs) {
    return solver->Solve(rhs);
  }

  ~BasicLinearSolverT() {
    if (solver) delete solver;
  }
};

using RBasicLinearSolver = BasicLinearSolverT<RMatrixT, double>;

using SymBasicLinearSolver = BasicLinearSolverT<SymRMatrixT, double>;

using AntisymBasicLinearSolver = BasicLinearSolverT<AntisymRMatrixT, double>;

using CBasicLinearSolver = BasicLinearSolverT<CMatrixT, double>;

using HermBasicLinearSolver = BasicLinearSolverT<HermCMatrixT, double>;

#endif // BASICOPERATOR_HPP
