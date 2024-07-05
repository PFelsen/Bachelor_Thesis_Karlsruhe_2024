#ifndef SPECTRUM_H
#define SPECTRUM_H

#include "CMatrix.hpp"

#include <memory>

typedef double Eigenvalue;

typedef RVector Eigenvalues;

typedef RVector Eigenvector;

typedef RMatrix Eigenvectors;

typedef std::complex<double> CEigenvalue;

typedef CVector CEigenvalues;

typedef CVector CEigenvector;

typedef CMatrix CEigenvectors;

/**
 * Solves the matrix eigenvalue problem A * x = lambda * x, where A is a real symmetric matrix.
 * Here the eigen vectors are calculated and stored in x
 */
void EVreal(const SymRMatrix &A, Eigenvalues &lambda, Eigenvectors &x);

/**
 * Solves the matrix eigenvalue problem A * x = lambda * x, where A is a real symmetric matrix.
 * Here no eigen vectors are calculated
 */
void EVreal(const SymRMatrix &A, Eigenvalues &lambda);

/**
 * Solves the matrix eigenvalue problem A * x = lambda * x, where A is a real matrix.
 * Here the eigen vectors are calculated and stored in x
 */
void EVreal(const RMatrix &A, CEigenvalues &lambda, CEigenvectors &x);

/**
 * Solves the matrix eigenvalue problem A * x = lambda * x, where A is a real matrix.
 * Here no eigen vectors are calculated
 */
void EVreal(const RMatrix &A, CEigenvalues &lambda);

/**
 * Solves the matrix eigenvalue problem A * x = lambda * B * x, where A is a real symmetric matrix
 * and B a real positive definit matrix. Here the eigen vectors are calculated and stored in x
 */
void EVreal(const SymRMatrix &A, const SymRMatrix &B, Eigenvalues &lambda, Eigenvectors &x);

/**
 * Solves the matrix eigenvalue problem A * x = lambda * B * x, where A is a real symmetric matrix
 * and B a real positive definit matrix. Here no eigen vectors are calculated
 */
void EVreal(const SymRMatrix &A, const SymRMatrix &B, Eigenvalues &lambda);

/**
 * Solves the matrix eigenvalue problem A * x = lambda * x, where A is a complex hermitian matrix.
 * Here the eigen vectors are calculated and stored in x
 */
void EVcomplex(const HermCMatrix &A, Eigenvalues &lambda, CEigenvectors &x);

/**
 * Solves the matrix eigenvalue problem A * x = lambda * x, where A is a complex hermitian matrix.
 * Here no eigen vectors are calculated
 */
void EVcomplex(const HermCMatrix &A, Eigenvalues &lambda);

/**
 * Solves the matrix eigenvalue problem A * x = lambda * x, where A is a complex matrix.
 * Here the eigen vectors are calculated and stored in x
 */
void EVcomplex(const CMatrix &A, CEigenvalues &lambda, CEigenvectors &x);

/**
 * Solves the matrix eigenvalue problem A * x = lambda * x, where A is a complex matrix.
 * Here no eigen vectors are calculated
 */
void EVcomplex(const CMatrix &A, CEigenvalues &lambda);

/**
 * Solves the matrix eigenvalue problem A * x = lambda * B * x, where A is a complex hermitian
 * matrix and B a complex positive definit matrix. Here the eigen vectors are calculated and stored
 * in x
 */
void EVcomplex(const HermCMatrix &A, const HermCMatrix &B, Eigenvalues &lambda, CEigenvectors &x);

/**
 * Solves the matrix eigenvalue problem A * x = lambda * B * x, where A is a complex hermitian
 * matrix and B a complex positive definit matrix. Here the eigen vectors are calculated and stored
 * in x
 */
void EVcomplex(const HermCMatrix &A, const HermCMatrix &B, Eigenvalues &lambda);

class Spectrum {
  bool computeEVec;
  Eigenvalues lambda;
  std::unique_ptr<Eigenvectors> e1 = nullptr;
  std::unique_ptr<CEigenvectors> e2 = nullptr;
public:
  explicit Spectrum(const SymRMatrix &A, bool computeEVec = false);

  explicit Spectrum(const HermCMatrix &A, bool computeEVec = false);

  Spectrum(const SymRMatrix &A, const SymRMatrix &B, bool computeEVec = false);

  Spectrum(const HermCMatrix &A, const HermCMatrix &B, bool computeEVec = false);

  Eigenvalue operator()(int i) const { return lambda[i]; }

  const Eigenvalues &operator()() const { return lambda; }

  Eigenvectors &getEigenvectors() const;

  CEigenvectors &getCEigenvectors() const;

  Eigenvalue operator()(int i, Eigenvector &w) const;

  Eigenvalue operator()(int i, CEigenvector &w) const;

  Eigenvalue min() const;

  double absmin() const;

  Eigenvalue max() const;

  double absmax() const;

  double cond() const;

  friend std::ostream &operator<<(std::ostream &s, const Spectrum &S) {
    return s << " min " << S.min() << " max " << S.max() << " absmin " << S.absmin() << " absmax "
             << S.absmax() << " cond " << S.cond() << endl;
  }
};

#endif // SPECTRUM_H
