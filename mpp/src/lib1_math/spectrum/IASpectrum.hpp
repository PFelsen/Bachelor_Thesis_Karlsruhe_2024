#ifndef IASPECTRUM_H
#define IASPECTRUM_H

#include "CMatrix.hpp"
#include "IAEigenvalueBounds.hpp"
#include "Spectrum.hpp"


/**
 * Encloses the eigenvalues of the matrix eigenvalue problem A * x = lambda * x,
 * where A is a real symmetric interval matrix.
 */
void verifiedEVreal(const IASymRMatrix &A, IAEigenvalueEnclosures &Lambda);

/**
 * Computes lower eigenvalue bounds of the matrix eigenvalue problem A * x = lambda * x,
 * where A is a real symmetric interval matrix.
 */
void verifiedEVreal(const IASymRMatrix &A, IALowerEigenvalueBounds &lambda);

/**
 * Computes upper eigenvalue bounds of the matrix eigenvalue problem A * x = lambda * x,
 * where A is a real symmetric interval matrix.
 */
void verifiedEVreal(const IASymRMatrix &A, IAUpperEigenvalueBounds &Lambda);

/**
 * Encloses the eigenvalues of the matrix eigenvalue problem A * x = lambda * B * x,
 * where A is a real symmetric interval matrix and B is a real positive definite interval matrix.
 */
void verifiedEVreal(const IASymRMatrix &A, const IASymRMatrix &B, IAEigenvalueEnclosures &Lambda);

/**
 * Computes lower eigenvalue bounds of the matrix eigenvalue problem A * x = lambda * B * x,
 * where A is a real symmetric interval matrix and B is a real positive definite interval matrix.
 */
void verifiedEVreal(const IASymRMatrix &A, const IASymRMatrix &B, IALowerEigenvalueBounds &lambda);

/**
 * Computes upper eigenvalue bounds of the matrix eigenvalue problem A * x = lambda * B * x,
 * where A is a real symmetric interval matrix and B is a real positive definite interval matrix.
 */
void verifiedEVreal(const IASymRMatrix &A, const IASymRMatrix &B, IAUpperEigenvalueBounds &Lambda);

/**
 * Encloses the eigenvalues of the matrix eigenvalue problem A * x = lambda * x,
 * where A is a complex hermitian interval matrix.
 */
void verifiedEVcomplex(const IAHermCMatrix &A, IAEigenvalueEnclosures &Lambda);

/**
 * Computes lower eigenvalue bounds of the matrix eigenvalue problem A * x = lambda * x,
 * where A is a complex hermitian interval matrix.
 */
void verifiedEVcomplex(const IAHermCMatrix &A, IALowerEigenvalueBounds &lambda);

/**
 * Computes upper eigenvalue bounds of the matrix eigenvalue problem A * x = lambda * x,
 * where A is a complex hermitian interval matrix.
 */
void verifiedEVcomplex(const IAHermCMatrix &A, IAUpperEigenvalueBounds &Lambda);

/**
 * Encloses the eigenvalues of the matrix eigenvalue problem A * x = lambda * B * x,
 * where A is a complex hermitian interval matrix and B is a complex positive definite interval
 * matrix.
 */
void verifiedEVcomplex(const IAHermCMatrix &A, const IAHermCMatrix &B,
                       IAEigenvalueEnclosures &Lambda);

/**
 * Computes lower eigenvalue bounds of the matrix eigenvalue problem A * x = lambda * B * x,
 * where A is a complex hermitian interval matrix and B is a complex positive definite interval
 * matrix.
 */
void verifiedEVcomplex(const IAHermCMatrix &A, const IAHermCMatrix &B,
                       IALowerEigenvalueBounds &lambda);

/**
 * Computes upper eigenvalue bounds of the matrix eigenvalue problem A * x = lambda * B * x,
 * where A is a complex hermitian interval matrix and B is a complex positive definite interval
 * matrix.
 */
void verifiedEVcomplex(const IAHermCMatrix &A, const IAHermCMatrix &B,
                       IAUpperEigenvalueBounds &Lambda);

#endif // EIGENVALUEMETHODS_H
