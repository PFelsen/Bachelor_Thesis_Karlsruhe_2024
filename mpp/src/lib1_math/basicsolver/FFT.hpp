#ifndef FFT_HPP
#define FFT_HPP

#include "CMatrix.hpp"
#include "CTensor.hpp"
#include "CVector.hpp"
#include "RMatrix.hpp"
#include "RTensor.hpp"
#include "RVector.hpp"

/* The whole namespace uses the most general FFT/IFFT, although there are specialized functions
 (r2r,r2c,r2c)
 */
namespace FFT {
void ComplexToComplexVector(const CVector &in, CVector &out);

void InvComplexToComplexVector(const CVector &in, CVector &out);

void RealToRealVector(const RVector &in, RVector &out);

void RealToComplexVector(const RVector &in, CVector &out);

void InvComplexToRealVector(const CVector &in, RVector &out);

void RealToRealMatrix(const RMatrix &in, RMatrix &out);

void ComplexToComplexMatrix(const CMatrix &in, CMatrix &out);

void RealToRealTensor(const RTensor &in, RTensor &out);

void ComplexToComplexTensor(const CTensor &in, CTensor &out);
}; // namespace FFT

#endif // FFT_HPP
