#ifndef ARGYRISELEMENT_HPP
#define ARGYRISELEMENT_HPP

#include "Element.hpp"

/*
 * Based on a paper by Dominguez and Sayas
 * "Algorithm 884: A Simple Matlab Implementation of the Argyris Element"
 */

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class ArgyrisBasicElementT : public ElementT<TT, sDim, tDim> {
  template<typename TT1, int sDim1, int tDim1>
  friend class BernsteinRangesT;
public:
  // corners contains the geometrical corners of domain (not the corners of all cells!)
  // Here all DoFs need to be set to zero
  static void H20BC(const std::list<Point> &corners, Vector &u);

  static void H20BC(const std::list<Point> &corners, Vector &u, const Cell &C);
protected:
  PointT<TT, sDim, tDim> normal[3]; // outer normal
  TT sign[3];

  TensorT<TT, 3> theta;
  TT C0[2][2];
  TT C1[3][3];
  TT C2[2][6];
  TT C3[2][6];
  TT C4[2][6];
  TT C5[3];

  ArgyrisBasicElementT(const VectorMatrixBase &base, const Cell &c);

  ArgyrisBasicElementT(const VectorMatrixBase &base,
                       const BasicElementT<TT, sDim, tDim> &baseElement);

  std::vector<TT> values(const PointT<TT, sDim, tDim> &z) const;

  std::vector<VectorFieldT<TT, sDim>> derivatives(const PointT<TT, sDim, tDim> &z) const;

  std::vector<SymTensorT<TT, sDim>> hessians(const PointT<TT, sDim, tDim> &z) const;
public:
  PointT<TT, sDim, tDim> OrientedNormal(int i) const { return sign[i] * normal[i]; }

  PointT<TT, sDim, tDim> OuterNormal(int i) const { return normal[i]; }

  int get_maxk(int i) const;

  int indexing_k(int i, int k, int m = 0) const;

  TT getC(int i, int j) const;

  TT getTheta(int i, int j) const { return theta[i][j]; }

  // TODO: Should be implemented in Transformation
  PointT<TT, sDim, tDim> GlobalToLocal(const PointT<TT, sDim, tDim> &global) const {
    TensorT<TT, sDim> R;
    R[0][0] = TT((this->GetCell())[1][0]) - (this->GetCell())[0][0];
    R[0][1] = TT((this->GetCell())[2][0]) - (this->GetCell())[0][0];
    R[1][0] = TT((this->GetCell())[1][1]) - (this->GetCell())[0][1];
    R[1][1] = TT((this->GetCell())[2][1]) - (this->GetCell())[0][1];
    return R.Invert() * (global - (this->GetCell())[0]);
  }
protected:
  void init();

  SymTensorT<TT, sDim> multiplyTheta(const SymTensorT<TT, sDim> &S) const;

  int indexingShape(int i, int k) const;
};

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class ArgyrisElementT : public ArgyrisBasicElementT<TT, sDim, tDim> {
protected:
  ShapeValues<TT> value;
  ShapeValues<VectorFieldT<TT, sDim>> gradient;
  ShapeValues<SymTensorT<TT, sDim>> hessian;
public:
  ArgyrisElementT(const VectorMatrixBase &base, const Cell &c);

  ArgyrisElementT(const VectorMatrixBase &base, const BasicElementT<TT, sDim, tDim> &baseElement);

  TT Value(int q, int i, int k) const { return Value(q, ShapeId{this->indexingShape(i, k)}); }

  TT Value(int q, const ShapeId &iter) const { return value[q][iter.id]; }

  TT Value(const PointT<TT, sDim, tDim> &z, int i, int k) const;

  TT Value(const PointT<TT, sDim, tDim> &z, const ShapeId &iter) const;

  TT Value(int q, const Vector &u, int m = 0) const;

  TT Value(const PointT<TT, sDim, tDim> &z, const Vector &u, int m = 0) const;

  const VectorFieldT<TT, sDim> &Derivative(int q, int i, int k) const {
    return Derivative(q, ShapeId{this->indexingShape(i, k)});
  }

  const VectorFieldT<TT, sDim> &Derivative(int q, const ShapeId &iter) const {
    return gradient[q][iter.id];
  }

  VectorFieldT<TT, sDim> Derivative(const PointT<TT, sDim, tDim> &z, int i, int k) const;

  VectorFieldT<TT, sDim> Derivative(const PointT<TT, sDim, tDim> &z, const ShapeId &iter) const;

  VectorFieldT<TT, sDim> Derivative(int q, const Vector &u, int m = 0) const;

  VectorFieldT<TT, sDim> Derivative(const PointT<TT, sDim, tDim> &z, const Vector &u,
                                    int m = 0) const;

  const SymTensorT<TT, sDim> &Hessian(int q, int i, int k) const {
    return Hessian(q, ShapeId{this->indexingShape(i, k)});
  }

  const SymTensorT<TT, sDim> &Hessian(int q, const ShapeId &iter) const {
    return hessian[q][iter.id];
  }

  SymTensorT<TT, sDim> Hessian(const PointT<TT, sDim, tDim> &z, int i, int k) const;

  SymTensorT<TT, sDim> Hessian(const PointT<TT, sDim, tDim> &z, const ShapeId &iter) const;

  SymTensorT<TT, sDim> Hessian(int q, const Vector &u, int m = 0) const;

  SymTensorT<TT, sDim> Hessian(const PointT<TT, sDim, tDim> &z, const Vector &u, int m = 0) const;

  TT Laplace(int q, int i, int k) const;

  TT Laplace(int q, const ShapeId &iter) const;

  TT Laplace(const PointT<TT, sDim, tDim> &z, int i, int k) const;

  TT Laplace(const PointT<TT, sDim, tDim> &z, const ShapeId &iter) const;

  TT Laplace(int q, const Vector &u, int m = 0) const;

  TT Laplace(const PointT<TT, sDim, tDim> &z, const Vector &u, int m = 0) const;

  friend std::ostream &operator<<(std::ostream &s, const ArgyrisElementT<TT, sDim, tDim> &E) {
    return s << E.shape << endl;
  }

  friend class TestArgyrisElement;
};

using ArgyrisElement = ArgyrisElementT<>;

#ifdef BUILD_IA

using IAArgyrisElement = ArgyrisElementT<IAInterval, SpaceDimension, TimeDimension>;

#endif

template<typename C, class D>
void applyC(const C (&C0)[2][2], const C (&C1)[3][3], const C (&C2)[2][6], const C (&C3)[2][6],
            const C (&C4)[2][6], const C (&C5)[3], D &c, const vector<D> &v, int i) {
  switch (i) {
  case 0:
    c = C2[0][0] * v[18] + C2[1][0] * v[20] + v[0];
    break;
  case 1:
    c = C2[0][1] * v[18] + C2[1][1] * v[20] + C0[0][0] * v[1] + C0[1][0] * v[2];
    break;
  case 2:
    c = C2[0][2] * v[18] + C2[1][2] * v[20] + C0[0][1] * v[1] + C0[1][1] * v[2];
    break;
  case 3:
    c = C2[0][3] * v[18] + C2[1][3] * v[20] + C1[0][0] * v[3] + C1[1][0] * v[4] + C1[2][0] * v[5];
    break;
  case 4:
    c = C2[0][4] * v[18] + C2[1][4] * v[20] + C1[0][1] * v[3] + C1[1][1] * v[4] + C1[2][1] * v[5];
    break;
  case 5:
    c = C2[0][5] * v[18] + C2[1][5] * v[20] + C1[0][2] * v[3] + C1[1][2] * v[4] + C1[2][2] * v[5];
    break;

  case 6:
    c = C3[0][0] * v[18] + C3[1][0] * v[19] + v[6];
    break;
  case 7:
    c = C3[0][1] * v[18] + C3[1][1] * v[19] + C0[0][0] * v[7] + C0[1][0] * v[8];
    break;
  case 8:
    c = C3[0][2] * v[18] + C3[1][2] * v[19] + C0[0][1] * v[7] + C0[1][1] * v[8];
    break;
  case 9:
    c = C3[0][3] * v[18] + C3[1][3] * v[19] + C1[0][0] * v[9] + C1[1][0] * v[10] + C1[2][0] * v[11];
    break;
  case 10:
    c = C3[0][4] * v[18] + C3[1][4] * v[19] + C1[0][1] * v[9] + C1[1][1] * v[10] + C1[2][1] * v[11];
    break;
  case 11:
    c = C3[0][5] * v[18] + C3[1][5] * v[19] + C1[0][2] * v[9] + C1[1][2] * v[10] + C1[2][2] * v[11];
    break;

  case 12:
    c = C4[0][0] * v[19] + C4[1][0] * v[20] + v[12];
    break;
  case 13:
    c = C4[0][1] * v[19] + C4[1][1] * v[20] + C0[0][0] * v[13] + C0[1][0] * v[14];
    break;
  case 14:
    c = C4[0][2] * v[19] + C4[1][2] * v[20] + C0[0][1] * v[13] + C0[1][1] * v[14];
    break;
  case 15:
    c = C4[0][3] * v[19] + C4[1][3] * v[20] + C1[0][0] * v[15] + C1[1][0] * v[16]
        + C1[2][0] * v[17];
    break;
  case 16:
    c = C4[0][4] * v[19] + C4[1][4] * v[20] + C1[0][1] * v[15] + C1[1][1] * v[16]
        + C1[2][1] * v[17];
    break;
  case 17:
    c = C4[0][5] * v[19] + C4[1][5] * v[20] + C1[0][2] * v[15] + C1[1][2] * v[16]
        + C1[2][2] * v[17];
    break;

  case 18:
    c = C5[0] * v[18];
    break;
  case 19:
    c = C5[1] * v[19];
    break;
  case 20:
    c = C5[2] * v[20];
    break;
  }
}

#ifdef BUILD_IA

using BernsteinCoeff = std::vector<IAInterval>;

template<typename T>
struct BernsteinDataT {
  T lambda;
  const Vector &vec;
};

using BernsteinData = BernsteinDataT<double>;

#ifdef BUILD_IA

using IABernsteinData = BernsteinDataT<IAInterval>;

#endif

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class BernsteinRangesT {
protected:
  // coefficients of fem solution in Bernstein basis
  BernsteinCoeff coeff{};
  // coefficients of gradient in Bernstein basis
  BernsteinCoeff coeffDx{};
  BernsteinCoeff coeffDy{};
  // coefficients of hessian in Bernstein basis
  BernsteinCoeff coeffDxx{};
  BernsteinCoeff coeffDxy{};
  BernsteinCoeff coeffDyy{};
public:
  BernsteinRangesT(const ArgyrisBasicElementT<TT, sDim, tDim> &element,
                   std::vector<BernsteinDataT<TT>> data);

  BernsteinRangesT(const ArgyrisBasicElementT<TT, sDim, tDim> &element, const Vector &vec);

  /// Encloses the range of the fem solution
  void Range(IAInterval &range) const;

  /// Encloses the range of the gradiant of the fem solution
  void RangeGradient(IAInterval &rangeDx, IAInterval &rangeDy) const;

  /// Encloses the range of the gradiant of the fem solution
  void RangeHessian(IAInterval &rangeDxx, IAInterval &rangeDxy, IAInterval &rangeDyy) const;

  /// Computes the abs max of the fem solution
  void AbsMax(double &absMax) const;

  /// Computes the abs max of the gradiant of the fem solution
  void AbsMaxGradient(double &absMaxDx, double &absMaxDy) const;

  /// Computes the abs max of the hessian of the fem solution
  void AbsMaxHessian(double &absMaxDxx, double &absMaxDxy, double &absMaxDyy) const;
};

using BernsteinRanges = BernsteinRangesT<>;

using IABernsteinRanges = BernsteinRangesT<IAInterval, SpaceDimension, TimeDimension>;

#endif // BUILD_IA

#endif // ARGYRISELEMENT_HPP
