#ifndef WCELEMENT_HPP
#define WCELEMENT_HPP

#include "Element.hpp"
#include "RVector.hpp"

class WCLaplaceElement : public Element {
  const Cell &c; // TODO: remove and use cell in Element!!
  const Shape *S_U;
  vector<vector<Scalar>> value_U;
  vector<vector<VectorField>> grad_U;
  vector<short> face_dofs[3][6];
  short nr_face_dofs[3][6];
  short dofs[6];
  int nState_;
  bool active[6];

  vector<Point> nodalPoints;
public:
  WCLaplaceElement(const VectorMatrixBase &base, const Cell &C) : Element(base, C), c(C) {
    const WCDiscretization &D = dynamic_cast<const WCDiscretization &>(base.GetDisc());

    S_U->NodalPoints(c, nodalPoints);
    nState_ = 0;
    int max_nodalpoints = 0;
    S_U = D.GetShape2(C, D.get_cell_deg_WC(C(), 0));
    nState_ += S_U->size();
    if (S_U->size() > max_nodalpoints) max_nodalpoints = S_U->size();
    value_U.resize(D.GetQuad(C).size(), vector<Scalar>(max_nodalpoints, 0.0));
    grad_U.resize(D.GetQuad(C).size(), vector<VectorField>(max_nodalpoints, zero));
    for (int f = 0; f < c.Faces(); ++f) {
      face_dofs[0][f] = D.get_face_dofs(c.Face(f), 0);
      face_dofs[1][f] = D.get_face_dofs(c.Face(f), 1);
      face_dofs[2][f] = D.get_face_dofs(c.Face(f), 2);
      nr_face_dofs[0][f] = face_dofs[0][f].size();
      nr_face_dofs[1][f] = face_dofs[1][f].size();
      nr_face_dofs[2][f] = face_dofs[2][f].size();
      dofs[f] = nr_face_dofs[0][f] + nr_face_dofs[1][f];
    }
    for (int q = 0; q < this->nQ(); ++q) {
      for (int i = 0; i < S_U->size(); ++i) {
        value_U[q][i] = (*(S_U))(this->LocalQPoint(q), i);
        grad_U[q][i] = GetTransformation(q) * S_U->LocalGradient(this->LocalQPoint(q), i);
      }
    }
  }

  int NodalPoints() const override { return nodalPoints.size(); }

  const Point &NodalPoint(int i) const override { return nodalPoints[i]; }

  Point GlobalToLocal(const Point &z) const {
    if (c.Corners() == 4) {
      if (c.dim() != 2) THROW("Not implemented!")
      Tensor T;
      T[0][0] = c[1][0] - c[0][0];
      T[0][1] = c[3][0] - c[0][0];
      T[1][0] = c[1][1] - c[0][1];
      T[1][1] = c[3][1] - c[0][1];
      if (c.dim() == 2) {
        T[0][2] = 0;
        T[1][2] = 0;
        T[2][0] = 0;
        T[2][1] = 0;
        T[2][2] = 1;
      } else {
        THROW("Not implemented")
      }
      Tensor IT = Invert(T);
      return IT * (z - c[0]);
    } else if (c.Corners() == 3) {
      Tensor T;
      T[0][0] = c[1][0] - c[0][0];
      T[0][1] = c[2][0] - c[0][0];
      T[1][0] = c[1][1] - c[0][1];
      T[1][1] = c[2][1] - c[0][1];
      if (c.dim() == 2) {
        T[0][2] = 0;
        T[1][2] = 0;
        T[2][0] = 0;
        T[2][1] = 0;
        T[2][2] = 1;
      } else {
        T[0][2] = c[3][0] - c[0][0];
        T[1][2] = c[3][1] - c[0][1];
        T[2][0] = c[1][2] - c[0][2];
        T[2][1] = c[2][2] - c[0][2];
        T[2][2] = c[3][2] - c[0][2];
      }
      Tensor IT = Invert(T);
      return IT * (z - c[0]);
    } else if (c.Corners() == 8) {
      Tensor T;
      T[0][0] = c[1][0] - c[0][0];
      T[0][1] = c[3][0] - c[0][0];
      T[0][2] = c[4][0] - c[0][0];
      T[1][0] = c[1][1] - c[0][1];
      T[1][1] = c[3][1] - c[0][1];
      T[1][2] = c[4][1] - c[0][1];
      T[2][0] = c[1][2] - c[0][2];
      T[2][1] = c[3][2] - c[0][2];
      T[2][2] = c[4][2] - c[0][2];
      Tensor IT = Invert(T);
      return IT * (z - c[0]);
    }
    THROW("Not implemented")
  }

  bool NonDirichletBnd(int i) const { return ((Bnd(i) != 1) && (Bnd(i) >= 0)); }

  bool FreeNeumannBnd(int i) const { return (Bnd(i) == 0); }

  bool NonFreeNeumannBnd(int i) const {
    return (Bnd(i) == 2 || Bnd(i) == 1016 || Bnd(i) == 1015 || Bnd(i) == 1201);
  }

  bool NeumannBnd(int i) const { return (Bnd(i) == 0) || (Bnd(i) == 2); }

  bool DirichletBnd(int i) const { return (Bnd(i) == 1 || Bnd(i) == 1000); }

  bool ContactBnd(int i) const { return (Bnd(i) == 3); }

  int nState() const { return nState_; }

  int nfTest(int f) const { return dofs[f] + nr_face_dofs[2][f]; }

  int nTest() const {
    int n = 0;
    for (int f = 0; f < c.Faces(); ++f)
      n += nfTest(f);
    return n;
  }

  int nU() const { return S_U->size(); }

  Scalar Value(int q, int i) const { return value_U[q][i]; }

  Scalar Value(int q, const RVector &UC) const {
    Scalar U = 0;
    for (int i = 0; i < nState(); ++i)
      U += UC[i] * Value(q, i);
    return U;
  }

  Scalar Value(const Point &lz, const RVector &UC) const {
    Scalar U = 0;
    for (int i = 0; i < nState(); ++i)
      U += UC[i] * (*(S_U))(lz, i);
    return U;
  }

  VectorField Derivative(const Point &lz, const RVector &UC) const {
    VectorField DU = zero;
    for (int i = 0; i < nState(); ++i)
      DU += UC[i] * (c.GetTransformation(lz) * S_U->LocalGradient(lz, i));
    return DU;
  }

  const VectorField Derivative(int q, int i) const {
    VectorField DU = zero;
    DU = grad_U[q][i];
    return DU;
  }

  VectorField Derivative(int q, const RVector &UC) const {
    VectorField DG = zero;
    for (int i = 0; i < nState(); ++i)
      DG += UC[i] * Derivative(q, i);
    return DG;
  }
};
#endif // WCELEMENT_HPP
