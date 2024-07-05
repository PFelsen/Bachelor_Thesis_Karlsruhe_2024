#include "LagrangeEllipticAssemble.hpp"
#include "ScalarElement.hpp"


const char *LagrangeEllipticAssemble::Name() const {
  return "LagrangeEllipticAssemble";
}

void LagrangeEllipticAssemble::Initialize(Vector &u) const {

  elementPool.Initialize(u);
  faceElementPool.Initialize(u);

  u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RowBndValues u_c(u, *c);
    if (!u_c.onBnd()) continue;

    const ScalarElement &elem = elementPool.Get(u, *c);
  
    for (int face = 0; face < c.Faces(); ++face) {
      if (u_c.bc(face) == 1) {
        for (int j = 0; j < u.NumberOfNodalPointsOnFace(*c, face); ++j) {
          int k = u.IdOfNodalPointOnFace(*c, face, j);
          u_c(k) = problem.Solution(elem.NodalPoint(k));
          u_c.D(k) = true;
        }
      }
    }
  }
  u.DirichletConsistent();
}

void LagrangeEllipticAssemble::Residual(const cell &c, const Vector &u, Vector &r) const {
  
  const ScalarElement &elem = elementPool.Get(u, *c);

  RowBndValues r_c(r, *c);
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    const Point &z = elem.QPoint(q);
    Tensor K = problem.Permeability(z);
    Scalar f = problem.Load(z, c);
    VectorField DU = elem.Derivative(q, u);
    VectorField K_DU = K * DU;
    for (int i = 0; i < elem.size(); ++i) {
      Scalar U_i = elem.Value(q, i);
      VectorField DU_i = elem.Derivative(q, i);
      r_c(i) += w * (K_DU * DU_i - f * U_i);
    }
  }
  if (!r_c.onBnd()) return;
  for (int f = 0; f < c.Faces(); ++f) {
    int bnd = r_c.bc(f);
    if (bnd != 2) continue;

    ScalarFaceElement faceElem(u, *c, f);

    RowValues r_f(r, faceElem);
    for (int q = 0; q < faceElem.nQ(); ++q) {
      double w = faceElem.QWeight(q);
      Scalar F = problem.Flux(faceElem.QPoint(q)) * faceElem.QNormal(q);
      for (int i = 0; i < faceElem.size(); ++i) {
        Scalar U_i = faceElem.Value(q, i);
        r_f(i) -= w * U_i * F;
      }
    }
  }
}

void LagrangeEllipticAssemble::Jacobi(const cell &c, const Vector &u, Matrix &A) const {
  
  const ScalarElement &elem = elementPool.Get(u, *c);

  RowEntries A_c(A, elem);
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    Tensor K = problem.Permeability(elem.QPoint(q));
    for (unsigned int i = 0; i < elem.size(); ++i) {
      VectorField Dv = K * elem.Derivative(q, i);
      for (unsigned int j = 0; j < elem.size(); ++j) {
        VectorField Dw = elem.Derivative(q, j);
        A_c(i, j) += w * Dv * Dw;
      }
    }
  }
}

FluxPair LagrangeEllipticAssemble::InflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] < 0) continue;

      const ScalarFaceElement &faceElem = faceElementPool.Get(u, *c, face);

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Tensor K = problem.Permeability(faceElem.QPoint(q));
        VectorField Q = K * faceElem.Derivative(q, u);
        Scalar F = Q * faceElem.QNormal(q);
        if (F > 0) { outflow += w * F; }
        else { inflow += w * F; }
      }
    }
  }
  return {PPM->SumOnCommSplit(inflow, u.CommSplit()),
          PPM->SumOnCommSplit(outflow, u.CommSplit())};
}

FluxPair LagrangeEllipticAssemble::PrescribedInflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] < 2) continue;

      const ScalarFaceElement &faceElem = faceElementPool.Get(u, *c, face);     

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Scalar F = problem.Flux(faceElem.QPoint(q)) *
                   faceElem.QNormal(q);
        if (F > 0) { outflow += w * F; }
        else { inflow += w * F; }
      }
    }
  }
  return {PPM->SumOnCommSplit(inflow, u.CommSplit()),
          PPM->SumOnCommSplit(outflow, u.CommSplit())};
}

FluxPair LagrangeEllipticAssemble::OutflowLeftRight(const Vector &u) const {
  double outflowLeft = 0;
  double outflowRight = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] != 1) continue;

      const ScalarFaceElement &faceElem = faceElementPool.Get(u, *c, face); 

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Tensor K = problem.Permeability(faceElem.QPoint(q));
        VectorField Q = K * faceElem.Derivative(q, u);
        Scalar F = Q * faceElem.QNormal(q);
        if (F > 0 && c()[0] <= 1) { outflowLeft += w * F; }
        else if (F > 0 && c()[0] >= 2) outflowRight += w * F;
      }
    }
  }
  return {PPM->SumOnCommSplit(outflowLeft, u.CommSplit()),
          PPM->SumOnCommSplit(outflowRight, u.CommSplit())};
}

double LagrangeEllipticAssemble::GoalFunctional(const Vector &u) const {
  double goal = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] != 1) continue;

      const ScalarFaceElement &faceElem = faceElementPool.Get(u, *c, face);  

      for (int q = 0; q < faceElem.nQ(); ++q) {
        Point z = faceElem.QPoint(q);
        if (z[0] < 0.25) continue;
        if (z[0] > 0.5) continue;
        double w = faceElem.QWeight(q);
        Tensor K = problem.Permeability(faceElem.QPoint(q));
        VectorField Q = K * faceElem.Derivative(q, u);
        Scalar F = Q * faceElem.QNormal(q);
        goal += w * F;
      }
    }
  }
  return PPM->SumOnCommSplit(goal, u.CommSplit());
}

void LagrangeEllipticAssemble::SetExactSolution(Vector &uEx) const {
  for (row r = uEx.rows(); r != uEx.rows_end(); ++r) {
    uEx(r, 0) = problem.Solution(r());
  }
}

void LagrangeEllipticAssemble::SetFlux(const Vector &u, Vector &flux) {
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const ScalarElement &elem = elementPool.Get(u, *c);

    VectorField F = zero;
    double area = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Tensor K = problem.Permeability(elem.QPoint(q));
      VectorField Q = K * elem.Derivative(q, u);
      F += w * Q;
      area += w;
    }
    F *= (1 / area);
    for (int d = 0; d < SpaceDimension; ++d) {
      flux(c(), d) = F[d];
    }
  }
}