#include "DGEllipticAssemble.hpp"
#include "DGElement.hpp"


const char *DGEllipticAssemble::Name() const {
  return "DGEllipticAssemble";
}

void DGEllipticAssemble::Initialize(Vector &u) const {

  elementPool.Initialize(u);
  faceElementPool.Initialize(u);


//  u.ClearDirichletFlags();
//    for (cell c = u.cells(); c != u.cells_end(); ++c) {
//        RowBndValues u_c(u, *c);
//        if (!u_c.onBnd()) continue;
//        DGElement elem(u, *c);
//       for (int face = 0; face < c.Faces(); ++face) {
//            if (u_c.bc(face) == 1) {
//              for (int j = 0; j < u.NodalPointsOnFace(*c, face); ++j) {
//                    int k = u.NodalPointOnFace(*c, face, j);
//                    u_c(k,0) = problem.Solution(elem.NodalPoint(k));
//                    u_c.D(k,0) = true;
//                }
//            }
//        }
//    }
//    u.DirichletConsistent();
}

void DGEllipticAssemble::Residual(const cell &c, const Vector &u, Vector &r) const {

  const DGElement &elem = elementPool.Get(u, *c);

  RowBndValues r_c(r, *c);
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    const Point &z = elem.QPoint(q);
    Tensor K = problem.Permeability(z);
    Scalar Load = problem.Load(z, c);
    VectorField K_DU = K * elem.Derivative(q, u);
    for (int i = 0; i < elem.NodalPoints(); ++i) {
      Scalar U_i = elem.Value(q, i);
      VectorField DU_i = elem.Derivative(q, i);
      r(c(), i) += w * (K_DU * DU_i - Load * U_i);
    }
  }
  for (int f = 0; f < c.Faces(); ++f) {
    Scalar scaledPenalty = (pow(degree,2)*penalty) / u.GetMesh().MaxMeshWidth();

    const DGFaceElement &faceElem = faceElementPool.Get(u, *c, f);

    if (r.OnBoundary(*c, f)) {
      int bnd = r_c.bc(f);
      if (bnd == 2) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          const Point &z = faceElem.QPoint(q);
          Scalar F = problem.Flux(z) * faceElem.QNormal(q);
          for (int i = 0; i < faceElem.NodalPoints(); ++i) {
            Scalar U_i = faceElem.Value(q, i);
            r(c(), i) -= w * U_i * F;
          }
        }
      } else if (bnd == 1) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          const Point &z = faceElem.QPoint(q);
          const Point &N = faceElem.QNormal(q);
          Tensor K = problem.Permeability(c());
          Scalar U_diff = faceElem.Value(q, u) - problem.Solution(z);
          VectorField K_DU = K * faceElem.Derivative(q, u);
          Scalar F = K_DU * N;
          for (int i = 0; i < faceElem.NodalPoints(); ++i) {
            Scalar phi_i = faceElem.Value(q, i);
            Scalar NDphi_i = (K * faceElem.Derivative(q, i)) * N;
            r(c(), i) -= w * (F * phi_i +
                              U_diff * (NDphi_i * sign - scaledPenalty * phi_i));
          }
        }
      }
    } else {
      cell cf = r.find_neighbour_cell(c, f);
      if (cf() < c()) continue;
      int f1 = r.find_neighbour_face_id(c.Face(f), cf);

      const DGFaceElement &otherFaceElem = faceElementPool.Get(u, *cf, f1);

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        const Point &N = faceElem.QNormal(q);
        Tensor K = problem.Permeability(c());
        Scalar U = faceElem.Value(q, u);
        VectorField K_DU = K * faceElem.Derivative(q, u);
        const Point &Qf_c = faceElem.QPoint(q);
        int q1 = faceElem.findQPointID(otherFaceElem, Qf_c);
        Tensor K_1 = problem.Permeability(cf());
        Scalar U_1 = otherFaceElem.Value(q1, u);
        VectorField K_DU_1 = K_1 * otherFaceElem.Derivative(q1, u);
        for (int i = 0; i < faceElem.NodalPoints(); ++i) {
          Scalar phi_i = faceElem.Value(q, i);
          Scalar NDphi_i = (K * faceElem.Derivative(q, i)) * N;
          r(c(), i) += w * ((scaledPenalty * U - 0.5 * (K_DU * N)) * phi_i
                            - 0.5 * U * NDphi_i * sign);
          r(c(), i) += w * ((-scaledPenalty * U_1 - 0.5 * (K_DU_1 * N)) * phi_i
                            + 0.5 * U_1 * NDphi_i * sign);
        }
        for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
          Scalar phi_j = otherFaceElem.Value(q1, j);
          Scalar NDphi_j = (K_1 * otherFaceElem.Derivative(q1, j)) * N;
          r(cf(), j) += w * ((0.5 * K_DU_1 * N + scaledPenalty * U_1) * phi_j
                             + 0.5 * U_1 * NDphi_j * sign);
          r(cf(), j) += w * ((0.5 * K_DU * N - scaledPenalty * U) * phi_j
                             - 0.5 * U * NDphi_j * sign);
        }
      }
    }
  }
}

void DGEllipticAssemble::Jacobi(const cell &c, const Vector &u, Matrix &A) const {

  const DGElement &elem = elementPool.Get(u, *c);

  DGRowEntries A_c(A, *c, *c);
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    Tensor K = problem.Permeability(elem.QPoint(q));
    for (int i = 0; i < elem.NodalPoints(); ++i) {
      VectorField K_DU_i = K * elem.Derivative(q, i);
      for (int j = 0; j < elem.NodalPoints(); ++j) {
        VectorField DU_j = elem.Derivative(q, j);
        A_c(i, j) += w * K_DU_i * DU_j;
      }
    }
  }
  BFParts bnd(u.GetMesh(), *c);
  for (int f = 0; f < c.Faces(); ++f) {
    Scalar s =  (pow(degree,2)*penalty) / u.GetMesh().MaxMeshWidth();

    const DGFaceElement &faceElem = faceElementPool.Get(u, *c, f);

    if (u.OnBoundary(*c, f)) {
      if (bnd[f] == 1) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          const Point &z = faceElem.QPoint(q);
          const Point &N = faceElem.QNormal(q);
          Tensor K = problem.Permeability(c());
          for (int i = 0; i < faceElem.NodalPoints(); ++i) {
            Scalar phi_i = faceElem.Value(q, i);
            Scalar NDphi_i = (K * faceElem.Derivative(q, i)) * N;
            for (int j = 0; j < faceElem.NodalPoints(); ++j) {
              Scalar phi_j = faceElem.Value(q, j);
              Scalar NDphi_j = (K * faceElem.Derivative(q, j)) * N;
              A_c(i, j) -= w * (NDphi_i * phi_j * sign +
                                phi_i * (NDphi_j - s * phi_j));
            }
          }
        }
      }
    } else {
      cell cf = u.find_neighbour_cell(c, f);
      if (cf() < c()) continue;
      DGRowEntries A_cf(A, *c, *cf);
      DGRowEntries A_fc(A, *cf, *c);
      DGRowEntries A_ff(A, *cf, *cf);
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);

      const DGFaceElement &otherFaceElem = faceElementPool.Get(u, *cf, f1);

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        const Point &N = faceElem.QNormal(q);
        Tensor K = problem.Permeability(c());
        // Tensor K = problem.Permeability(*cf());
        const Point &Qf_c = faceElem.QPoint(q);
        int q1 = otherFaceElem.findQPointID(faceElem, Qf_c);
        Tensor K_1 = problem.Permeability(cf());
        for (int i = 0; i < faceElem.NodalPoints(); ++i) {
          Scalar phi_i = faceElem.Value(q, i);
          Scalar NDphi_i = (K * faceElem.Derivative(q, i)) * N;
          for (int j = 0; j < faceElem.NodalPoints(); ++j) {
            Scalar phi_j = faceElem.Value(q, j);
            Scalar NDphi_j = (K * faceElem.Derivative(q, j)) * N;
            A_c(i, j) += w * (-0.5 * NDphi_i * phi_j * sign
                              - 0.5 * phi_i * NDphi_j
                              + s * phi_i * phi_j);
          }
          for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
            Scalar phi_j = otherFaceElem.Value(q1, j);
            Scalar NDphi_j = (K_1 * otherFaceElem.Derivative(q1, j)) * N;
            A_cf(i, j) += w * (0.5 * NDphi_i * phi_j * sign
                               - phi_i * 0.5 * NDphi_j
                               - s * phi_i * phi_j);
            A_fc(j, i) += w * (0.5 * NDphi_i * phi_j
                               - phi_i * 0.5 * NDphi_j * sign
                               - s * phi_i * phi_j);
          }
        }
        for (int i = 0; i < otherFaceElem.NodalPoints(); ++i) {
          Scalar phi_i = otherFaceElem.Value(q1, i);
          Scalar NDphi_i = (K_1 * otherFaceElem.Derivative(q1, i)) * N;
          for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
            Scalar phi_j = otherFaceElem.Value(q1, j);
            Scalar NDphi_j = (K_1 * otherFaceElem.Derivative(q1, j)) * N;
            A_ff(i, j) += w * (0.5 * NDphi_i * phi_j * sign
                               + phi_i * 0.5 * NDphi_j
                               + s * phi_i * phi_j);
          }
        }
      }
    }
  }
}

FluxPair DGEllipticAssemble::InflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] < 0) continue;

      const DGFaceElement &faceElem = faceElementPool.Get(u, *c, face);

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

FluxPair DGEllipticAssemble::PrescribedInflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] < 2) continue;

      const DGFaceElement &faceElem = faceElementPool.Get(u, *c, face);

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Scalar FN = problem.Flux(faceElem.QPoint(q)) * faceElem.QNormal(q);
        if (FN > 0) { outflow += w * FN; }
        else { inflow += w * FN; }
      }
    }
  }
  return {PPM->SumOnCommSplit(inflow, u.CommSplit()),
          PPM->SumOnCommSplit(outflow, u.CommSplit())};
}

FluxPair DGEllipticAssemble::OutflowLeftRight(const Vector &u) const {
  double outflowLeft = 0;
  double outflowRight = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] != 1) continue;

      const DGFaceElement &faceElem = faceElementPool.Get(u, *c, face);

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

double DGEllipticAssemble::GoalFunctional(const Vector &u) const {
  double goal = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] != 1) continue;

      const DGFaceElement &faceElem = faceElementPool.Get(u, *c, face);

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

void DGEllipticAssemble::SetExactSolution(Vector &u) const {
  u = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    row r = u.find_row(c());

    const DGElement &elem = elementPool.Get(u, *c);

    for (int j = 0; j < elem.NodalPoints(); ++j)
      u(r, j) = problem.Solution(elem.NodalPoint(j));
  }
  u.Accumulate();
}

void DGEllipticAssemble::SetFlux(const Vector &u, Vector &flux) {}

void DGEllipticAssemble::SetPressure(const Vector &u, Vector &p) const {
  p.Clear();
  Vector e(p);
  for (cell c=p.cells(); c!=p.cells_end(); ++c) {
    const Shape& S = u.GetDisc().GetShape(*c);
    row r = u.find_row(c());
    for (int k = 0; k < c.Corners(); ++k) {
      p(c.Corner(k),0) += u(r,k);
      e(c.Corner(k),0) += 1;      
    }
  }
  p.Accumulate();
  e.Accumulate();
  for (row r=p.rows(); r!=p.rows_end(); ++r)
    p(r,0) /= e(r,0);
}
