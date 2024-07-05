#include "EGEllipticAssemble.hpp"
#include "MixedEGElement.hpp"

const char *EGEllipticAssemble::Name() const {
  return "EGEllipticAssemble";
}

void EGEllipticAssemble::Initialize(Vector &u) const {

  elementPool.Initialize(u);
  faceElementPool.Initialize(u);

  u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RowBndValues u_c(u, *c);
    if (!u_c.onBnd()) continue;

    const MixedEGElement &elem = elementPool.Get(u, *c);

    for (int face = 0; face < c.Faces(); ++face) {
      if (u_c.bc(face) == 1) {
        for (int j = 0; j < u.NumberOfNodalPointsOnFace(*c, face); ++j) {
          int k = u.IdOfNodalPointOnFace(*c, face, j);
          u_c(k, 0) = problem.Solution(elem.NodalPoint(k));
          u_c.D(k, 0) = true;
        }
      }
    }
  }
  u.DirichletConsistent();

/*
  elementPool.Initialize(u);
  faceElementPool.Initialize(u);

  u.ClearDirichletFlags();
  ExchangeBuffer EB(u.CommSplit());
  Point z = FindDirichlet(u);
  if (z != infty)
    for (int q = 0; q < PPM->Size(u.CommSplit()); ++q)
      EB.Send(q) << z;
  EB.Communicate();
  z = infty;
  for (int q = 0; q < PPM->Size(u.CommSplit()); ++q) {
    while (EB.Receive(q).size() < EB.Receive(q).Size()) {
      EB.Receive(q) >> z;
      break;
    }
    if (z != infty) break;
  }
  row r = u.find_row(z);
  if (r != u.rows_end()) {
    u(r, 0) = problem.Solution(z);
    u.D(r, 0) = true;
  }
  u.DirichletConsistent();*/
}

void EGEllipticAssemble::Residual(const cell &c, const Vector &u, Vector &r) const {

  const MixedEGElement &elem = elementPool.Get(u, *c);

  MixedRowBndValues r_c(r, *c);
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    const Point &z = elem.QPoint(q);
    Tensor K = problem.Permeability(z);
    Scalar f = problem.Load(z, c);
    VectorField K_DU = K * elem.Derivative(q, u);
    for (int i = 0; i < elem.ValueSize(); ++i) {
      Scalar U_i = elem.Value(q, i);
      VectorField DU_i = elem.Derivative(q, i);
      r_c(elem.ValueIndex(), i) += w * (K_DU * DU_i - f * U_i);
    }
    r_c(elem.PenaltyIndex(), 0) +=
        w * (K_DU * elem.PenaltyGradient(q) - f * elem.PenaltyValue(q));
  }

  for (int f = 0; f < c.Faces(); ++f) {
    Scalar s = (pow(degree,2)*penalty)/  u.GetMesh().MaxMeshWidth();

    const MixedEGFaceElement &faceElem = faceElementPool.Get(u, *c, f);

    if (r.OnBoundary(*c, f)) {
      int bnd = r_c.bc(f);
      if (bnd == 2) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          const Point &z = faceElem.QPoint(q);
          Point N = faceElem.QNormal(q);
          Scalar F = problem.Flux(z) * N;
          for (int i = 0; i < faceElem.ValueSize(); ++i) {
            Scalar phi_i = faceElem.Value(q, i);
            r_c(faceElem.ValueIndex(), i) -= w * phi_i * F;
          }
          r_c(faceElem.PenaltyIndex(), 0) -= w * faceElem.PenaltyValue(q) * F;
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


          Scalar theta = faceElem.PenaltyValue(q);
          Scalar Dtheta = (K * faceElem.PenaltyGradient(q)) * N;
          for (int i = 0; i < faceElem.ValueSize(); ++i) {
            Scalar phi_i = faceElem.Value(q, i);
            Scalar NDphi_i = (K * faceElem.Derivative(q, i)) * N;
            r_c(faceElem.ValueIndex(), i)
                -= w * (F * phi_i
                        - U_diff *
                          (sign * NDphi_i
                           + s * (K * N) * N * phi_i));
          }
          r_c(faceElem.PenaltyIndex(), 0)
              -= w * (F * theta
                      - U_diff * (sign * Dtheta
                                  + s * (K * N) * N * theta));
        }
      }
    } else {
      cell cf = r.find_neighbour_cell(c, f);
      if (cf() < c()) continue;
      int f1 = r.find_neighbour_face_id(c.Face(f), cf);

      const MixedEGFaceElement &otherFaceElem = faceElementPool.Get(u, *cf, f1);

      MixedRowBndValues r_cf(r, *cf);
      Tensor K = problem.Permeability(c()); // kappa_+
      Tensor K_1 = problem.Permeability(cf()); //kappa_-
      for (int q = 0; q < faceElem.nQ(); ++q) {
        const Point &z = faceElem.QPoint(q);
        int q1 = faceElem.FindQPointID(otherFaceElem, z);
        double w = faceElem.QWeight(q);
        const Point &N = faceElem.QNormal(q);
        double K_e = 2 * ((K * N) * N) * ((K_1 * N) * N) / ((K * N) * N + (K_1 * N) * N); // kappa_e

        Scalar K_DU = (faceElem.Derivative(q, u)) * N; // kappa grad v * N
        Scalar U = faceElem.Value(q, u); // v
        Scalar K_DU_1 = (otherFaceElem.Derivative(q1, u)) * N; // kappa grad u_nghbr
        Scalar U_1 = otherFaceElem.Value(q1, u); // u_nghbr

        for (int i = 0; i < faceElem.ValueSize(); ++i) {
          Scalar phi_i = faceElem.Value(q, i); // w
          Scalar NDphi_i = faceElem.Derivative(q, i) * N; // grad w
          r_c(faceElem.ValueIndex(), i) +=
              w * K_e * (s * U * phi_i - 0.5 * K_DU * phi_i + 0.5 * sign * U * NDphi_i);
          r_c(faceElem.ValueIndex(), i) +=
              w * K_e * (-s * U_1 * phi_i - 0.5 * K_DU_1 * phi_i - 0.5 * sign * U_1 * NDphi_i);
        }
        for (int j = 0; j < otherFaceElem.ValueSize(); ++j) {
          Scalar phi_j_1 = otherFaceElem.Value(q1, j); // w_nghbr
          Scalar NDphi_j_1 = otherFaceElem.Derivative(q1, j) * N; // grad w_nghbr
          r_cf(otherFaceElem.ValueIndex(), j) +=
              w * K_e * (s * U_1 * phi_j_1 + 0.5 * K_DU_1 * phi_j_1 - 0.5 * sign * U_1 * NDphi_j_1);
          r_cf(otherFaceElem.ValueIndex(), j) +=
              w * K_e * (-s * U * phi_j_1 + 0.5 * K_DU * phi_j_1 + 0.5 * sign * U * NDphi_j_1);
        }

        Scalar theta_i = faceElem.PenaltyValue(q);
        r_c(faceElem.PenaltyIndex(), 0) +=
            w * K_e * (s * U * theta_i - 0.5 * K_DU * theta_i /*+ 0.5 * sign * U * NDphi_i*/);
        r_c(faceElem.PenaltyIndex(), 0) += w * K_e * (-s * U_1 * theta_i - 0.5 * K_DU_1 *
                                                                           theta_i /*- 0.5 * sign * U_1 * NDphi_i*/);

        Scalar theta_j = otherFaceElem.PenaltyValue(q1);
        r_cf(otherFaceElem.PenaltyIndex(), 0) += w * K_e * (s * U_1 * theta_j + 0.5 * K_DU_1 *
                                                                                theta_j /*- 0.5 * sign * U_1 * NDphi_j_1*/);
        r_cf(otherFaceElem.PenaltyIndex(), 0) +=
            w * K_e * (-s * U * theta_j + 0.5 * K_DU * theta_j /*+ 0.5 * sign * U * NDphi_j_1*/);

      }
    }
  }
}

void EGEllipticAssemble::Jacobi(const cell &c, const Vector &u, Matrix &A) const {

  const MixedEGElement &elem = elementPool.Get(u, *c);

  MixedRowEntries A_c(A, *c);
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    Point z = elem.QPoint(q);
    Tensor K = problem.Permeability(z);
    for (int i = 0; i < elem.ValueSize(); ++i) {
      VectorField K_DU_i = K * elem.Derivative(q, i);
      for (int j = 0; j < elem.ValueSize(); ++j) {
        VectorField DU_j = elem.Derivative(q, j);
        A_c(elem.ValueIndex(), elem.ValueIndex(), i, j) += w * K_DU_i * DU_j;
      }
    }
  }
  BFParts bnd(u.GetMesh(), *c);
  for (int f = 0; f < c.Faces(); ++f) {
    Scalar s = (pow(degree,2)*penalty)/  u.GetMesh().MaxMeshWidth();

    const MixedEGFaceElement &faceElem = faceElementPool.Get(u, *c, f);

    if (elem.Bnd(f) == 1) {
      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        const Point &z = faceElem.QPoint(q);
        const Point &N = faceElem.QNormal(q);
        Tensor K = problem.Permeability(c());

        Scalar theta = faceElem.PenaltyValue(q);
        Scalar Dtheta = (K * faceElem.PenaltyGradient(q)) * N; // == zero
        for (int i = 0; i < faceElem.ValueSize(); ++i) {
          Scalar phi_i = faceElem.Value(q, i);
          Scalar NDphi_i = (K * faceElem.Derivative(q, i)) * N;
          for (int j = 0; j < faceElem.ValueSize(); ++j) {
            Scalar phi_j = faceElem.Value(q, j);
            Scalar NDphi_j = (K * faceElem.Derivative(q, j)) * N;
            A_c(faceElem.ValueIndex(), faceElem.ValueIndex(), i, j)
                -= w * (phi_i * NDphi_j
                        - phi_j * (
                sign * NDphi_i
                + s * (K * N) * N * phi_i));
          }
          A_c(faceElem.ValueIndex(), faceElem.PenaltyIndex(), i, 0)
              -= w * (Dtheta * phi_i
                      - theta * (sign * NDphi_i + s * (K * N) * N * phi_i));

          A_c(faceElem.PenaltyIndex(), faceElem.ValueIndex(), 0, i)
              -= w * (NDphi_i * theta
                      - phi_i * (sign * Dtheta +
                                 s * (K * N) * N * theta));

        }

        A_c(faceElem.PenaltyIndex(), faceElem.PenaltyIndex(), 0, 0)
            -= w * (/*Dtheta * theta*/   -theta *
                                         (/*sign * Dtheta +*/ s * (K * N) *
                                                              N * theta));
      }
    } else if (elem.Bnd(f) == -1) {
      cell cf = u.find_neighbour_cell(c, f);
      if (cf() < c()) continue;
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);

      const MixedEGFaceElement &otherFaceElem = faceElementPool.Get(u, *cf, f1);

      MixedRowEntries A_cf(A, *c, *cf);
      MixedRowEntries A_fc(A, *cf, *c);
      MixedRowEntries A_f(A, *cf);

      Tensor K = problem.Permeability(c());
      Tensor K_1 = problem.Permeability(cf());
      for (int q = 0; q < faceElem.nQ(); ++q) {
        const Point &z = faceElem.QPoint(q);
        int q1 = otherFaceElem.FindQPointID(faceElem, z);
        double w = faceElem.QWeight(q);
        const Point &N = faceElem.QNormal(q);
        double K_e = 2 * ((K * N) * N) * ((K_1 * N) * N) / ((K * N) * N + (K_1 * N) * N);

        Scalar theta = faceElem.PenaltyValue(q);
        Scalar Dtheta = (faceElem.PenaltyGradient(q)) * N; // == zero


        Scalar theta_1 = otherFaceElem.PenaltyValue(q1);
        Scalar Dtheta_1 = (otherFaceElem.PenaltyGradient(q1)) * N; // == zero

        for (int i = 0; i < faceElem.ValueSize(); ++i) {
          Scalar phi_i = faceElem.Value(q, i);
          Scalar NDphi_i = faceElem.Derivative(q, i) * N;
          for (int j = 0; j < faceElem.ValueSize(); ++j) {
            Scalar phi_j = faceElem.Value(q, j);
            Scalar NDphi_j = faceElem.Derivative(q, j) * N;
            A_c(faceElem.ValueIndex(), faceElem.ValueIndex(), i, j) +=
                w * K_e * (s * phi_i * phi_j - 0.5 * phi_i * NDphi_j
                           + 0.5 * sign * NDphi_i * phi_j);
          }

          A_c(faceElem.ValueIndex(), faceElem.PenaltyIndex(), i, 0)
              += w * K_e * (s * phi_i * theta /*- 0.5 * phi_i * Dtheta*/
                            + 0.5 * sign * NDphi_i * theta);

          A_c(faceElem.PenaltyIndex(), faceElem.ValueIndex(), 0, i)
              += w * K_e * (s * theta * phi_i - 0.5 * theta * NDphi_i
              /*+ 0.5 * sign * Dtheta * phi_i*/);

          for (int j = 0; j < otherFaceElem.ValueSize(); ++j) {
            Scalar phi_j = otherFaceElem.Value(q1, j);
            Scalar NDphi_j = otherFaceElem.Derivative(q1, j) * N;
            A_cf(faceElem.ValueIndex(), otherFaceElem.ValueIndex(), i, j)
                += w * K_e * (-s * phi_i * phi_j
                              - phi_i * 0.5 * NDphi_j
                              - 0.5 * NDphi_i * phi_j * sign);


            A_fc(otherFaceElem.ValueIndex(), faceElem.ValueIndex(), j, i) +=
                w * K_e * (-s * phi_i * phi_j + 0.5 * NDphi_i * phi_j
                           + phi_i * 0.5 * NDphi_j * sign);
          }
          A_cf(faceElem.ValueIndex(), otherFaceElem.PenaltyIndex(), i, 0)
              += w * K_e * (-s * phi_i * theta_1
                            - phi_i * 0.5 * Dtheta_1
                            - 0.5 * NDphi_i * theta_1 * sign);


          A_fc(otherFaceElem.PenaltyIndex(), faceElem.ValueIndex(), 0, i) +=
              w * K_e * (-s * phi_i * theta_1 + 0.5 * NDphi_i * theta_1
                  /*+ phi_i * 0.5 * Dtheta_1 * sign*/);

        }
        for (int i = 0; i < otherFaceElem.ValueSize(); ++i) {
          Scalar phi_i = otherFaceElem.Value(q1, i);
          Scalar NDphi_i = otherFaceElem.Derivative(q1, i) * N;


          A_cf(faceElem.PenaltyIndex(), otherFaceElem.ValueIndex(), 0, i)
              += w * K_e * (-s * theta * phi_i
                            - theta * 0.5 * NDphi_i
              /*- 0.5 * Dtheta * phi_i * sign*/);

          A_fc(otherFaceElem.ValueIndex(), faceElem.PenaltyIndex(), i, 0) +=
              w * K_e * (-s * theta * phi_i + 0.5 * Dtheta * phi_i
                         + theta * 0.5 * NDphi_i * sign);

          for (int j = 0; j < otherFaceElem.ValueSize(); ++j) {
            Scalar phi_j = otherFaceElem.Value(q1, j);
            Scalar NDphi_j = otherFaceElem.Derivative(q1, j) * N;
            A_f(otherFaceElem.ValueIndex(), otherFaceElem.ValueIndex(), i, j) +=
                w * K_e * (s * phi_i * phi_j + phi_i * 0.5 * NDphi_j
                           - 0.5 * NDphi_i * phi_j * sign);

          }
          A_f(otherFaceElem.ValueIndex(), otherFaceElem.PenaltyIndex(), i, 0) +=
              w * K_e * (s * phi_i * theta_1 + phi_i * 0.5 * Dtheta_1
                         - 0.5 * NDphi_i * theta_1 * sign);

          A_f(otherFaceElem.PenaltyIndex(), otherFaceElem.ValueIndex(), 0, i) +=
              w * K_e * (s * theta_1 * phi_i + theta_1 * 0.5 * NDphi_i
                  /*- 0.5 * Dtheta_1 * phi_i * sign*/);

        }


        A_c(faceElem.PenaltyIndex(), faceElem.PenaltyIndex(), 0, 0)
            += w * K_e * (s * theta * theta /*- 0.5 * theta * Dtheta
                          + 0.5 * sign * Dtheta * theta*/);

        A_cf(faceElem.PenaltyIndex(), otherFaceElem.PenaltyIndex(), 0, 0)
            += w * K_e * (-s * theta * theta_1
            /*- theta * 0.5 * Dtheta_1
            - 0.5 * Dtheta * theta_1 * sign*/);

        A_fc(otherFaceElem.PenaltyIndex(), faceElem.PenaltyIndex(), 0, 0) +=
            w * K_e * (-s * theta * theta_1
                /*+ 0.5 * Dtheta * theta_1
                           + theta * 0.5 * Dtheta_1 * sign*/);

        A_f(otherFaceElem.PenaltyIndex(), otherFaceElem.PenaltyIndex(), 0, 0) +=
            w * K_e * (s * theta_1 * theta_1
                /*+ theta_1 * 0.5 * Dtheta_1
                           - 0.5 * Dtheta_1 * theta_1 * sign*/);
      }
    }
  }
}

FluxPair EGEllipticAssemble::InflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] < 0) continue;

      const MixedEGFaceElement &faceElem = faceElementPool.Get(u, *c, face);

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

FluxPair EGEllipticAssemble::PrescribedInflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] < 2) continue;

      const MixedEGFaceElement &faceElem = faceElementPool.Get(u, *c, face);

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

FluxPair EGEllipticAssemble::OutflowLeftRight(const Vector &u) const {
  double outflowLeft = 0;
  double outflowRight = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] != 1) continue;

      const MixedEGFaceElement &faceElem = faceElementPool.Get(u, *c, face);

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

double EGEllipticAssemble::GoalFunctional(const Vector &u) const {
  double goal = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] != 1) continue;

      const MixedEGFaceElement &faceElem = faceElementPool.Get(u, *c, face);

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

void EGEllipticAssemble::SetExactSolution(Vector &uEx) const {
  uEx = 0;
  for (cell c = uEx.cells(); c != uEx.cells_end(); ++c) {

    const MixedEGElement &elem = elementPool.Get(uEx, *c);

    for (int i = 0; i < elem.NodalPoints() - 1; ++i)
      uEx(elem.NodalPoint(i), 0) = problem.Solution(elem.NodalPoint(i));
    uEx(c(), 0) = problem.Solution(c()) - elem.Value(c.LocalCenter(), uEx);
  }
  uEx.MakeAdditive();
  uEx.Collect();

}

void EGEllipticAssemble::SetFlux(const Vector &u, Vector &flux) {
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const MixedEGElement &elem = elementPool.Get(u, *c);

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

Point EGEllipticAssemble::FindDirichlet(Vector &u) const {
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RowBndValues u_c(u, *c);
    if (!u_c.onBnd()) continue;

    const MixedEGElement &elem = elementPool.Get(u, *c);

    for (int i = 0; i < c.Faces(); ++i) {
      if (u_c.bc(i) == 1) {
        int k = u.IdOfNodalPointOnFace(*c, i, 0);
        return elem[k]();
      }
    }
  }
  return infty;
}

void EGEllipticAssemble::SetPressure(const Vector &u, Vector &p) const {
  for (row r=p.rows(); r!=p.rows_end(); ++r)
    p(r,0) = u(r(),0);
}
