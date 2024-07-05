#include "DGReactionAssemble.hpp"
#include "Plotting.hpp"
#include "DGElement.hpp"


void DGReactionAssemble::Initialize(Vector &u) const {
  elementPool.Initialize(u);
  faceElementPool.Initialize(u);

  /*u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RowBndValues u_c(u, *c);
    if (!u_c.onBnd()) continue;

    const DGElement &elem = elementPool.Get(u, *c);

    for (int f = 0; f < c.Faces(); ++f) {
      if (u_c.bc(f) == -1) continue;

      const DGFaceElement &faceElem = faceElementPool.Get(u, *c, f);

      Scalar F = problem.FaceConvection(c, f, faceElem.QNormal(0), faceElem.QPoint(0));
      if (F > -1e-6) continue;
      for (int j = 0; j < u.NumberOfNodalPointsOnFace(*c, f); ++j) {
        int k = u.IdOfNodalPointOnFace(*c, f, j);
        u_c(k) = problem.Concentration(Time(), elem.NodalPoint(k));
        u_c.D(k) = true;
      }
    }
  }
  u.DirichletConsistent();*/
}

double DGReactionAssemble::Residual(const Vector &u, Vector &r) const {
  r = 0;
  const Vector &u_old = PreviousSolution();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RowBndValues r_c(r, *c);

    const DGElement &elem = elementPool.Get(u, *c);

    double kappa = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point z = elem.QPoint(q);
      Tensor K = problem.Diffusion(z);
      kappa = max(kappa, norm(K));
      Scalar f = problem.Load(z);
      Scalar U = elem.Value(q, u);
      VectorField DU = elem.Derivative(q, u);
      VectorField K_DU = K * DU;
      Scalar Uold = elem.Value(q, u_old);
      Scalar R = problem.Reaction(z, U);
      VectorField B = problem.CellConvection(*c, z);
      for (int i = 0; i < elem.NodalPoints(); ++i) {
        Scalar U_i = elem.Value(q, i);
        VectorField DU_i = elem.Derivative(q, i);
        r(c(), i) += w * (K_DU * DU_i + (B * DU) * U_i
          + (1.0 / StepSize()) * (U - Uold) * U_i - (R + f) * U_i);
      }
    }
    for (int f = 0; f < c.Faces(); ++f) {
      Scalar s = penalty * c.FaceArea(f);

      const DGFaceElement &faceElem = faceElementPool.Get(u, *c, f);

      if (r.OnBoundary(*c, f)) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Point z = faceElem.QPoint(q);
          Point N = faceElem.QNormal(q);
          Scalar BN = problem.FaceConvection(*c, f, N, z);
          Tensor K = problem.Diffusion(z);
          Scalar U_diff = faceElem.Value(q, u) - problem.Concentration(Time(), z);
          VectorField K_DU = K * faceElem.Derivative(q, u);
          if (BN < 0) {
            Scalar F = K_DU * N;
            for (int i = 0; i < faceElem.NodalPoints(); ++i) {
              Scalar phi_i = faceElem.Value(q, i);
              Scalar NDphi_i = (K * faceElem.Derivative(q, i)) * N;
              r(c(), i) -= w * (F * phi_i + BN * U_diff * phi_i
                                + U_diff * (NDphi_i - s * phi_i));
            }
          } else {
            Scalar F = problem.Flux(z) * N;
            for (int i = 0; i < faceElem.NodalPoints(); ++i) {
              Scalar phi_i = faceElem.Value(q, i);
              r(c(), i) -= w * F * phi_i;
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
          Point z = faceElem.QPoint(q);
          Point N = faceElem.QNormal(q);
          Scalar BN = problem.FaceConvection(*c, f, N, z);
          Tensor K = problem.Diffusion(z);
          Scalar U = faceElem.Value(q, u);
          VectorField K_DU = K * faceElem.Derivative(q, u);
          Point Qf_c = faceElem.QPoint(q);
          int q1 = faceElem.findQPointID(otherFaceElem, Qf_c);
          Tensor K_1 = problem.Diffusion(z);
          Scalar U_1 = otherFaceElem.Value(q1, u);
          VectorField K_DU_1 = K_1 * otherFaceElem.Derivative(q1, u);
          for (int i = 0; i < faceElem.NodalPoints(); ++i) {
            Scalar phi_i = faceElem.Value(q, i);
            Scalar NDphi_i = (K * faceElem.Derivative(q, i)) * N;
            r(c(), i) += w * ((s * U - 0.5 * (K_DU * N)) * phi_i - 0.5 * U * NDphi_i);
            r(c(), i) += w * ((-s * U_1 - 0.5 * (K_DU_1 * N)) * phi_i
                              + 0.5 * U_1 * NDphi_i);
            if (BN < 0) {
              r(c(), i) -= w * BN * U * phi_i;
              r(c(), i) += w * BN * U_1 * phi_i;
            }
          }
          for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
            Scalar phi_j = otherFaceElem.Value(q1, j);
            Scalar NDphi_j = (K_1 * otherFaceElem.Derivative(q1, j)) * N;
            r(cf(), j) += w * ((0.5 * K_DU_1 * N + s * U_1) * phi_j
                               + 0.5 * U_1 * NDphi_j);
            r(cf(), j) += w * ((0.5 * K_DU * N - s * U) * phi_j
                               - 0.5 * U * NDphi_j);
            if (BN > 0) {
              r(cf(), j) += w * BN * U_1 * phi_j;
              r(cf(), j) -= w * BN * U * phi_j;
            }
          }
        }
      }
    }
  }
  r.Collect();
  return r.norm();
}

void DGReactionAssemble::Jacobi(const Vector &u, Matrix &A) const {
  A = 0;
  for (cell c = A.cells(); c != A.cells_end(); ++c) {
    int sd = c.Subdomain();

    const DGElement &elem = elementPool.Get(u, *c);

    DGRowEntries A_c(A, *c, *c);
    double kappa = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point z = elem.QPoint(q);
      Tensor K = problem.Diffusion(z);
      kappa = max(kappa, norm(K));
      Scalar U = elem.Value(q, u);
      Scalar DR = problem.DerivativeReaction(z, U);
      VectorField B = problem.CellConvection(*c, z);
      for (int i = 0; i < elem.NodalPoints(); ++i) {
        Scalar U_i = elem.Value(q, i);
        VectorField DU_i = elem.Derivative(q, i);
        VectorField K_DU_i = K * DU_i;
        for (int j = 0; j < elem.NodalPoints(); ++j) {
          Scalar U_j = elem.Value(q, j);
          VectorField DU_j = elem.Derivative(q, j);
          A_c(i, j) += w * (K_DU_i * DU_j
            + (B * DU_j) * U_i
            + (1.0 / StepSize()) * U_i * U_j
            - DR * U_i * U_j);
        }
      }
    }
    BFParts bnd(u.GetMesh(), *c);
    for (int f = 0; f < c.Faces(); ++f) {
      Scalar s = penalty * c.FaceArea(f);

      const DGFaceElement &faceElem = faceElementPool.Get(u, *c, f);

      if (u.OnBoundary(*c, f)) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Point z = faceElem.QPoint(q);
          Point N = faceElem.QNormal(q);
          Scalar BN = problem.FaceConvection(*c, f, N, z);
          if (BN < 0) {
            Tensor K = problem.Diffusion(z);
            for (int i = 0; i < faceElem.NodalPoints(); ++i) {
              Scalar phi_i = faceElem.Value(q, i);
              Scalar NDphi_i = (K * faceElem.Derivative(q, i)) * N;
              for (int j = 0; j < faceElem.NodalPoints(); ++j) {
                Scalar phi_j = faceElem.Value(q, j);
                Scalar NDphi_j = (K * faceElem.Derivative(q, j)) * N;
                A_c(i, j) -= w * (NDphi_i * phi_j + BN * phi_j * phi_i
                                  + phi_i * (NDphi_j - s * phi_j));
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
          Point z = faceElem.QPoint(q);
          Point N = faceElem.QNormal(q);
          Scalar BN = problem.FaceConvection(*c, f, N, z);
          Tensor K = problem.Diffusion(z);
          Point Qf_c = faceElem.QPoint(q);
          int q1 = faceElem.findQPointID(otherFaceElem, Qf_c);
          Tensor K_1 = problem.Diffusion(z);
          for (int i = 0; i < faceElem.NodalPoints(); ++i) {
            Scalar phi_i = faceElem.Value(q, i);
            Scalar NDphi_i = (K * faceElem.Derivative(q, i)) * N;
            for (int j = 0; j < faceElem.NodalPoints(); ++j) {
              Scalar phi_j = faceElem.Value(q, j);
              Scalar NDphi_j =
                  (K * faceElem.Derivative(q, j)) * N;
              A_c(i, j) += w * (-0.5 * NDphi_i * phi_j
                                - 0.5 * phi_i * NDphi_j
                                + s * phi_i * phi_j);
              if (BN < 0)
                A_c(i, j) -= w * BN * phi_j * phi_i;
            }
            for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
              Scalar phi_j = otherFaceElem.Value(q1, j);
              Scalar NDphi_j = (K_1 * otherFaceElem.Derivative(q1, j)) * N;
              A_cf(i, j) += w * (0.5 * NDphi_i * phi_j
                                 - phi_i * 0.5 * NDphi_j - s * phi_i * phi_j);
              A_fc(j, i) += w * (0.5 * NDphi_i * phi_j
                                 - phi_i * 0.5 * NDphi_j - s * phi_i * phi_j);
              if (BN < 0)
                A_cf(i, j) += w * BN * phi_j * phi_i;
              else
                A_fc(j, i) -= w * BN * phi_j * phi_i;
            }
          }
          for (int i = 0; i < otherFaceElem.NodalPoints(); ++i) {
            Scalar phi_i = otherFaceElem.Value(q1, i);
            Scalar NDphi_i =
                (K_1 * otherFaceElem.Derivative(q1, i)) * N;
            for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
              Scalar phi_j = otherFaceElem.Value(q1, j);
              Scalar NDphi_j =
                  (K_1 * otherFaceElem.Derivative(q1, j)) * N;

              A_ff(i, j) += w * (0.5 * NDphi_i * phi_j
                                 + phi_i * 0.5 * NDphi_j
                                 + s * phi_i * phi_j);

              if (BN > 0) A_ff(i, j) += w * BN * phi_j * phi_i;
            }
          }
        }
      }
    }
  }
}

RatePair DGReactionAssemble::InflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    int sd = c.Subdomain();
    for (int f = 0; f < c.Faces(); ++f) {
      if (bnd[f] == -1) continue;

      const DGFaceElement &faceElem = faceElementPool.Get(u, *c, f);

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Scalar F = problem.FaceConvection(*c, f, faceElem.QNormal(q),
                                          faceElem.QPoint(q));
        Scalar U = faceElem.Value(q, u);
        if (F > 0) outflow += w * F * U;
        else inflow += w * F * U;
      }
    }
  }
  return {PPM->SumOnCommSplit(inflow, u.CommSplit()),
          PPM->SumOnCommSplit(outflow, u.CommSplit())};
}

void DGReactionAssemble::FinishTimeStep(const Vector &u) {
  PlotIteration(u);
  PrintIteration(u);
}

void DGReactionAssemble::SetInitialValue(Vector &u) {
  u = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    row r = u.find_row(c());

    const DGElement &elem = elementPool.Get(u, *c);

    double lambda = 1e-9;
    for (int j = 0; j < elem.NodalPoints(); ++j)
      u(r, j) = problem.Concentration(0, lambda * c() + (1 - lambda) * elem.NodalPoint(j));
  }
  u.Accumulate();
}

double DGReactionAssemble::Mass(const Vector &u) const {
  Scalar e = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const DGElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar U = elem.Value(q, u);
      e += w * U;
    }
  }
  return PPM->SumOnCommSplit(e, u.CommSplit());
}

void DGReactionAssemble::PlotIteration(const Vector &u) const {
  // Todo SampleSolution tmp(*cellDisc, GetProblem().CurrentID(), "C");
  if (plotting == 0 || Step() % plotting != 0) return;
  auto cellDisc = std::make_shared<const DGDiscretization>(u.GetDisc().GetMeshes(), 0, 1);
  Vector tmp(0.0, cellDisc);
  for (cell c = tmp.cells(); c != tmp.cells_end(); ++c) {
    DGElement elem(tmp, *c);
    Scalar U = 0.0;
    double a = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      U += w * elem.Value(q, u);
      a += w;
    }
    U *= (1 / a);
    row r = u.find_row(c());
    tmp(r)[0] = U;
  }
  std::string filename = "C." + to_string(Step());
  mpp::plot(filename) << tmp << mpp::endp;
}
