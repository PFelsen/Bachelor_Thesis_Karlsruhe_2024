#include "PGReactionAssemble.hpp"
#include "ScalarElement.hpp"


void PGReactionAssemble::Initialize(Vector &u) const {
  elementPool.Initialize(u);
  faceElementPool.Initialize(u);

  u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RowBndValues u_c(u, *c);
    if (!u_c.onBnd()) continue;

    const ScalarElement &elem = elementPool.Get(u, *c);

    for (int i = 0; i < c.Faces(); ++i) {
      if (u_c.bc(i) == -1) continue;
      ScalarFaceElement FE(u, *c, i);
      Scalar F = problem.FaceConvection(*c, i, FE.QNormal(0), FE.QPoint(0));
      if (F > -1e-6) continue;
      for (int j = 0; j < u.NumberOfNodalPointsOnFace(*c, i); ++j) {
        int k = u.IdOfNodalPointOnFace(*c, i, j);
        u_c(k) = problem.Concentration(Time(), elem.NodalPoint(k));
        u_c.D(k) = true;
      }
    }
  }
  u.DirichletConsistent();
}

double PGReactionAssemble::Residual(const Vector &u, Vector &r) const {
  r = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const ScalarElement &elem = elementPool.Get(u, *c);

    RowBndValues r_c(r, *c);
    double delta_h = Delta(c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar U = elem.Value(q, u);
      Scalar Uold = elem.Value(q, PreviousSolution());
      Tensor K = problem.Diffusion(elem.QPoint(q));
      Scalar f = problem.Load(elem.QPoint(q));
      VectorField C = problem.CellConvection(*c, elem.QPoint(q));
      Scalar R = problem.Reaction(elem.QPoint(q), U);
      VectorField DU = elem.Derivative(q, u);
      for (unsigned int i = 0; i < elem.size(); ++i) {
        Scalar U_i = elem.Value(q, i);
        VectorField DU_i = elem.Derivative(q, i);
        r_c(i) += w * ((K * DU) * DU_i + U_i * (C * DU)
                       + (1.0 / StepSize()) * (U - Uold) * U_i - (R + f) * U_i
                       + delta_h * ((C * DU) + (1.0 / StepSize()) * (U - Uold) - (R + f))
                         * (C * DU_i));
      }
      BFParts bnd(u.GetMesh(), *c);
      if (bnd.onBnd()) continue;
      for (int f = 0; f < c.Faces(); ++f) {
        if (bnd[f] == -1) continue;

        const ScalarFaceElement &faceElem = faceElementPool.Get(u, *c, f);

        RowValues r_f(r, faceElem);
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Scalar F = problem.FaceConvection(*c, f, faceElem.QNormal(q), faceElem.QPoint(q));
          if (F < 1e-6) continue;
          Scalar U = 0;
          for (unsigned int i = 0; i < faceElem.size(); ++i) {
            Scalar U_i = faceElem.Value(q, i);
            U += u(faceElem[i], 0) * U_i;
          }
          for (unsigned int j = 0; j < faceElem.size(); ++j) {
            Scalar U_j = faceElem.Value(q, j);
            r_f(j) += w * F * U * U_j;
          }
        }
      }
    }
  }
  r.ClearDirichletValues();
  r.Collect();
  return r.norm();
}

void PGReactionAssemble::Jacobi(const Vector &u, Matrix &A) const {
  A = 0;
  for (cell c = A.cells(); c != A.cells_end(); ++c) {

    const ScalarElement &elem = elementPool.Get(u, *c);

    RowEntries A_c(A, elem);
    double delta_h = Delta(c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar U = elem.Value(q, u);
      Tensor K = problem.Diffusion(elem.QPoint(q));
      VectorField C = problem.CellConvection(*c, elem.QPoint(q));
      Scalar DR = problem.DerivativeReaction(elem.QPoint(q), U);
      for (unsigned int i = 0; i < elem.size(); ++i) {
        VectorField DU_i = elem.Derivative(q, i);
        Scalar U_i = elem.Value(q, i);
        for (unsigned int j = 0; j < elem.size(); ++j) {
          VectorField DU_j = elem.Derivative(q, j);
          Scalar U_j = elem.Value(q, j);
          A_c(i, j) += w * ((K * DU_i) * DU_j + (C * DU_j) * U_i
                            + (1.0 / StepSize()) * U_i * U_j - DR * U_i * U_j
                            + delta_h * ((C * DU_j) + (1.0 / StepSize()) * U_j - DR * U_j) *
                              (C * DU_i));
        }
      }
      BFParts bnd(u.GetMesh(), *c);
      if (bnd.onBnd()) continue;
      for (int f = 0; f < c.Faces(); ++f) {
        if (bnd[f] == -1) continue;

        const ScalarFaceElement &faceElem = faceElementPool.Get(u, *c, f);

        RowEntries A_f(A, faceElem);
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Scalar F = problem.FaceConvection(*c, f, faceElem.QNormal(q), faceElem.QPoint(q));
          if (F < 1e-6) continue;
          for (unsigned int i = 0; i < faceElem.size(); ++i) {
            Scalar U_i = faceElem.Value(q, i);
            for (unsigned int j = 0; j < faceElem.size(); ++j) {
              Scalar U_j = faceElem.Value(q, j);
              A_f(i, j) += w * F * U_i * U_j;
            }
          }
        }
      }
    }
  }
  A.ClearDirichletValues();
}

RatePair PGReactionAssemble::InflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    int sd = c.Subdomain();
    for (int f = 0; f < c.Faces(); ++f) {
      if (bnd[f] == -1) continue;

      const ScalarFaceElement &faceElem = faceElementPool.Get(u, *c, f);

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
  inflow = PPM->SumOnCommSplit(inflow, u.GetMesh().CommSplit());
  outflow = PPM->SumOnCommSplit(outflow, u.GetMesh().CommSplit());
  if (inflow < 1e-7)
    inflow = 0.0;
  if (outflow < 1e-7)
    outflow = 0.0;
  return {inflow, outflow};
}

double PGReactionAssemble::Mass(const Vector &u) const {
  Scalar e = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const ScalarElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar U = elem.Value(q, u);
      e += w * U;
    }
  }
  return PPM->SumOnCommSplit(e, u.GetMesh().CommSplit());
}

void PGReactionAssemble::SetInitialValue(Vector &u) {
  u = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const ScalarElement &elem = elementPool.Get(u, *c);

    RowBndValues u_c(u, *c);
    for (unsigned int i = 0; i < elem.size(); ++i)
      u_c(i) = problem.Concentration(0, elem[i]());
  }
}

void PGReactionAssemble::FinishTimeStep(const Vector &u) {
  PlotIteration(u);
  PrintIteration(u);
}

void PGReactionAssemble::PrintIteration(const Vector &u) const {
  if (verbose == 0 || Step() % verbose != 0) return;
  std::pair<double, double> rate = InflowOutflow(u);
  mout.PrintIteration(
      verbose,
      PrintIterEntry("n", Step(), 5, 1),
      PrintIterEntry("t", Time(), 13, 1),
      PrintIterEntry("Mass", Mass(u), 13, 1),
      PrintIterEntry("IFR", rate.first, 13, 1),
      PrintIterEntry("OFR", rate.second, 13, 1)
  );
}

void PGReactionAssemble::PlotIteration(const Vector &u) const {
  // Todo SampleSolution tmp(u, GetProblem().CurrentID(), "C");
  if (plotting == 0 || Step() % plotting != 0) return;
  std::string filename = "C." + to_string(Step());
  mpp::plot(filename) << u << mpp::endp;
}

void PGReactionAssemble::PrintInfo() const {
  int verbose = 1;
  mout.PrintInfo("Assemble", verbose,
                 PrintInfoEntry("Name", Name()),
//                 PrintInfoEntry("Problem", problem.Name()),
                 PrintInfoEntry("Convection", problem.GetConvection()),
                 PrintInfoEntry("Diffusion", problem.GetDiffusion()),
                 PrintInfoEntry("Reaction", problem.GetReaction()),
                 PrintInfoEntry("Delta", delta),
                 PrintInfoEntry("Discretization", disc->DiscName()));
}

double PGReactionAssemble::Delta(const cell &c) const {
  double h = diam(c);
  Tensor K = problem.Diffusion(c());
  VectorField C = problem.CellConvection(*c, c());
  double Pe = 0.5 * h * norm(C) / norm(K);
  if (Pe > 1) return h * delta;
  return h * h * delta / norm(K);
}

double PGReactionAssemble::diam(const cell &c) const {
  double h = 0;
  for (int i = 0; i < c.Corners(); ++i)
    for (int j = i + 1; j < c.Corners(); ++j)
      h = max(h, dist(c[i], c[j]));
  return h;
}
