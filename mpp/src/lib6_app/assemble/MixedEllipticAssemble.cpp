#include "MixedEllipticAssemble.hpp"
#include "LagrangeEllipticAssemble.hpp"
#include "Newton.hpp"
#include "RTElement.hpp"
#include "RTLagrangeElement.hpp"
#include "ScalarElement.hpp"

//TODO: Extend FaceElement parts to use higher order elements

const char *MixedEllipticAssemble::Name() const {
  return "MixedEllipticAssemble";
}

void MixedEllipticAssemble::Initialize(Vector &u) const {
  elementPool.Initialize(u);
  faceElementPool.Initialize(u);

  u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      const RTFaceElementT<> &faceElem = faceElementPool.Get(u, *c, face);
      RowBndValues uf_c(u, *c);
      if (bnd[face] == 2) {
        VectorField F = problem.Flux(c.Face(face));
        VectorField N = faceElem.Normal();
        uf_c(face) = F * N * faceElem.Area();
        uf_c.D(face) = true;
      } else if (bnd[face] == 0) {
        uf_c(face) = 0.0;
        uf_c.D(face) = true;
      }
    }
  }
  u.DirichletConsistent();
}

double MixedEllipticAssemble::Energy(const Vector &u) const {
  double energy = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Tensor IK = Invert(problem.Permeability(elem.QPoint(q)));
      VectorField Q = elem.VelocityField(q, u);
      energy += w * IK * Q * Q;
    }
  }
  return sqrt(PPM->SumOnCommSplit(energy, u.GetMesh().CommSplit()));
}

void MixedEllipticAssemble::Residual(const cell &c, const Vector &u, Vector &r) const {

  const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

  MixedRowBndValues r_c(r, *c);
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    Scalar f = problem.Load(elem.QPoint(q), c);
    VectorField Q = elem.VelocityField(q, u);
    Tensor IK = Invert(problem.Permeability(elem.QPoint(q)));
    Scalar divQ = elem.VelocityFieldDivergence(q, u);
    Scalar p = elem.PressureValue(q, u);
    for (int i = 0; i < elem.VelocitySize(); ++i) {
      for (int k = 0; k < elem.VelocityMaxk(i); ++k) {
        Velocity Q_i = elem.VelocityField(q, i, k);
        Scalar divQ_i = elem.VelocityFieldDivergence(q, i, k);
        r_c(elem.VelocityIndex(), i, k) += w * (((IK * Q) * Q_i) + divQ_i * p);
      }
    }
    for (unsigned int j = 0; j < elem.PressureSize(); ++j) {
      Scalar p_j = elem.PressureValue(q, j);
      r_c(elem.PressureIndex(), j, 0) += w * (divQ - f) * p_j;
    }
  }
  if (!r_c.onBnd()) return;
  for (int face = 0; face < c.Faces(); ++face) {
    if (r_c.bc(face) == 1) {

      const RTFaceElementT<> &faceElem = faceElementPool.Get(u, *c, face);

      RowBndValues rf_c(r, *c);
      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Scalar P = problem.Solution(faceElem.QPoint(q));
//        for (int i = 0; i < faceElem.size(); i++) {
        Velocity face_Q_i = faceElem.VelocityField(q, face);
        double uin = face_Q_i * faceElem.QNormal(q);
        rf_c(face) -= w * uin * P;
//        }
      }
    }
  }
}

void MixedEllipticAssemble::Jacobi(const cell &c, const Vector &u, Matrix &A) const {

  const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

  MixedRowEntries A_c(A, *c, elem);
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    Tensor IK = Invert(problem.Permeability(elem.QPoint(q)));
    for (int i = 0; i < elem.VelocitySize(); ++i) {
      for (int k = 0; k < elem.VelocityMaxk(i); ++k) {
        VectorField Q_i = elem.VelocityField(q, i, k);
        VectorField IKQ_i = IK * Q_i;
        Scalar divQ_i = elem.VelocityFieldDivergence(q, i, k);
        for (int j = 0; j < elem.VelocitySize(); ++j) {
          for (int l = 0; l < elem.VelocityMaxk(j); ++l) {
            Velocity Q_j = elem.VelocityField(q, j, l);
            A_c(elem.VelocityIndex(), elem.VelocityIndex(), i, j, k, l) += w * (IKQ_i * Q_j);
          }
        }
        for (unsigned int j = 0; j < elem.PressureSize(); ++j) {
          Scalar p_j = elem.PressureValue(q, j);
          A_c(elem.PressureIndex(), elem.VelocityIndex(), j, i, 0, k) += w * divQ_i * p_j;
          A_c(elem.VelocityIndex(), elem.PressureIndex(), i, j, k, 0) += w * divQ_i * p_j;
        }
      }
    }
  }
}

double MixedEllipticAssemble::EnergyError(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Tensor IK = Invert(problem.Permeability(elem.QPoint(q)));
      VectorField Q = elem.VelocityField(q, u);
      VectorField F = problem.Flux(elem.QPoint(q));
      VectorField Diff = Q - F;
      err += w * IK * Diff * Diff;
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.GetMesh().CommSplit()));
}

double MixedEllipticAssemble::L2(const Vector &u) const {
  double l2 = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar p = elem.PressureValue(q, u);
      l2 += w * p * p;
    }
  }
  return sqrt(PPM->SumOnCommSplit(l2, u.GetMesh().CommSplit()));
}

double MixedEllipticAssemble::H1(const Vector &u) const {
  double l2 = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Tensor IK = Invert(problem.Permeability(elem.QPoint(q)));
      VectorField Q = elem.VelocityField(q, u);
      Scalar p = elem.PressureValue(q, u);
      VectorField IK_Q = IK * Q;
      l2 += w * IK_Q * IK_Q + w * p * p;
    }
  }
  return sqrt(PPM->SumOnCommSplit(l2, u.GetMesh().CommSplit()));
}

double MixedEllipticAssemble::L2Error(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar p = elem.PressureValue(q, u);
      Scalar Sol = problem.Solution(elem.QPoint(q));
      err += w * (p - Sol) * (p - Sol);
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.GetMesh().CommSplit()));
}

double MixedEllipticAssemble::L2CellAvgError(const Vector &u) const {
  double err = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    double w = 0;
    Scalar p = 0;
    Scalar Sol = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      w += elem.QWeight(q);
      p += elem.PressureValue(q, u);
      Sol += problem.Solution(elem.QPoint(q));
    }
    err += w * (p - Sol) * (p - Sol);
  }
  return sqrt(PPM->SumOnCommSplit(err, u.GetMesh().CommSplit()));
}

double MixedEllipticAssemble::MaxError(const Vector &u) const {
  double err = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      Scalar p = elem.PressureValue(q, u);
      err = max(err, abs(p - problem.Solution(elem.QPoint(q))));
    }
  }
  return PPM->Max(err, u.GetMesh().CommSplit());
}

double MixedEllipticAssemble::FluxError(const Vector &u) const {
  double flux_error = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] == 2) {

        const RTFaceElementT<> &faceElem = faceElementPool.Get(u, *c, face);

        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          VectorField F = problem.Flux(faceElem.QPoint(q));
          Scalar Diff = (faceElem.VelocityField(q, u) - F) * faceElem.QNormal(q);
          flux_error += w * Diff * Diff;
        }
      } else if (bnd[face] == 0) {

        const RTFaceElementT<> &faceElem = faceElementPool.Get(u, *c, face);

        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Scalar FN = faceElem.VelocityField(q, u) * faceElem.QNormal(q);
          flux_error += w * FN * FN;
        }
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(flux_error, u.GetMesh().CommSplit()));
}

double MixedEllipticAssemble::FaceError(const Vector &u) const {
  double face_error = 0;  // For P0 without any pressure values on face
  return sqrt(PPM->SumOnCommSplit(face_error, u.GetMesh().CommSplit()));
}

FluxPair MixedEllipticAssemble::InflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] < 0) continue;

      const RTFaceElementT<> &faceElem = faceElementPool.Get(u, *c, face);

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Scalar FN = faceElem.VelocityField(q, u) * faceElem.QNormal(q);
        if (FN > 0) { outflow += w * FN; }
        else { inflow += w * FN; }
      }
    }
  }
  return {PPM->SumOnCommSplit(inflow, u.GetMesh().CommSplit()),
          PPM->SumOnCommSplit(outflow, u.GetMesh().CommSplit())};
}

FluxPair MixedEllipticAssemble::PrescribedInflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] < 2) continue;

      const RTFaceElementT<> &faceElem = faceElementPool.Get(u, *c, face);

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Scalar FN = problem.Flux(faceElem.QPoint(q)) *
                    faceElem.QNormal(q);
        if (FN > 0) { outflow += w * FN; }
        else { inflow += w * FN; }
      }
    }
  }
  return {PPM->SumOnCommSplit(inflow, u.GetMesh().CommSplit()),
          PPM->SumOnCommSplit(outflow, u.GetMesh().CommSplit())};
}

FluxPair MixedEllipticAssemble::OutflowLeftRight(const Vector &u) const {
  double outflowLeft = 0;
  double outflowRight = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] != 1) continue;

      const RTFaceElementT<> &faceElem = faceElementPool.Get(u, *c, face);

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        VectorField Q = faceElem.VelocityField(q, u);
        Scalar F = Q * faceElem.QNormal(q);
        if (F > 0 && c()[0] <= 1) outflowLeft += w * F;
        else if (F > 0 && c()[0] >= 2) outflowRight += w * F;
      }
    }
  }
  return {PPM->SumOnCommSplit(outflowLeft, u.GetMesh().CommSplit()),
          PPM->SumOnCommSplit(outflowRight, u.GetMesh().CommSplit())};
}

double MixedEllipticAssemble::GoalFunctional(const Vector &u) const {
  double goal = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] != 1) continue;

      const RTFaceElementT<> &faceElem = faceElementPool.Get(u, *c, face);

      for (int q = 0; q < faceElem.nQ(); ++q) {
        Point z = faceElem.QPoint(q);
        if (z[0] < 0.25) continue;
        if (z[0] > 0.5) continue;
        double w = faceElem.QWeight(q);
        VectorField Q = faceElem.VelocityField(q, u);
        Scalar F = Q * faceElem.QNormal(q);
        goal += w * F;
      }
    }
  }
  return PPM->SumOnCommSplit(goal, u.GetMesh().CommSplit());
}

void MixedEllipticAssemble::SetExactSolution(Vector &u) const {
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    vector<Point> nodalPoints = u.GetNodalPoints(elem.PressureIndex(), *c);
    MixedRowValues u_c(u, *c, elem);
    for (int i = 0; i < elem.PressureSize(); i++) {
      u_c(elem.PressureIndex(), i) = problem.Solution(nodalPoints[i]);
    }
    for (int face = 0; face < c.Faces(); ++face) {

      const RTFaceElementT<> &faceElem = faceElementPool.Get(u, *c, face);

      RowValues uf_c(u, faceElem);
      VectorField N = faceElem.Normal();
      vector<Point> nodalPointsOnFace = u.GetNodalPointsOnFace(0, *c, face);
      for (int j = 0; j < faceElem.size(); j++) {
        VectorField F = problem.Flux(nodalPointsOnFace[j]);
        uf_c(j) = faceElem.Area() * F * N;
      }
    }
  }
}

void MixedEllipticAssemble::SetPressure(const Vector &u, Vector &p) const {
  for (row r = p.rows(); r != p.rows_end(); r++) {
    row ru = u.find_row(r());
    for (int i = 0; i < r.n(); i++) {
      p(r, i) = u(ru, i);
    }
  }
};

void MixedEllipticAssemble::SetFlux(const Vector &u, Vector &flux) {
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    VectorField F = zero;
    double area = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField Q = elem.VelocityField(q, u);
      F += w * Q;
      area += w;
    }
    F *= (1 / area);
    for (int d = 0; d < SpaceDimension; ++d) {
      flux(c(), d) = F[d];
    }
  }
}

double MixedEllipticAssemble::DualPrimalError(const Vector &u) const {
  mout.StartBlock("Dual Primal Error");
  mout << endl << "Start" << endl;
  LagrangeEllipticAssemble lagrangeAssemble(problem, 1);
  Vector u_lagrange(0.0, lagrangeAssemble.GetSharedDisc(), u.Level());
  u_lagrange.SetAccumulateFlag(true);
  NewtonMethod(lagrangeAssemble, u_lagrange, true);

  double err = 0;
  for (cell c = u_lagrange.cells(); c != u_lagrange.cells_end(); ++c) {
    RTLagrangeElement rt0LagrangeElem(u, *c);
    ScalarElement scalarElem(u_lagrange, *c);
    for (int q = 0; q < rt0LagrangeElem.nQ(); ++q) {
      double w = rt0LagrangeElem.QWeight(q);
      Tensor K = problem.Permeability(rt0LagrangeElem.QPoint(q));
      Tensor IK = Invert(problem.Permeability(rt0LagrangeElem.QPoint(q)));
      VectorField K_DU = K * scalarElem.Derivative(q, u_lagrange);
      VectorField Q = rt0LagrangeElem.VelocityField(q, u);
      VectorField Diff = Q - K_DU;
      err += w * (IK * Diff) * Diff;
    }
  }
  mout.EndBlock();
  return sqrt(PPM->SumOnCommSplit(err, u.GetMesh().CommSplit()));
}

void MixedEllipticAssemble::SetPressureFlux(const Vector &u, Vector &p, Vector &flux) const {
  SetPressure(u, p);

  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);
    VectorField F = zero;
    double area = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField Q = elem.VelocityField(q, u);
      F += w * Q;
      area += w;
    }
    F *= (1 / area);
    for (int d = 0; d < SpaceDimension; ++d) {
      flux(c(), d) = F[d];
    }
  }
};

VectorField MixedEllipticAssemble::EvaluateCellFlux(const Vector &flux,
                                                    const Cell &c) const {
  RTLagrangeElementT<> elem(flux, c);

  VectorField F = zero;
  double area = 0;
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    area += w;
    F += w * elem.VelocityField(q, flux);
  }
  F *= (1 / area);
  return F;
}

double MixedEllipticAssemble::EvaluateNormalFlux(const Vector &flux, const Cell &c, int i) const {
  RTLagrangeElement elem(flux, c);
  RTFaceElement faceElem(flux, c, i);
  return flux(elem[i], 0) * elem.Sign(i) / faceElem.Area();
}

void MixedEllipticAssemble::SetNormalFlux(const Vector &u, Vector &flux) {
  flux = u;
}
