#include "HybridEllipticAssemble.hpp"
#include "LagrangeEllipticAssemble.hpp"
#include "BasicSolver.hpp"
#include "Newton.hpp"
#include "ScalarElement.hpp"


const char *HybridEllipticAssemble::Name() const {
  return "HybridEllipticAssemble";
}

void HybridEllipticAssemble::Initialize(Vector &u) const {
  elementPool.Initialize(u);
  faceElementPool.Initialize(u);

  if (problem.HasExactSolution())
    SetExactSolution(u);
  
  u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RowBndValues u_c(u, *c);
    if (!u_c.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (u_c.bc(face) == 1) {
        u_c(face) = problem.Solution(c.Face(face));
        u_c.D(face) = true;
      }
    }
  }
  u.DirichletConsistent();
}

double HybridEllipticAssemble::Energy(const Vector &u) const {
  return 0.0;
}

void HybridEllipticAssemble::Residual(const cell &c, const Vector &u, Vector &r) const {
  int Faces = c.Faces();
  RMatrix A(Faces);
  RVector B(Faces), R(Faces);
  LocalProblem(c, u, A, B, R);
  RBasicLinearSolver solver(A);
  RVector IAB = solver.Solve(B);
  Scalar BIAB = B * IAB;

  RVector RU(R);
  RVector UU(R);

  const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

  for (int face = 0; face < Faces; ++face)
    UU[face] = u(elem[face], 0);
  for (int face = 0; face < Faces; ++face)
    RU[face] *= UU[face];
  RVector IARU = solver.Solve(RU);
  Scalar BIARU = B * IARU;
  Scalar P = BIARU / BIAB;

  B *= P;
  RU -= B;
  RVector FF = solver.Solve(RU);
  RowBndValues r_c(r, *c);
  for (int face = 0; face < Faces; ++face)
    r_c(face) += R[face] * FF[face];
  if (!r_c.onBnd()) return;
  for (int face = 0; face < Faces; ++face) {
    int bnd = r_c.bc(face);
    if (bnd != 2) continue;
    RTFaceElementT<> faceElem(u, *c, face);
    for (int q = 0; q < faceElem.nQ(); ++q) {
      double w = faceElem.QWeight(q);
      Scalar F = problem.Flux(faceElem.QPoint(q)) * faceElem.QNormal(q);
      r_c(face) -= w * F;
    }
  }
}

void HybridEllipticAssemble::Jacobi(const cell &c, const Vector &u, Matrix &J) const {
  int Faces = c.Faces();
  RMatrix A(Faces);
  RVector B(Faces), R(Faces);
  LocalProblem(c, u, A, B, R);
  RMatrix IA = invert(A);
  RVector IAB = IA * B;
  Scalar BIAB = B * IAB;


  const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

  RowEntries J_c(J, elem);
  for (int face_i = 0; face_i < Faces; ++face_i)
    for (int face_j = 0; face_j < Faces; ++face_j)
      J_c(face_i, face_j) += R[face_i] * (IA[face_i][face_j] -
                                          IAB[face_i] * IAB[face_j] / BIAB) * R[face_j];
}

void HybridEllipticAssemble::LocalProblem(const cell &c, const Vector &u,
                                          RMatrix &A, RVector &B, RVector &R) const {
  A = 0, B = 0, R = 0;

  const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    Tensor IK = Invert(problem.Permeability(elem.QPoint(q)));
    for (int i = 0; i < elem.VelocitySize(); ++i) {
      for (int k = 0; k < elem.VelocityMaxk(i); ++k) {
        Velocity IK_U_i = IK * elem.VelocityField(q, i, k);
        Scalar divU_i = elem.VelocityFieldDivergence(q, i, k);
        B[i] += w * divU_i;
        for (int j = 0; j < c.Faces(); ++j) {
          Velocity U_j = elem.VelocityField(q, j, k);
          A[i][j] += w * IK_U_i * U_j;
        }
      }
    }
  }
  for (int face = 0; face < c.Faces(); ++face) {
    RTFaceElementT<> faceElem(u, *c, face);
    for (int q = 0; q < faceElem.nQ(); ++q) {
      double w = faceElem.QWeight(q);
      for (int i = 0; i < faceElem.size(); i++) {
        Scalar face_U_i = faceElem.VelocityField(q, i) * faceElem.QNormal(q);
        R[face] += w * face_U_i;
      }
    }
  }
}

std::pair<double, double> HybridEllipticAssemble::InflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] < 0) continue;
      RTFaceElementT<> faceElem(u, *c, face);
      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Scalar F = FaceFlux(c, u, faceElem, q, face) *
          faceElem.QNormal(q);
        if (F > 0) outflow += w * F;
        else inflow += w * F;
      }
    }
  }
  return {PPM->SumOnCommSplit(inflow, u.CommSplit()),
          PPM->SumOnCommSplit(outflow, u.CommSplit())};
}

FluxPair HybridEllipticAssemble::OutflowLeftRight(const Vector &p) const {
  double outflowLeft = 0;
  double outflowRight = 0;
  for (cell c = p.cells(); c != p.cells_end(); ++c) {
    BFParts bnd(p.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] < 0) continue;

      const RTFaceElementT<> &faceElem = faceElementPool.Get(p, *c, face);

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Scalar F = FaceFlux(c, p, faceElem, q, face) *
                   faceElem.QNormal(q);
        if (F > 0 && c()[0] <= 1) outflowLeft += w * F;
        else if (F > 0 && c()[0] >= 2) outflowRight += w * F;
      }
    }
  }
  return {PPM->SumOnCommSplit(outflowLeft, p.CommSplit()),
          PPM->SumOnCommSplit(outflowRight, p.CommSplit())};
}

//void HybridEllipticAssemble::SetPressure(const Vector &u, Vector &p) const {
//  for (cell c = u.cells(); c != u.cells_end(); ++c) {
//    int N = c.Faces();
//    RMatrix A(N);
//    RVector B(N);
//    RVector R(N);
//    LocalProblem(c, u, A, B, R);
//    RBasicLinearSolver solver(A);
//    RVector IAB = solver.Solve(B);
//    Scalar BIAB = B * IAB;
//    RVector RU(R);
//    RVector UU(R);
//    RTLagrangeElementT<> elem(u, *c);
//    for (int i = 0; i < c.Faces(); ++i)
//      UU[i] = u(elem[i], 0);
//    for (int i = 0; i < c.Faces(); ++i)
//      RU[i] *= UU[i];
//    RVector IARU = solver.Solve(RU);
//    Scalar BIARU = B * IARU;
//    p(c(), 0) = BIARU / BIAB;
//  }
//}

void HybridEllipticAssemble::SetFlux(const Vector &u, Vector &flux) {
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    int N = c.Faces();
    RMatrix A(N);
    RVector B(N), R(N);
    LocalProblem(c, u, A, B, R);
    RBasicLinearSolver solver(A);
    RVector IAB = solver.Solve(B);
    Scalar BIAB = B * IAB;
    RVector RU(R);
    RVector UU(R);

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    for (int i = 0; i < c.Faces(); ++i)
      UU[i] = u(elem[i], 0);
    for (int i = 0; i < c.Faces(); ++i)
      RU[i] *= UU[i];
    RVector IARU = solver.Solve(RU);
    Scalar BIARU = B * IARU;
    Scalar P = BIARU / BIAB;
    B *= P;
    RU -= B;
    RVector FF = solver.Solve(RU);
    VectorField F = zero;
    double area = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      area += w;
      for (int i = 0; i < elem.VelocitySize(); ++i) {
        for (int k = 0; k < elem.VelocityMaxk(i); ++k) {
          VectorField U_i = elem.VelocityField(q, i, k);
          F += w * FF[i] * U_i;
        }
      }
    }
    F *= (1 / area);
    for (int d = 0; d < SpaceDimension; ++d)
      flux(c(), d) = F[d];
  }
}

void HybridEllipticAssemble::SetPressureFlux(const Vector &u, Vector &p, Vector &flux) const {
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    int N = c.Faces();
    RMatrix A(N);
    RVector B(N), R(N);
    LocalProblem(c, u, A, B, R);
    RBasicLinearSolver solver(A);
    RVector IAB = solver.Solve(B);
    Scalar BIAB = B * IAB;
    RVector RU(R);
    RVector UU(R);

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    double p_c = 0;
    for (int i = 0; i < c.Faces(); ++i) {
      UU[i] = u(elem[i], 0);
      p_c += UU[i];
    }
    p(c(),0) = p_c / c.Faces();    
    for (int i = 0; i < c.Faces(); ++i)
      RU[i] *= UU[i];
    RVector IARU = solver.Solve(RU);
    Scalar BIARU = B * IARU;
    Scalar P = BIARU / BIAB;
    B *= P;
    RU -= B;
    RVector FF = solver.Solve(RU);
    VectorField F = zero;
    double area = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      area += w;
      for (int i = 0; i < elem.VelocitySize(); ++i) {
        for (int k = 0; k < elem.VelocityMaxk(i); ++k) {
          VectorField U_i = elem.VelocityField(q, i, k);
          F += w * FF[i] * U_i;
        }
      }
    }
    F *= (1 / area);
    for (int d = 0; d < SpaceDimension; ++d)
      flux(c(), d) = F[d];
  }
}

void HybridEllipticAssemble::SetExactSolution(Vector &u) const {
  for (face f = u.faces(); f != u.faces_end(); ++f)
    u(f(), 0) = problem.Solution(f());
}

//void SetSolution (const Vector& u, Vector& p) {
//  for (cell c = u.cells(); c != u.cells_end(); ++c) {
//    int N = c.Faces();
//    SmallMatrix A (N);
//    SmallVector B (N);
//    SmallVector R (N);
//    LocalProblem (c, u, A, B, R);
//    InverseSmallMatrix IA (A);
//    SmallVector IAB (IA, B);
//    Scalar  BIAB = B * IAB;
//    SmallVector RU (R);
//    SmallVector UU (R);
//    RT0_P0Element E (disc, u, *c);
//    for (int i = 0; i < c.Faces(); ++i)
//      UU[i] = u (E[i], 0);
//    for (int i = 0; i < c.Faces(); ++i)
//      RU[i] *= UU[i];
//    SmallVector IARU (IA, RU);
//    Scalar BIARU = B * IARU;
//    p (c(), 0) = BIARU / BIAB;
//  }
//}

void HybridEllipticAssemble::SetNormalFlux(const Vector &u, Vector &flux) {
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    int Faces = c.Faces();
    RMatrix A(Faces);
    RVector B(Faces), R(Faces);
    LocalProblem(c, u, A, B, R);
    RBasicLinearSolver solver(A);
    RVector IAB = solver.Solve(B);
    Scalar BIAB = B * IAB;
    RVector RU(R);
    RVector UU(R);

    const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

    for (int face = 0; face < c.Faces(); ++face)
      UU[face] = u(elem[face], 0);
    for (int face = 0; face < c.Faces(); ++face)
      RU[face] *= UU[face];
    RVector IARU = solver.Solve(RU);
    Scalar BIARU = B * IARU;
    Scalar P = BIARU / BIAB;
    B *= P;
    RU -= B;
    RVector F = solver.Solve(RU);
    for (unsigned int i = 0; i < elem.VelocitySize(); ++i)
      flux(elem[i], 0) = F[i];
  }
}

VectorField HybridEllipticAssemble::EvaluateCellFlux(const Vector &flux,
                                                     const Cell &c) const {
  if constexpr (DebugLevel > 0) {
    if (flux.find_cell(c()) == flux.cells_end()) {
      mout << PPM->Proc() << " cell " << c() << "not found" << endl;
      //std::cout << PPM->Proc() << "space_center:" << DOUT(c.Center()) << endl;
      THROW("HybridEllipticAssemble::EvaluateCellFlux: cell not found.")
    }
  }

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

double HybridEllipticAssemble::EvaluateNormalFlux(const Vector &flux, const Cell &c, int i) const {
  if constexpr (DebugLevel > 0) {
    if (flux.find_cell(c()) == flux.cells_end()) {
      mout << PPM->Proc() << " cell " << c() << "not found" << endl;
      //std::cout << PPM->Proc() << "space_center:" << DOUT(c.Center()) << endl;
      THROW("HybridEllipticAssemble::EvaluateNormalFlux: cell not found.")
    }
  }

  RTLagrangeElementT<> elem(flux, c);
  RTFaceElementT<> faceElem(flux, c, i);
  return flux(elem[i], 0) * elem.Sign(i) / faceElem.Area();
}


VectorField HybridEllipticAssemble::EvaluateNormal(const Vector &flux,
                                                   const cell &c, int i) const {


  const RTLagrangeElementT<> &elem = elementPool.Get(flux, *c);

  return elem.OuterNormal(i);
}

Scalar HybridEllipticAssemble::Value(const cell &c, const Vector &u,
                                     RTLagrangeElementT<> &E) const {
  int Faces = c.Faces();
  RMatrix A(Faces);
  RVector B(Faces), R(Faces);
  LocalProblem(c, u, A, B, R);
  RBasicLinearSolver solver(A);
  RVector IAB = solver.Solve(B);
  Scalar BIAB = B * IAB;
  RVector RU(R);
  RVector UU(R);

  const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

  for (int face = 0; face < c.Faces(); ++face)
    UU[face] = u(elem[face], 0);
  for (int face = 0; face < c.Faces(); ++face)
    RU[face] *= UU[face];
  RVector IARU = solver.Solve(RU);
  Scalar BIARU = B * IARU;
  return BIARU / BIAB;
}

Scalar HybridEllipticAssemble::FaceValue(const cell &c, const Vector &u,
                                         RTFaceElementT<> &FE) const {
  int Faces = c.Faces();
  RMatrix A(Faces);
  RVector B(Faces), R(Faces);
  LocalProblem(c, u, A, B, R);
  RBasicLinearSolver solver(A);
  RVector IAB = solver.Solve(B);
  Scalar BIAB = B * IAB;
  RVector RU(R);
  RVector UU(R);

  const RTLagrangeElementT<> &elem = elementPool.Get(u, *c);

  for (int face = 0; face < c.Faces(); ++face)
    UU[face] = u(elem[face], 0);
  for (int face = 0; face < Faces; ++face)
    RU[face] *= UU[face];
  RVector IARU = solver.Solve(RU);
  Scalar BIARU = B * IARU;
  return BIARU / BIAB;
}

VectorField HybridEllipticAssemble::Flux(const cell &c, const Vector &u,
                                         RTLagrangeElementT<> &elem, int q) const {
  int Faces = c.Faces();
  RMatrix A(Faces);
  RVector B(Faces), R(Faces);
  LocalProblem(c, u, A, B, R);
  RBasicLinearSolver solver(A);
  RVector IAB = solver.Solve(B);
  Scalar BIAB = B * IAB;
  RVector RU(R);
  RVector UU(R);
  for (int face = 0; face < Faces; ++face)
    UU[face] = u(elem[face], 0);
  for (int face = 0; face < Faces; ++face)
    RU[face] *= UU[face];
  RVector IARU = solver.Solve(RU);
  Scalar BIARU = B * IARU;
  Scalar P = BIARU / BIAB;
  B *= P;
  RU -= B;
  RVector FF = solver.Solve(RU);
  VectorField F = zero;
  for (int face = 0; face < Faces; ++face) {
    VectorField U_i = elem.VelocityField(q, face, 0);
    F += FF[face] * U_i;
  }
  return F;
}

VectorField HybridEllipticAssemble::FaceFlux(const cell &c, const Vector &u,
                                             const RTFaceElementT<> &faceElem,
                                             int q, int face) const {
  RTLagrangeElementT<> elem(u, *c);
  int Faces = c.Faces();
  RMatrix A(Faces);
  RVector B(Faces), R(Faces);
  LocalProblem(c, u, A, B, R);
  RBasicLinearSolver solver(A);
  RVector IAB = solver.Solve(B);
  Scalar BIAB = B * IAB;
  RVector RU(R);
  RVector UU(R);
  for (int f = 0; f < Faces; ++f)
    UU[f] = u(elem[f], 0);
  for (int f = 0; f < Faces; ++f)
    RU[f] *= UU[f];
  RVector IARU = solver.Solve(RU);
  Scalar BIARU = B * IARU;
  Scalar P = BIARU / BIAB;
  B *= P;
  RU -= B;
  RVector FF = solver.Solve(RU);
  return FF[face] * faceElem.VelocityField(q, face);
}

double HybridEllipticAssemble::EnergyError(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RTLagrangeElementT<> elem(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point z = elem.QPoint(q);
      Tensor IK = Invert(problem.Permeability(z));
      VectorField diff = (Flux(c, u, elem, q) -
        problem.Flux(z));
      err += w * (IK * diff) * diff;
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

double HybridEllipticAssemble::L2(const Vector &u) const {
  double l2 = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RTLagrangeElementT<> elem(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar U = Value(c, u, elem);
      l2 += w * U * U;
    }
  }
  return sqrt(PPM->SumOnCommSplit(l2, u.CommSplit()));
}

double HybridEllipticAssemble::H1(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RTLagrangeElementT<> elem(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar U = Value(c, u, elem);
      Tensor IK = Invert(problem.Permeability(elem.QPoint(q)));
      VectorField DU = IK * Flux(c, u, elem, q);
      err += w * U * U + w * DU * DU;
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

double HybridEllipticAssemble::L2Error(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RTLagrangeElementT<> elem(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar U = Value(c, u, elem);
      Scalar Sol = problem.Solution(elem.QPoint(q));
      err += w * (U - Sol) * (U - Sol);
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

double HybridEllipticAssemble::L2CellAvgError(const Vector &u) const {
  double err = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RTLagrangeElementT<> elem(u, *c);
    double w = 0;
    Scalar U = 0;
    Scalar Sol = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      w += elem.QWeight(q);
      U += Value(c, u, elem);
      Sol += problem.Solution(elem.QPoint(q));
    }
    err += w * (U - Sol) * (U - Sol);
  }
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

double HybridEllipticAssemble::MaxError(const Vector &u) const {
  double err = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RTLagrangeElementT<> elem(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      Scalar U = Value(c, u, elem);
      err = max(err, abs(U - problem.Solution(elem.QPoint(q))));
    }
  }
  return PPM->Max(err, u.CommSplit());
}

double HybridEllipticAssemble::FluxError(const Vector &u) const {
  double flux_error = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] == 2) {
        RTFaceElement faceElem(u, *c, face);
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          VectorField G = problem.Flux(faceElem.QPoint(q));
          Scalar F = (FaceFlux(c, u, faceElem, q, face) - G) * faceElem.QNormal(q);
          flux_error += w * F * F;
        }
      } else if (bnd[face] == 0) {
        RTFaceElement faceElem(u, *c, face);
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Scalar F = FaceFlux(c, u, faceElem, q, face) * faceElem.QNormal(q);
          flux_error += w * F * F;
        }
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(flux_error, u.CommSplit()));
}

double HybridEllipticAssemble::FaceError(const Vector &u) const {
  double face_error = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    for (int face = 0; face < c.Faces(); ++face) {
      RTFaceElement faceElem(u, *c, face);
      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Scalar U = FaceValue(c, u, faceElem);
        Scalar Sol = problem.Solution(faceElem.QPoint(q));
        face_error += w * (U - Sol) * (U - Sol);
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(face_error, u.CommSplit()));
}

FluxPair HybridEllipticAssemble::PrescribedInflowOutflow(const Vector &u) const {
  double inflow = 0;
  double outflow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] < 2) continue;
      RTFaceElement faceElem(u, *c, face);
      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Scalar F = problem.Flux(faceElem.QPoint(q)) * faceElem.QNormal(q);
        if (F > 0) { outflow += w * F; }
        else { inflow += w * F; }
      }
    }
  }
  return {PPM->SumOnCommSplit(inflow, u.CommSplit()),
          PPM->SumOnCommSplit(outflow, u.CommSplit())};
}

double HybridEllipticAssemble::GoalFunctional(const Vector &u) const {
  double goal = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] != 1) continue;
      RTFaceElement faceElem(u, *c, face);
      for (int q = 0; q < faceElem.nQ(); ++q) {
        Point z = faceElem.QPoint(q);
        if (z[0] < 0.25) continue;
        if (z[0] > 0.5) continue;
        double w = faceElem.QWeight(q);
        Scalar F = FaceFlux(c, u, faceElem, q, face) * faceElem.QNormal(q);
        goal += w * F;
      }
    }
  }
  return PPM->SumOnCommSplit(goal, u.CommSplit());
}

double HybridEllipticAssemble::DualPrimalError(const Vector &u) const {
  mout.StartBlock("Dual Primal Error");
  mout << endl <<"Start" << endl;
  LagrangeEllipticAssemble lagrangeAssemble(problem, 1);
  Vector u_lagrange(0.0, lagrangeAssemble.GetSharedDisc(), u.Level());
  u_lagrange.SetAccumulateFlag(true);
  NewtonMethod(lagrangeAssemble, u_lagrange, true);

  double err = 0;
  for (cell c = u_lagrange.cells(); c != u_lagrange.cells_end(); ++c) {
    RTLagrangeElementT<> rtElem(u, *c);
    ScalarElement scalarElem(u_lagrange, *c);
    for (int q = 0; q < rtElem.nQ(); ++q) {
      double w = rtElem.QWeight(q);
      Tensor K = problem.Permeability(rtElem.QPoint(q));
      Tensor IK = Invert(problem.Permeability(rtElem.QPoint(q)));
      VectorField K_DU = K * scalarElem.Derivative(q, u_lagrange);
      VectorField U = Flux(c, u, rtElem, q) - K_DU;
      err += w * (IK * U) * U;
    }
  }
  mout.EndBlock();
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}
