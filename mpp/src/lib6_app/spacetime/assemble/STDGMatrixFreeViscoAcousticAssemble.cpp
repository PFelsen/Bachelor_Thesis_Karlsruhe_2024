#include "STDGMatrixFreeViscoAcousticAssemble.hpp"

#include <tuple>

#include "DebuggingTools.hpp"
#include "LagrangeDiscretization.hpp"
#include "LinearSolver.hpp"
#include "STDGDGViscoAcousticElement.hpp"
#include "STDGDGViscoAcousticFaceElement.hpp"
#include "ScalarElement.hpp"
#include "SymRMatrix.hpp"
#include "systems/CalculationInfo.hpp"
#include "systems/IterableSystem.hpp"

inline Scalar __internal_VFLUX_CF(const VectorField &V, const VectorField &TN,
                                  bool bnd, int bnd_id = 2) {
  Scalar TNV = TN * V;
  if (!bnd)
    return TNV;
  else {
    if (bnd_id == 1)
      return TNV;  // Dirichlet bnd
    if (bnd_id == 3)
      return -TNV;  // Robin bnd
    return -TNV;    // Neumann bnd
  }
}

inline Scalar __internal_PFLUX_CF(const double &P, bool bnd, int bnd_id = 2) {
  if (!bnd)
    return P;
  else {
    if (bnd_id == 1)
      return -P;  // Dirichlet bnd
    if (bnd_id == 3)
      return -P;  // Robin bnd
    return P;     // Neumann bnd
  }
}

inline void CalculateFirstLoopIntern(CalculationInfo &info) {
  auto &element = info.element;
  auto problem = info.problem;
  const auto &kappaInverse = info.kappaInverse;
  const auto &kappaTauInverse = info.kappaTauInverse;
  auto &rowEntries = info.rowEntries;
  const Cell &c = *info.cellIterator;
  if (problem->HasRHS()) {
    const auto &cellPoint = info.point;
    auto &vector = info.vector;
    for (auto &[index, weight, point] : quadratureIterableWithPoint(info)) {
      const Scalar generalPressure = problem->F_Pressure(point.t(), c, point);
      const VectorField &generalVelocity = problem->F_Velocity(point.t(), c, point);
      const DampingVector &dampingVector = problem->F_DampingPressure(point.t(), c, point);
      for (auto &[out, pressure] :
           pressureIterable<true>(index, info,
                                  VectorIterator(vector, cellPoint))) {
        int __variableIInternal = element.variable(out.i);  // Todo remove
        if (__variableIInternal > 1) {
          *out += weight * dampingVector[__variableIInternal - 2] * pressure;
        } else if (__variableIInternal == 1) {
          *out += weight * generalPressure * pressure;
        }
      }
      for (auto &[out, velocity] :
           velocityIterable<true>(index, info,
                                  VectorIterator(vector, cellPoint))) {
        int __variableIInternal = element.variable(out.i);  // Todo remove
        if (__variableIInternal == 0) {
          *out += weight * generalVelocity * velocity;
        }
      }
    }
  }

  for (auto [index, weight] : quadratureIterable(info)) {
    Date __debugStartPressure;
    for (auto [next, pressure] : pressureIterable<true>(index, info)) {
      int __variableIInternal = element.variable(next.i);  // TODO Remove
      for (auto [out, derivativePressure] :
           derivativePressureIterable<true>(index, info, next)) {
        int __variableJInternal = element.variable(out.j);  // TODO Remove
        if (__variableIInternal != __variableJInternal)
          continue;
        *out += weight * kappaInverse[__variableJInternal] *
                (derivativePressure * pressure);
      }
      for (auto [out, pressureOther] :
           pressureIterable<true>(index, info, next)) {
        int __variableJInternal = element.variable(out.j);  // TODO Remove
        if (__variableJInternal <= 1 ||
            __variableIInternal != __variableJInternal)
          continue;
        *out += weight * kappaTauInverse[__variableJInternal] *
                (pressure * pressureOther);
      }
      for (auto [out, divergentVelocity] :
           divergenceVelocityIterable<true>(index, info, next)) {
        *out -= weight * (divergentVelocity * pressure);
      }
    }

    const double rho = info.rho;
    for (auto &[next, velocity] : velocityIterable<true>(index, info)) {
      for (auto &[out, gradientPressure] :
           gradientPressureIterable<true>(index, info, next)) {
        *out -= weight * (gradientPressure * velocity);
      }
      for (auto &[out, derivativeVelocity] :
           derivativeVelocityIterable<true>(index, info, next)) {
        *out += weight * rho * (derivativeVelocity * velocity);
      }
    }
  };
}

inline void CalculateSecondLoopIntern(CalculationInfo &info) {
  const auto &point = info.point;
  const auto &cellIterator = info.cellIterator;
  auto problem = info.problem;
  const double z_K = sqrt(problem->Rho(*cellIterator, point) * problem->Kappa(*cellIterator, point));
  auto &vector = info.vector;
  auto &matrix = info.matrix;
  const auto &mesh = matrix.GetMesh();
  auto disc = info.discretization;
  auto &rowEntries = info.rowEntries;

  std::vector<Scalar> boundaryCache;
  for (size_t faceID = 0; faceID < cellIterator.Faces() - 2; ++faceID) {
    CalculateFaceInfo faceInfo(info, faceID, matrix, z_K);

    if (faceInfo.boundaryID <= 2 && faceInfo.boundaryID >= 1) {
      boundaryCache.resize(faceInfo.quadratureSize);
      double pressureAlpha = 1;
      double velocityAlpha = faceInfo.z_K;
      if (faceInfo.boundaryID == 1) {
        pressureAlpha = 1.0 / faceInfo.z_K;
        velocityAlpha = 1;
        for (auto &[index, weight, quadraturePoint] :
             quadratureIterableWithPoint(faceInfo)) {
          boundaryCache[index] =
              weight *
              (info.problem->ut_Pressure(quadraturePoint.t(), quadraturePoint,
                                         *info.cellIterator) +
               info.problem
                   ->ut_DampingPressure(quadraturePoint.t(), quadraturePoint, *info.cellIterator)
                   .sum());
        }
      } else {
        for (auto &[index, weight, quadraturePoint] :
             quadratureIterableWithPoint(faceInfo)) {
          const VectorField &Nq = faceInfo.element.qNormal[index];
          boundaryCache[index] =
              weight * (Nq * info.problem->ut_Velocity(quadraturePoint.t(), quadraturePoint,
                                                       *info.cellIterator));
        }
      }

      for (size_t index = 0; index < faceInfo.quadratureSize; index++) {
        for (auto &[out, pressure] :
             pressureIterable(index, faceInfo, VectorIterator(vector, point))) {
          *out += (pressureAlpha * pressure) * boundaryCache[index];
        }
        for (auto &[out, flux] :
             fluxIterable(index, faceInfo, VectorIterator(vector, point))) {
          *out += (velocityAlpha * flux) * boundaryCache[index];
        }
      }
    }

    auto &faceRowEntries = faceInfo.rowEntries;
    for (auto &[quadrature, weight, point] :
         quadratureIterableWithPoint(faceInfo)) {
      int q1 = find_q_id(point, faceInfo.secondElementInfo.element);
      const VectorField &Nq =
          faceInfo.element.qNormal[quadrature];  // TODO remove
      for (auto [nextPack, generalFluxVelocity] :
           fluxIterable(quadrature, faceInfo,
                        zipOutput(&rowEntries, &faceRowEntries))) {
        auto &[next, nextFace] = nextPack;
        for (auto [out, currentFluxVelocity] :
             fluxIterable(quadrature, faceInfo, next)) {
          *out += weight * faceInfo.alpha[1] * generalFluxVelocity *
                  currentFluxVelocity;
        }
        for (auto [j, pressure] :
             pressureIterable(quadrature, faceInfo, next)) {
          *j += weight * faceInfo.alpha[3] * generalFluxVelocity * pressure;
        }

        for (auto &[j, velocity] :
             velocityIterable(q1, faceInfo.secondElementInfo, nextFace)) {
          double VFlux_cf_V1j =
              __internal_VFLUX_CF(velocity, Nq, faceInfo.isOnBoundary,
                                  faceInfo.boundaryID);  // TODO remove
          *j -= weight * faceInfo.alpha[1] * VFlux_cf_V1j * generalFluxVelocity;
        }
        for (auto [j, pressureFace] :
             pressureIterable(q1, faceInfo.secondElementInfo, nextFace)) {
          double PFlux_cf_P1j =
              __internal_PFLUX_CF(pressureFace, faceInfo.isOnBoundary,
                                  faceInfo.boundaryID);  // TODO remove
          *j -= weight * faceInfo.alpha[3] * PFlux_cf_P1j * generalFluxVelocity;
        }
      }

      for (auto [nextPack, generalPressure] :
           pressureIterable(quadrature, faceInfo,
                            zipOutput(&rowEntries, &faceRowEntries))) {
        auto &[next, nextFace] = nextPack;
        for (auto [out, currentFluxVelocity] :
             fluxIterable(quadrature, faceInfo, next)) {
          *out += weight * faceInfo.alpha[2] * generalPressure *
                  currentFluxVelocity;
        }
        for (auto [out, pressure] :
             pressureIterable(quadrature, faceInfo, next)) {
          *out += weight * faceInfo.alpha[0] * generalPressure * pressure;
        }

        for (auto &[out, velocity] :
             velocityIterable(q1, faceInfo.secondElementInfo, nextFace)) {
          double VFlux_cf_V1j =
              __internal_VFLUX_CF(velocity, Nq, faceInfo.isOnBoundary,
                                  faceInfo.boundaryID);
          *out -= weight * faceInfo.alpha[2] * VFlux_cf_V1j * generalPressure;
        }
        for (auto [out, pressureFace] :
             pressureIterable(q1, faceInfo.secondElementInfo, nextFace)) {
          double PFlux_cf_P1j =
              __internal_PFLUX_CF(pressureFace, faceInfo.isOnBoundary,
                                  faceInfo.boundaryID);
          *out -= weight * faceInfo.alpha[0] * PFlux_cf_P1j * generalPressure;
        }
      }
    }
  }
}

inline void CalculateThirdLoopIntern(CalculationInfo &info) {
  auto &matrix = info.matrix;
  const auto &cellIterator = info.cellIterator;
  const auto &problem = info.problem;
  const auto &discretization = info.discretization;
  auto &rowEntries = info.rowEntries;
  const auto rho = info.rho;
  cell c_prev = matrix.find_previous_cell(cellIterator);
  SpaceTimeViscoAcousticDGTFaceElement felem(*discretization, matrix,
                                             cellIterator,
                                             cellIterator.Faces() - 2,
                                             problem->nL(), "dg");

  DGRowEntries M_c_prev(matrix, *cellIterator, *c_prev, false);

  ElementWeightCacheable cacheable(felem, info.weightFunction, rowEntries);

  for (auto &[quadrature, weight, point] :
       quadratureIterableWithPoint(cacheable)) {
    for (auto &[next, velocity] : velocityIterable(quadrature, cacheable)) {
      int var_i = felem.variable(next.i);
      for (auto &[out, otherVelocity] :
           velocityIterable(quadrature, cacheable, next)) {
        int var_j = felem.variable(out.j);
        if (var_i != var_j)
          continue;
        *out += weight * rho * velocity * otherVelocity;
      }
    }
    for (auto &[next, pressure] : pressureIterable(quadrature, cacheable)) {
      int var_i = felem.variable(next.i);
      for (auto &[out, otherPressure] :
           pressureIterable(quadrature, cacheable, next)) {
        int var_j = felem.variable(out.j);
        if (var_i != var_j)
          continue;
        *out += weight * info.kappaInverse[var_j] * otherPressure * pressure;
      }
    }
  }

  if (cellIterator.min() == 0.0) {
    for (auto &[quadrature, weight, point] :
         quadratureIterableWithPoint(cacheable)) {
      VectorField generalVelocity = problem->ut_Velocity(point.t(), point, *cellIterator);
      for (auto &[out, velocity] :
           velocityIterable(quadrature, cacheable, VectorIterator(info.vector, info.point))) {
        int var_i = felem.variable(out.i);
        if (var_i != 0)
          continue;
        *out += weight * rho * generalVelocity * velocity;
      }

      DampingVector dampingVector = problem->ut_DampingPressure(point.t(), point, *cellIterator);
      double generalPressure = problem->ut_Pressure(point.t(), point, *cellIterator);
      for (auto [i, pressure] : pressureIterable(quadrature, cacheable, size_t{0})) {
        int var_i = felem.variable(i);
        if (var_i >= 2) {
          info.vector(info.point, i) += weight * dampingVector[var_i - 2] * pressure;
        } else if (var_i == 1) {
          info.vector(info.point, i) += weight * info.kappaInverse[var_i] * generalPressure * pressure;
        }
      }
    }
  } else if (cellIterator != c_prev) {
    SpaceTimeViscoAcousticDGTFaceElement felem_1(*discretization, matrix,
                                                 c_prev, c_prev.Faces() - 1,
                                                 problem->nL(), "dg",
                                                 cellIterator);
    ElementHolder holder(felem_1, M_c_prev);
    for (auto &[quadrature, weight, point] : quadratureIterableWithPoint(cacheable)) {
      OutputIterator<false> output(&M_c_prev);
      for (auto &[next, velocity] : velocityIterable(quadrature, cacheable, output)) {
        int var_i = felem_1.variable(next.i);
        for (auto &[out, otherVelocity] : velocityIterable(quadrature, holder, next)) {
          int var_j = felem_1.variable(out.j);
          if (var_i != var_j)
            continue;
          *out -= weight * rho * otherVelocity * velocity;
        }
      }
      for (auto [next, pressure] :
           pressureIterable(quadrature, cacheable, output)) {
        int var_i = felem_1.variable(next.i);
        for (auto [out, otherPressure] :
             pressureIterable(quadrature, holder, next)) {
          int var_j = felem_1.variable(out.j);
          if (var_i != var_j)
            continue;
          *out -= weight * info.kappaInverse[var_j] * otherPressure * pressure;
        }
      }
    }
  }
}

using namespace std;

void STDGMatrixFreeViscoAcousticAssemble::MassMatrix(Matrix &M) const {
  M = 0.0;
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(M, c, prob->nL());
    DGRowEntries M_c(M, *c, *c, false);
    double rho = prob->Rho(*c, c());
    vector<double> kappaInv(2 + prob->nL());
    vector<double> kappaTauInv(2 + prob->nL());
    kappaInv[1] = 1.0 / prob->Kappa_i(*c, c(), 0);
    for (int j = 2; j < 2 + prob->nL(); j++) {
      kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j - 1);
      kappaTauInv[j] = kappaInv[j] / prob->Tau_i(*c, c(), j - 1);
    }
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      for (int i = 0; i < elem.i_dimension(); ++i) {
        int var_i = elem.variable(i);
        VectorField V_i = elem.Velocity(q, i);
        double P_i = elem.Pressure(q, i);
        for (int j = 0; j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          if (var_j != var_i)
            continue;
          VectorField V_j = elem.Velocity(q, j);
          Scalar P_j = elem.Pressure(q, j);
          M_c(i, j) += w * (rho * (V_j * V_i) + kappaInv[var_j] * (P_j * P_i));
        }
      }
    }
  }
}

double STDGMatrixFreeViscoAcousticAssemble::MhalfNorm(const Vector &u) const {
  double normsq = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); c++) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());

    //double rho = prob->Rho(c());
    double rho_inv = 1.0 / prob->Rho(*c, c());
    double kappa = prob->Kappa(*c, c());
    double kappa_inv = 1.0 / prob->Kappa(*c, c());
    for (int q = 0; q < elem.nQ(); q++) {
      Point QP = elem.QPoint(q);
      double w = elem.QWeight(q);

      double dtP = elem.DtPressure(q, u);
      VectorField gradP = elem.GradPressure(q, u);
      VectorField dtV = elem.DtVelocity(q, u);
      double divV = elem.DivVelocity(q, u);

      double P = kappa * dtP - divV;
      VectorField V = rho_inv * dtV - gradP;

      normsq += w * (kappa * P * P + rho_inv * V * V);
    }
  }
  return sqrt(PPM->SumOnCommSplit(normsq, u.CommSplit()));
}

double STDGMatrixFreeViscoAcousticAssemble::MhalfInvLNorm(
    const Vector &u) const {
  double normsq = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); c++) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());

    double rho = prob->Rho(*c, c());
    double rho_inv = 1.0 / prob->Rho(*c, c());
    double kappa = prob->Kappa(*c, c());
    double kappa_inv = 1.0 / prob->Kappa(*c, c());
    for (int q = 0; q < elem.nQ(); q++) {
      Point QP = elem.QPoint(q);
      double w = elem.QWeight(q);

      double dtP = elem.DtPressure(q, u);
      VectorField gradP = elem.GradPressure(q, u);
      VectorField dtV = elem.DtVelocity(q, u);
      double divV = elem.DivVelocity(q, u);

      double P = kappa * dtP - divV;
      VectorField V = rho_inv * dtV - gradP;

      normsq += w * (kappa_inv * P * P + rho * V * V);
    }
  }
  return sqrt(PPM->SumOnCommSplit(normsq, u.CommSplit()));
}

void STDGMatrixFreeViscoAcousticAssemble::System(Matrix &M, Vector &RHS) const {
  Date Start;
  mout.PrintInfo("STDGMatrixFreeViscoAcousticAssemble", verbose,
                 PrintInfoEntry<int>("Problem Size", M.pSize(), 2),
                 PrintInfoEntry<Date>("start assemble", Start, 2));
  Time t_cell;
  Time t_face;
  const double T = M.GetMesh().GetEndTime();
  std::function<double(Point)> d_T = [T](Point QP) { return (T - QP.t()); };
  static constexpr auto constOne_func = [](Point QP) { return 1.0; };
  bool WeightedAssemble = false;
  Config::Get("WeightedAssemble", WeightedAssemble);
  std::function<double(Point)> weight_func =
      WeightedAssemble ? d_T : constOne_func;

  M = 0;

  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    Date Start_cell;

    CalculationInfo info(M, prob, RHS, c, weight_func, disc);
    auto &elem = info.element;
    const double rho = info.rho;
    auto &M_c = info.rowEntries;

    if (Point global_src; prob->HasRHS() && prob->hasPointSource(global_src)) {
      if (c.PointInCell(global_src)) {
        for (int i = 0; i < elem.i_dimension(); ++i) {
          RHS(c(), i) += elem.PressureGlobal(global_src, i);
        }
      }
    }

    CalculateFirstLoopIntern(info);

    t_cell += Date() - Start_cell;
    Date Start_face;
    CalculateSecondLoopIntern(info);

    CalculateThirdLoopIntern(info);
    t_face += Date() - Start_face;
  }
  t_cell.Max();
  t_face.Max();

  Time t_cell_max = PPM->Max(t_cell.t);
  Time t_cell_min = PPM->Min(t_cell.t);
  Time t_face_max = PPM->Max(t_face.t);
  Time t_face_min = PPM->Min(t_face.t);

  mout.PrintInfo(Name(), verbose,
                 PrintInfoEntry<Time>("min cell assemble time", t_cell_min),
                 PrintInfoEntry<Time>("max cell assemble time", t_cell_max),
                 PrintInfoEntry<Time>("min face assemble time", t_face_min),
                 PrintInfoEntry<Time>("max face assemble time", t_face_max),
                 PrintInfoEntry<Time>("assemble time", Date() - Start),
                 PrintInfoEntry<Date>("finish assemble", Date()));
}

void STDGMatrixFreeViscoAcousticAssemble::SystemAddDoubleD(Matrix &M) const {
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(M, c, prob->nL());
    DGRowEntries M_c(M, *c, *c, false);
    const double rho = prob->Rho(*c, c());
//    vector<double> kappaTauInv(2 + prob->nL());
//    for (int j = 0; j < 2 + prob->nL(); j++) {
//      kappaTauInv[j] = 1.0 / (prob->Kappa(c(), j) * prob->Tau_i(c(), j));
//    }
    vector<double> kappaInv(2 + prob->nL());
    vector<double> kappaTauInv(2 + prob->nL());
    kappaInv[1] = 1.0 / prob->Kappa_i(*c, c(), 0);
    for (int j = 2; j < 2 + prob->nL(); j++) {
      kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j - 1);
      kappaTauInv[j] = kappaInv[j] / prob->Tau_i(*c, c(), j - 1);
    }
    for (int q = 0; q < elem.nQ(); ++q) {
      Point QP = elem.QPoint(q);
      double t = QP.t();
      double w = elem.QWeight(q);

      VectorField V_rhs = prob->F_Velocity(t, *c, QP);
      Scalar P_rhs = prob->F_Pressure(t, *c, QP);

      for (int i = 0; i < elem.i_dimension(); ++i) {
        int var_i = elem.variable(i);
        if (var_i < 2)
          continue;
        double P_i = elem.Pressure(q, i);
        for (int j = 0; j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          if (var_i != var_j)
            continue;
          Scalar P_j = elem.Pressure(q, j);
          M_c(i, j) += 2 * w * kappaTauInv[var_j] * (P_j * P_i);
        }
      }
    }
  }
}

pair<double, double> STDGMatrixFreeViscoAcousticAssemble::DiscNorm(
    const Matrix &L, const Vector &u) const {
  Vector Lu(L * u);
  Matrix M(u);
  MassMatrix(M);
  Vector Mu(M * u);

  Preconditioner *PC = GetPC("PointBlockGaussSeidel");
  LinearSolver *S = GetLinearSolver("GMRES", PC, "NormSolver");
  (*S)(M);

  Vector mLu(u);
  mLu = (*S) * Lu;

  double discNormW = u * Mu;
  double discNormV = Lu * mLu;

  // mout << "discrete norm ||u-U||_Wh  = " << scientific << sqrt(discNormW) <<
  // endl; mout << "discrete norm ||u-U||_Vh  = " << scientific <<
  // sqrt(discNormW + discNormV) << endl;

  // double theta = -1.0;
  // Config::Get("theta", theta);
  // Vector V(L.GetVector());
  // if (theta > 0) System(L, V);

  return pair<double, double>(std::sqrt(discNormW),
                              std::sqrt(discNormW + discNormV));
}

double STDGMatrixFreeViscoAcousticAssemble::L1Error(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField V = elem.Velocity(q, u);
      Scalar P = elem.Pressure(q, u);
      DampingVector PV = elem.DampingPressure(q, u);
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      V -= prob->ut_Velocity(QP.t(), QP, *c);
      P -= prob->ut_Pressure(QP.t(), QP, *c);
      PV -= prob->ut_DampingPressure(QP.t(), QP, *c);
      nrm += w * (absNorm(V) + abs(P) + absNorm(PV));
    }
  }
  return PPM->SumOnCommSplit(nrm, u.CommSplit());
}

void STDGMatrixFreeViscoAcousticAssemble::printConformReconstructionL2Error(
    Vector &u) const {
  if (u.dim() != 2 || prob->nL() != 0)
    Exit("check");

  double nrmW = 0;
  double nrmConformReconstr = 0;
  double nrmWConformReconstr = 0;

  double nrmConformReconstrAlternative = 0;
  double nrmWConformReconstrAlternative = 0;

  vector<vector<double>> ci(6);
  ci[0].resize(2);
  ci[0][0] = 0.0;
  ci[0][1] = 1.0;
  ci[1].resize(3);
  ci[1][0] = 0.0;
  ci[1][1] = 1.0 / 3.0;
  ci[1][2] = 1.0;
  ci[2].resize(4);
  ci[2][0] = 0.0;
  ci[2][1] = 0.4 - sqrt(0.06);
  ci[2][2] = 0.4 + sqrt(0.06);
  ci[2][3] = 1.0;
  ci[3].resize(5);
  ci[3][0] = 0.0;
  ci[3][1] = 0.088588;
  ci[3][2] = 0.409467;
  ci[3][3] = 0.7876595;
  ci[3][4] = 1.0;
  ci[4].resize(6);
  ci[4][0] = 0.0;
  ci[4][1] = 0.057104196114518;
  ci[4][2] = 0.276843013638124;
  ci[4][3] = 0.583590432368917;
  ci[4][4] = 0.860240135656219;
  ci[4][5] = 1.0;
  ci[5].resize(7);
  ci[5][0] = 0.0;
  ci[5][1] = 0.039809857051469;
  ci[5][2] = 0.198013417873608;
  ci[5][3] = 0.437974810247386;
  ci[5][4] = 0.695464273353636;
  ci[5][5] = 0.901464914201174;
  ci[5][6] = 1.0;

  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    row r = u.find_row(c());
    vector<Scalar> c_nodal_value;
    DegreePair deg = u.GetDoF().GetDegree(*c);
    if (c.min() == 0.0) {
      vector<Point> c_nodal = u.GetDoF().GetNodalPoints(*c);
      c_nodal_value.resize(u.GetDoF().get_m(c()));

      int cnt = 0;
      for (COMPONENT comp : prob->GetComponents()) {
        for (int i = 0; i < u.GetDoF().NodalPointsLocal(deg.space); ++i) {
          // Initial Value
          c_nodal_value[cnt] = prob->ut(0.0, c_nodal[cnt].CopyWithT(0.0), *c, comp);
          cnt++;
        }
      }
    }

    const Scalar *u_r = u(r);

    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField V = 0.0;
      Scalar P = 0.0;
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      for (int j = 0; j < elem.j_dimension(); ++j) {
        VectorField phi_V_j = elem.Velocity(q, j);
        double phi_P_j = elem.Pressure(q, j);
        P += u_r[j] * phi_P_j;
        V += u_r[j] * phi_V_j;
      }

      V -= prob->ut_Velocity(QP.t(), QP, *c);
      P -= prob->ut_Pressure(QP.t(), QP, *c);
      nrmW += w * (prob->Rho(*c, QP) * V * V + P * P / prob->Kappa(*c, QP));
    }

    cell c_prev = c;
    if (c.min() != 0) {
      face ff = u.find_face(c.Face(c.Faces() - 2));
      c_prev = u.find_cell(ff.Left());
      if (c_prev == u.cells_end())
        c_prev = u.find_overlap_cell(ff.Left());
    }
    row r_prev = u.find_row(c_prev());
    DegreePair deg_c_prev = u.GetDoF().GetDegree(*c_prev);
    const Scalar *u_r_prev = u(r_prev);

    const Quadrature TimeQ(GetQuadrature("Qint7"));
    double timedet = c.max() - c.min();
    const Quadrature SpaceQ(disc->GetCellQuad(deg.space));
    const Shape &spaceShape(disc->GetCellShape(deg.space));
    const Shape &timeShape(disc->GetTimeShape(deg.time));

    const int dofPerTdeg = r.n() / (deg.time + 1);

    int shift = deg_c_prev.time * r_prev.n() / (deg_c_prev.time + 1);

    for (int tq = 0; tq < TimeQ.size(); ++tq) {
      double time = c.min() + timedet * TimeQ.QPoint(tq)[0];
      double timeWeight = timedet * TimeQ.Weight(tq);
      for (int sq = 0; sq < SpaceQ.size(); ++sq) {
        Transformation T = c.SpaceCell().GetTransformation(SpaceQ.QPoint(sq));
        double space_qWeight = T.Det() * SpaceQ.Weight(sq);

        Point QP =
            c.SpaceCell().LocalToGlobal(SpaceQ.QPoint(sq)).CopyWithT(time);

        VectorField V = 0.0;
        Scalar P = 0.0;
        for (int pi = 0; pi < prob->Dim(); pi++) {
          COMPONENT comp = prob->GetComponents()[pi];
          for (int si = 0; si < spaceShape.size(); si++) {
            if (c.min() == 0) {
              vector<Point> c_nodal = u.GetDoF().GetNodalPoints(*c);
              for (int i = 0; i < c_nodal.size(); ++i)
                c_nodal[i] = c_nodal[i].CopyWithT(0.0);

              double tmp = prob->ut(c_nodal[si].t(), c_nodal[si], *c, comp) * spaceShape(sq, si);
              ;
              for (int tii = 1; tii <= deg.time + 1; tii++)
                tmp *= (1.0 - (TimeQ.QPoint(tq)[0] / ci[deg.time][tii]));
              if (pi < c.dim())
                V[pi] += tmp;
              else
                P += tmp;
            }
            if (c.min() != 0.0) {
              double tmp = u_r_prev[shift + pi * spaceShape.size() + si] *
                           spaceShape(sq, si);
              for (int tii = 1; tii <= deg.time + 1; tii++)
                tmp *= (1.0 - (TimeQ.QPoint(tq)[0] / ci[deg.time][tii]));
              if (pi < c.dim())
                V[pi] += tmp;
              else
                P += tmp;
            }
            for (int ti = 1; ti < deg.time + 1; ti++) {
              double tmp = 0.0;
              for (int tii = 1; tii <= deg.time + 1; tii++) {
                tmp +=
                    u_r[(tii - 1) * dofPerTdeg + pi * spaceShape.size() + si] *
                    spaceShape(sq, si) *
                    timeShape(Point(ci[deg.time][ti]), tii - 1);
              }
              for (int tii = 0; tii <= deg.time + 1; tii++) {
                if (tii == ti)
                  continue;
                tmp *= (TimeQ.QPoint(tq)[0] - ci[deg.time][tii]) /
                       (ci[deg.time][ti] - ci[deg.time][tii]);
              }
              if (pi < c.dim())
                V[pi] += tmp;
              else
                P += tmp;
            }
            for (int ti = deg.time + 1; ti <= deg.time + 1; ti++) {
              double tmp =
                  u_r[deg.time * dofPerTdeg + pi * spaceShape.size() + si] *
                  spaceShape(sq, si);
              for (int tii = 0; tii < deg.time + 1; tii++)
                tmp *= (TimeQ.QPoint(tq)[0] - ci[deg.time][tii]) /
                       (1.0 - ci[deg.time][tii]);
              if (pi < c.dim())
                V[pi] += tmp;
              else
                P += tmp;
            }
          }
        }

        V -= prob->ut_Velocity(QP.t(), QP, *c);
        P -= prob->ut_Pressure(QP.t(), QP, *c);

        double w = timeWeight * space_qWeight;

        nrmConformReconstr += w * (V * V + P * P);
        nrmWConformReconstr +=
            w * (V * V * prob->Rho(*c, QP) + P * P / prob->Kappa(*c, QP));
      }
    }

    // conform reconstruction equidistant
    for (int q = 0; q < elem.nQ(); ++q) {
      RVector X_conformReconstruction(0.0, prob->Dim());

      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      Point localPT;
      transformPointGlobalToLocal(QP, localPT, c);
      const Shape &timeShapeUp(disc->GetTimeShape(deg.time + 1));
      if (c.min() == 0) {
        vector<Point> c_nodal = u.GetDoF().GetNodalPoints(*c);

        int j = 0;
        int NX = 0;
        for (int p = 0; p < prob->Dim(); ++p) {
          COMPONENT comp = prob->GetComponents()[p];
          const Shape &shape(disc->GetCellShape(deg.space));
          for (int n = 0; n < u.GetDoF().NodalPointsLocal(deg.space, p); n++) {
            X_conformReconstruction[p] += shape(localPT, n) *
                                          timeShapeUp(Point(localPT.t()), 0) *
                                          prob->ut(c_nodal[j++].t(), c_nodal[j++], *c, comp);
          }
          NX += u.GetDoF().NodalPointsLocal(deg.space, p);
        }
      } else {
        // row r_prev = u.find_row(c_prev());
        // int prev_time_deg = u.GetDoF().get_time_deg(c_prev);
        // int shift = prev_time_deg * r_prev.n() / (prev_time_deg + 1);

        int NX = 0;
        // int prev_deg = u.GetDoF().get_space_deg(c_prev);

        for (int p = 0; p < prob->Dim(); ++p) {
          int nploc = u.GetDoF().NodalPointsLocal(deg_c_prev.space, p);
          const Shape &shape(disc->GetCellShape(deg_c_prev.space));
          for (int n = 0; n < nploc; n++) {
            X_conformReconstruction[p] += shape(localPT, n) *
                                          timeShapeUp(Point(localPT.t()), 0) *
                                          u(r_prev)[shift + NX + n];
          }
          NX += nploc;
        }
      }

      int NX = 0;
      // int deg = u.GetDoF().get_space_deg(c);

      for (int p = 0; p < prob->Dim(); ++p) {
        const int nploc = u.GetDoF().NodalPointsLocal(deg.space, p);
        // const Shape &shape(disc->GetCellShape(deg.space));
        // const Shape &timeShape(disc->GetTimeShape(deg.time));
        for (int n = 0; n < nploc; n++) {
          double SSvalue = spaceShape(localPT, n);
          for (int s = 0; s <= deg.time; ++s) {
            // double TSvalue = 0.0;
            for (int ss = 0; ss <= deg.time; ++ss) {
              int shift = ss * r.n() / (deg.time + 1);
              double ts = (s + 1) / (double(deg.time + 1));
              double TSvalue = timeShapeUp(Point(localPT.t()), s + 1) *
                               timeShape(Point(ts), ss);
              X_conformReconstruction[p] +=
                  SSvalue * TSvalue * u(r)[shift + NX + n];
            }
          }
        }
        NX += nploc;
      }

      VectorField V =
          VectorField{X_conformReconstruction[0], X_conformReconstruction[1]};
      V -= prob->ut_Velocity(QP.t(), QP, *c);
      Scalar P = X_conformReconstruction[c.dim()] - prob->ut_Pressure(QP.t(), QP, *c);
      if (prob->Dim() != 2 || prob->nL() != 0)
        Exit("check");
      nrmConformReconstrAlternative += w * (V * V + P * P);
      nrmWConformReconstrAlternative +=
          w * (prob->Rho(*c, QP) * V * V + P * P / prob->Kappa(*c, QP));
    }
  }
  vout(10) << "NormError in W                      = "
           << sqrt(PPM->SumOnCommSplit(nrmW, u.CommSplit())) << endl;
  vout(10) << "NormConfReconstrErrorRadau in L2    = "
           << sqrt(PPM->SumOnCommSplit(nrmConformReconstr, u.CommSplit()))
           << endl;
  vout(10) << "NormConfReconstrErrorRadau in W     = "
           << sqrt(PPM->SumOnCommSplit(nrmWConformReconstr, u.CommSplit()))
           << endl;
  vout(10) << "NormConfReconstrErrorEquidist in L2 = "
           << sqrt(PPM->SumOnCommSplit(nrmConformReconstrAlternative,
                                       u.CommSplit()))
           << endl;
  vout(10) << "NormConfReconstrErrorEquidist in W  = "
           << sqrt(PPM->SumOnCommSplit(nrmWConformReconstrAlternative,
                                       u.CommSplit()))
           << endl;
}

double STDGMatrixFreeViscoAcousticAssemble::DGError(const Vector &u) const {
  if (prob->nL() != 0)
    Exit("check");
  double nrm = 0;
  double h = u.GetMesh().MaxMeshWidth();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    double rho = prob->Rho(*c, c());
    double rho_inv = 1.0 / prob->Rho(*c, c());
    double kappa = prob->Kappa(*c, c());
    double kappaInv = 1.0 / kappa;
    double Z = sqrt(rho * kappa);
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    for (int q = 0; q < elem.nQ(); ++q) {
      Point QP = elem.QPoint(q);
      double t = QP.t();
      VectorField V = rho * elem.DtVelocity(q, u) - elem.GradPressure(q, u);
      double P = kappaInv * elem.DtPressure(q, u) - elem.DivVelocity(q, u);
      double w = elem.QWeight(q);
      V -= prob->F_Velocity(t, *c, QP);
      P -= prob->F_Pressure(t, *c, QP);
      nrm += h * w * (rho_inv * V * V + kappa * P * P);
    }
    for (int f = 0; f < c.Faces() - 2; ++f) {
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, f, prob->nL());
      if (u.OnBoundary(*c, f)) {
        int bnd_id = prob->BndID(c.Face(f));
        for (int q = 0; q < felem.nQ(); ++q) {
          VectorField Nq = felem.QNormal(q);
          Point QP = felem.QPoint(q);
          double w = felem.QWeight(q);
          if (bnd_id == 1) {
            double P = felem.Pressure(q, u);
            P -= prob->ut_Pressure(QP.t(), QP, *c);
            nrm += w * (1 / Z) * P * P;
          } else if (bnd_id == 2) {
            double Vn = felem.Velocity(q, u) * Nq;
            Vn -= prob->ut_Velocity(QP.t(), QP, *c) * Nq;
            nrm += w * Z * Vn * Vn;
          }
        }
        continue;
      }
      cell cf = u.find_neighbour_cell(c, f);
      if (c() < cf())
        continue;
      double Zf = sqrt(prob->Rho(*cf, cf()) * prob->Kappa(*cf, cf()));
      double Zinv = 1 / (Z + Zf);
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1,
                                                   prob->nL());
      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        VectorField Nq = felem.QNormal(q);
        Point QP = felem.QPoint(q);
        double w = felem.QWeight(q);
        double Vn = (felem.Velocity(q, u) - felem_1.Velocity(q1, u)) * Nq;
        double P = felem.Pressure(q, u) - felem_1.Pressure(q1, u);
        nrm += w * Zinv * (Z * Zf * Vn * Vn + P * P);
      }
    }
    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 2,
                                               prob->nL(), "dg");
    if (c.min() == u.t(0)) {
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        Point QP = felem.QPoint(q);
        VectorField V = felem.Velocity(q, u);
        double P = felem.Pressure(q, u);
        V -= prob->ut_Velocity(QP.t(), QP, *c);
        P -= prob->ut_Pressure(QP.t(), QP, *c);
        nrm += w * (rho * V * V + kappaInv * P * P);
      }
    } else if (u.has_previous_cell(c)) {
      cell c_prev = u.find_previous_cell(c);
      if (c_prev == u.cells_end() ||
          c_prev == u.overlap_end()) {
        std::stringstream s;
        s << c();
        Warning("A)) Cell has no previous cell!: " + s.str());
      } else if (c_prev == c) {
        std::stringstream s;
        s << c() << c_prev();
        Warning("B)) Cell has no previous cell!: " + s.str());
      }
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, c_prev,
                                                   c_prev.Faces() - 1,
                                                   prob->nL(), "dg", c);
      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        VectorField Nq = felem.QNormal(q);
        Point QP = felem.QPoint(q);
        double w = felem.QWeight(q);
        VectorField V = felem.Velocity(q, u) - felem_1.Velocity(q1, u);
        double P = felem.Pressure(q, u) - felem_1.Pressure(q1, u);
        nrm += w * (rho * V * V + kappaInv * P * P);
      }
    }
    if (c.max() == u.GetMesh().GetEndTime()) {
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 1,
                                                 prob->nL(), "dg");
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        Point QP = felem.QPoint(q);
        VectorField V = felem.Velocity(q, u);
        double P = felem.Pressure(q, u);
        V -= prob->ut_Velocity(QP.t(), QP, *c);
        P -= prob->ut_Pressure(QP.t(), QP, *c);
        nrm += w * (rho * V * V + kappaInv * P * P);
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGMatrixFreeViscoAcousticAssemble::L2Error(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    double rho = prob->Rho(*c, c());
    double kappaInv = 1.0 / prob->Kappa(*c, c());
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField V = elem.Velocity(q, u);
      Scalar P = elem.Pressure(q, u);
      DampingVector PV = elem.DampingPressure(q, u);
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      V -= prob->ut_Velocity(QP.t(), QP, *c);
      P -= prob->ut_Pressure(QP.t(), QP, *c);
      PV -= prob->ut_DampingPressure(QP.t(), QP, *c);
      const double pvScalarProduct = (PV * PV).sum();
      nrm += w * (rho * V * V + kappaInv * P * P + pvScalarProduct);
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGMatrixFreeViscoAcousticAssemble::LInfError(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField V = elem.Velocity(q, u);
      Scalar P = elem.Pressure(q, u);
      DampingVector DV = elem.DampingPressure(q, u);
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      V -= prob->ut_Velocity(QP.t(), QP, *c);
      P -= prob->ut_Pressure(QP.t(), QP, *c);
      DV -= prob->ut_DampingPressure(QP.t(), QP, *c);
      nrm = max(nrm, maxNorm(V));
      nrm = max(nrm, abs(P));
      nrm = max(nrm, maxNorm(DV));
    }
  }
  return PPM->Max(nrm, u.CommSplit());
};

double STDGMatrixFreeViscoAcousticAssemble::GNError(const Vector &) const {
  return -1.0;
}

double STDGMatrixFreeViscoAcousticAssemble::EnergyNorm(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    row r_c = u.find_row(c());
    const double *u_c = u(r_c);
    const double rho = prob->Rho(*c, c());
    vector<double> kappaInv(2 + prob->nL());
    kappaInv[1] = 1.0 / prob->Kappa_i(*c, c(), 0);
    for (int j = 2; j < 2 + prob->nL(); j++) {
      kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j - 1);
    }
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      for (int i = 0; i < elem.i_dimension(); ++i) {
        VectorField V = elem.Velocity(q, u);
        double P = elem.Pressure(q, u);
        nrm += w * (rho * V * V + kappaInv[0] * P * P);
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGMatrixFreeViscoAcousticAssemble::VNorm(const Matrix &L,
                                                  const Vector &u) const {
  Vector Lu = L * u;
  double normLu = L1Norm(Lu);
  double normSpaceIntegral = 0;
  double normSpaceTimeIntegral = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    row r_c = u.find_row(c());
    int c_deg = u.GetDoF().GetDegree(*c).space;
    const double rho = prob->Rho(*c, c());
    vector<double> kappaInv(2 + prob->nL());
    vector<double> kappaTauInv(2 + prob->nL());
    kappaInv[1] = 1.0 / prob->Kappa_i(*c, c(), 0);
    for (int j = 2; j < 2 + prob->nL(); j++) {
      kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j - 1);
      kappaTauInv[j] = kappaInv[j] / prob->Tau_i(*c, c(), j - 1);
    }
    const double z_K = sqrt(prob->Rho(*c, c()) * prob->Kappa(*c, c()));
    for (int f = 0; f < c.Faces() - 2; ++f) {
      double face_norm = 0;
      bool bnd = u.OnBoundary(*c, f);
      if (bnd)
        continue;  // TODO check this continue;
      int bnd_id = bnd ? prob->BndID(c.Face(f)) : -1;
      cell cf = bnd ? c : u.find_neighbour_cell(c, f);
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, f, prob->nL());
      int cf_deg = u.GetDoF().GetDegree(*cf).space;
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);
      row r_cf = u.find_row(cf());
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1,
                                                   prob->nL());
      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        double w = felem.QWeight(q);
        VectorField Nq = felem.QNormal(q);
        Point QP = felem.QPoint(q);

        double face_jump_p =
            felem.Pressure(q, u(r_c)) - felem_1.Pressure(q1, u(r_cf));
        double face_jump_v =
            Nq * (felem.Velocity(q, u(r_c)) - felem_1.Velocity(q1, u(r_cf)));
        normSpaceIntegral += w * (abs(face_jump_p) + abs(face_jump_v));
      }
    }

    cell c_prev = u.find_previous_cell(c);
    row r_cp = u.find_row(c_prev());
    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 2,
                                               prob->nL(), "dg");
    int c_prev_deg = u.GetDoF().GetDegree(*c_prev).space;
    SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, c_prev,
                                                 c_prev.Faces() - 1, prob->nL(),
                                                 "dg", c);

    if (c == c_prev) {
      continue;  // TODO check this continue
    }

    for (int q = 0; q < felem.nQ(); ++q) {
      double w = felem.QWeight(q);
      int q1 = find_q_id(felem.QPoint(q), felem_1);

      double face_jump_p =
          felem.Pressure(q, u(r_c)) - felem_1.Pressure(q1, u(r_cp));
      VectorField face_jump_v =
          felem.Velocity(q, u(r_c)) - felem_1.Velocity(q1, u(r_cp));
      double face_jump_v_value = abs(face_jump_v[0]) + abs(face_jump_v[1]);
      normSpaceIntegral +=
          w * (kappaInv[0] * abs(face_jump_p) + rho * face_jump_v_value);
    }
  }

  return PPM->SumOnCommSplit(normLu + normSpaceIntegral +
                                 normSpaceTimeIntegral / 2.0,
                             u.CommSplit());
}

double STDGMatrixFreeViscoAcousticAssemble::L1Norm(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField V = elem.Velocity(q, u);
      Scalar P = elem.Pressure(q, u);
      double w = elem.QWeight(q);
      nrm += w * (absNorm(V) + abs(P));
    }
  }
  return PPM->SumOnCommSplit(nrm, u.CommSplit());
}

double STDGMatrixFreeViscoAcousticAssemble::L2Norm(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    double rho = prob->Rho(*c, c());
    double kappaInv = 1.0 / prob->Kappa(*c, c());
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField V = elem.Velocity(q, u);
      Scalar P = elem.Pressure(q, u);
      double w = elem.QWeight(q);
      nrm += w * (rho * V * V + kappaInv * P * P);
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGMatrixFreeViscoAcousticAssemble::L2SpaceNormAtTime(
    const Vector &u, double time) const {
  double nrm = 0.0;
  if (time == u.GetMesh().GetEndTime()) {
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != time)
        continue;
      double rho = prob->Rho(*c, c());
      double kappaInv = 1.0 / prob->Kappa(*c, c());
      SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL(), true);
      double localtime = c.TimeCell().GlobalToLocal(time)[0];
      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        double w = elem.SpaceQWeight(q);
        Point evalP = elem.SpaceQPoint(q).CopyWithT(localtime);
        VectorField V = elem.VelocityLocal(evalP, u);
        double P = elem.PressureLocal(evalP, u);
        nrm += w * (kappaInv * P * P + rho * V * V);
      }
    }
  } else
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (!(c.min() == 0 && time == 0) || time <= c.min() || c.max() < time)
        continue;
      double rho = prob->Rho(*c, c());
      double kappaInv = 1.0 / prob->Kappa(*c, c());
      SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL(), true);
      double localtime = c.TimeCell().GlobalToLocal(time)[0];
      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        double w = elem.SpaceQWeight(q);
        Point evalP = elem.SpaceQPoint(q).CopyWithT(localtime);
        VectorField V = elem.VelocityLocal(evalP, u);
        double P = elem.PressureLocal(evalP, u);
        nrm += w * (kappaInv * P * P + rho * V * V);
      }
    }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGMatrixFreeViscoAcousticAssemble::L2SpaceNormAtTimeError(
    const Vector &u, double time) const {
  double nrm = 0.0;
  const double T = u.GetMesh().GetEndTime();
  if (time == T) {
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != T)
        continue;
      double rho = prob->Rho(*c, c());
      double kappaInv = 1.0 / prob->Kappa(*c, c());
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 1,
                                                 prob->nL(), "dg");
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        Point QP = felem.QPoint(q);
        VectorField V = felem.Velocity(q, u);
        double P = felem.Pressure(q, u);
        V -= prob->ut_Velocity(QP.t(), QP, *c);
        P -= prob->ut_Pressure(QP.t(), QP, *c);
        nrm += w * (rho * V * V + kappaInv * P * P);
      }
    }
    return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
  }
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    if (c.min() == time) {
      double rho = prob->Rho(*c, c());
      double kappaInv = 1.0 / prob->Kappa(*c, c());
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 2,
                                                 prob->nL(), "dg");
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        Point QP = felem.QPoint(q);
        VectorField V = felem.Velocity(q, u);
        double P = felem.Pressure(q, u);
        V -= prob->ut_Velocity(QP.t(), QP, *c);
        P -= prob->ut_Pressure(QP.t(), QP, *c);
        nrm += w * (rho * V * V + kappaInv * P * P);
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGMatrixFreeViscoAcousticAssemble::DGSemiNorm(const Vector &u) const {
  double nrm = 0.0;

  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    row r_c = u.find_row(c());
    int c_deg = u.GetDoF().GetDegree(*c).space;
    const double rho = prob->Rho(*c, c());
    vector<double> kappaInv(2 + prob->nL());
    vector<double> kappaTauInv(2 + prob->nL());
    kappaInv[1] = 1.0 / prob->Kappa_i(*c, c(), 0);
    for (int j = 2; j < 2 + prob->nL(); j++) {
      kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j - 1);
      kappaTauInv[j] = kappaInv[j] / prob->Tau_i(*c, c(), j - 1);
    }
    const double z_K = sqrt(prob->Rho(*c, c()) * prob->Kappa(*c, c()));
    for (int f = 0; f < c.Faces() - 2; ++f) {
      double face_norm = 0;
      bool bnd = u.OnBoundary(*c, f);
      if (bnd)
        continue;  // TODO check this continue;
      int bnd_id = bnd ? prob->BndID(c.Face(f)) : -1;
      cell cf = bnd ? c : u.find_neighbour_cell(c, f);
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, f, prob->nL());
      int cf_deg = u.GetDoF().GetDegree(*cf).space;
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);
      row r_cf = u.find_row(cf());
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1,
                                                   prob->nL());

      const double z_Kf = sqrt(prob->Rho(*cf, cf()) * prob->Kappa(*cf, cf()));
      const double alpha1 = 1.0 / (z_K + z_Kf);
      const double alpha2 = z_Kf * z_K * alpha1;
      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        double w = felem.QWeight(q);
        VectorField Nq = felem.QNormal(q);
        Point QP = felem.QPoint(q);

        double face_jump_p =
            felem.Pressure(q, u(r_c)) - felem_1.Pressure(q1, u(r_cf));
        double face_jump_v =
            Nq * (felem.Velocity(q, u(r_c)) - felem_1.Velocity(q1, u(r_cf)));
        face_jump_p = alpha1 * face_jump_p * face_jump_p;
        face_jump_v = alpha2 * face_jump_v * face_jump_v;
        nrm += w * (face_jump_p + face_jump_v);
      }
    }
    cell c_prev = u.find_previous_cell(c);
    row r_cp = u.find_row(c_prev());
    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 2,
                                               prob->nL(), "dg");
    int c_prev_deg = u.GetDoF().GetDegree(*c_prev).space;
    SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, c_prev,
                                                 c_prev.Faces() - 1, prob->nL(),
                                                 "dg", c);

    if (c == c_prev) {
      continue;  // TODO check this continue
    }

    for (int q = 0; q < felem.nQ(); ++q) {
      double w = felem.QWeight(q);
      int q1 = find_q_id(felem.QPoint(q), felem_1);

      double face_jump_p =
          felem.Pressure(q, u(r_c)) - felem_1.Pressure(q1, u(r_cp));
      VectorField face_jump_v =
          felem.Velocity(q, u(r_c)) - felem_1.Velocity(q1, u(r_cp));
      face_jump_p = kappaInv[0] * face_jump_p * face_jump_p;
      double face_jump_v_value = rho * face_jump_v * face_jump_v;
      nrm += w * (face_jump_p + face_jump_v_value);
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()) / 2);
};

double STDGMatrixFreeViscoAcousticAssemble::DGNormError(const Vector &u) {
  std::shared_ptr<AcousticProblem> temp_prob = prob;
  prob = CreateAcousticProblemShared("Zero");
  Matrix M_zero(u);
  Vector rhs_zero(0.0, u);
  System(M_zero, rhs_zero);
  prob = temp_prob;

  double normsq = 0;
  Vector Lu = M_zero * u;
  double buu = Lu * u;

  double h = u.GetMesh().MaxMeshWidth();

  for (cell c = u.cells(); c != u.cells_end(); c++) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());

    double rho = prob->Rho(*c, c());
    double rho_inv = 1.0 / prob->Rho(*c, c());
    double kappa = prob->Kappa(*c, c());
    double kappa_inv = 1.0 / prob->Kappa(*c, c());
    for (int q = 0; q < elem.nQ(); q++) {
      Point QP = elem.QPoint(q);
      double t = QP.t();
      double w = elem.QWeight(q);

      double dtP = elem.DtPressure(q, u);
      VectorField gradP = elem.GradPressure(q, u);
      VectorField dtV = elem.DtVelocity(q, u);
      double divV = elem.DivVelocity(q, u);

      double P = rho * dtP - divV - prob->F_Pressure(t, *c, QP);
      VectorField V = kappa_inv * dtV - gradP - prob->F_Velocity(t, *c, QP);

      normsq += w * (rho_inv * P * P + kappa * V * V);
    }
  }
  return sqrt(buu + PPM->SumOnCommSplit(h * normsq, u.CommSplit()));
}

double STDGMatrixFreeViscoAcousticAssemble::DGNorm(const Vector &u,
                                                   const Matrix &L) const {
  double normsq = 0;
  Vector Lu = L * u;
  double buu = Lu * u;

  for (cell c = u.cells(); c != u.cells_end(); c++) {
    double h = max(0.0, dist(c[0], c[c.Corners() - 1]));
    for (int cs = 0; cs < c.Corners() - 1; ++cs)
      h = max(h, dist(c[cs], c[cs + 1]));

    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());

    double rho = prob->Rho(*c, c());
    double rho_inv = 1.0 / prob->Rho(*c, c());
    double kappa = prob->Kappa(*c, c());
    double kappa_inv = 1.0 / prob->Kappa(*c, c());
    for (int q = 0; q < elem.nQ(); q++) {
      double w = elem.QWeight(q);

      double dtP = elem.DtPressure(q, u);
      VectorField gradP = elem.GradPressure(q, u);
      VectorField dtV = elem.DtVelocity(q, u);
      double divV = elem.DivVelocity(q, u);

      double P = rho * dtP - divV;
      VectorField V = kappa_inv * dtV - gradP;

      normsq += h * w * (rho_inv * P * P + kappa * V * V);
    }
  }
  return sqrt(buu + PPM->SumOnCommSplit(normsq, u.CommSplit()));
}

void STDGMatrixFreeViscoAcousticAssemble::DualRHS_linear(Vector &RHS,
                                                         const Vector &U,
                                                         Point a,
                                                         Point b) const {
  RHS = 0;

  if (a.t() != b.t()) {  // volume evaluation
    Exit("TODO")
  } else {  // evaluation at some time point
    for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {
      if (c.max() != a.t())
        continue;

      SpaceTimeViscoAcousticDGTElement elem(RHS, c, prob->nL());

      row r_RHS = RHS.find_row(c());

      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        double w = elem.SpaceQWeight(q);
        Point QP = elem.QPoint(q).CopyWithT(0.0);

        if ((a[0] > QP[0]) || (QP[0] > b[0]))
          continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1]))
          continue;

        VectorField j_V;
        for (int i = 0; i < c.dim(); ++i) {
          COMPONENT comp = prob->GetComponents()[i];
          j_V[i] = prob->ut_dual(QP, comp, a, b);
        }
        vector<Scalar> j_P(2 + prob->nL());
        for (int k = 0; k < 1 + prob->nL(); ++k) {
          COMPONENT comp = prob->GetComponents()[c.dim() + k];
          j_P[k + 1] = prob->ut_dual(QP, comp, a, b);
        }
        int time_deg = U.GetDoF().GetDegree(*c).time;

        int j_start = time_deg * (elem.j_dimension() / (time_deg + 1));
        Scalar scale = 0.0;

        for (int j = j_start; j < elem.j_dimension(); ++j)
          scale += elem.Velocity(q, j)[0];

        for (int j = j_start; j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          VectorField V_j = elem.Velocity(q, j) / scale;
          Scalar P_j = elem.Pressure(q, j) / scale;
          RHS(r_RHS)[j] -= w * (j_P[var_j] * P_j + j_V * V_j);
        }
      }
    }
  }
}

void STDGMatrixFreeViscoAcousticAssemble::DualRHS_quadratic(Vector &RHS,
                                                            const Vector &U,
                                                            Point a,
                                                            Point b) const {
  Exit("TODO") Exit("Not implemented");
}

double STDGMatrixFreeViscoAcousticAssemble::MminushalfLuMinusF(
    const cell &c, const Vector &u) const {
  double normsq = 0;
  SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
  double rho = prob->Rho(*c, c());
  double rho_inv = 1.0 / prob->Rho(*c, c());
  double kappa = prob->Kappa(*c, c());
  double kappa_inv = 1.0 / prob->Kappa(*c, c());
  for (int q = 0; q < elem.nQ(); q++) {
    Point QP = elem.QPoint(q);
    double t = QP.t();
    double w = elem.QWeight(q);

    double dtP = elem.DtPressure(q, u);
    VectorField gradP = elem.GradPressure(q, u);
    VectorField dtV = elem.DtVelocity(q, u);
    double divV = elem.DivVelocity(q, u);

    double P = rho * dtP - divV - prob->F_Pressure(t, *c, QP);
    VectorField V = kappa_inv * dtV - gradP - prob->F_Velocity(t, *c, QP);

    normsq += w * (rho_inv * P * P + kappa * V * V);
  }

  return sqrt(normsq);
}

double STDGMatrixFreeViscoAcousticAssemble::MhalfLu_exactMinusF(
    LevelPair levels) const {
  Vector res(disc, levels);
  res.Clear();
  for (cell c = res.cells(); c != res.cells_end(); c++) {
    double normsq = 0.0;
    SpaceTimeViscoAcousticDGTElement elem(res, c, prob->nL());
    double rho = prob->Rho(*c, c());
    double rho_inv = 1.0 / prob->Rho(*c, c());
    double kappa = prob->Kappa(*c, c());
    double kappa_inv = 1.0 / prob->Kappa(*c, c());
    for (int q = 0; q < elem.nQ(); q++) {
      Point QP = elem.QPoint(q);
      double t = QP.t();
      double w = elem.QWeight(q);

      double Mdtut_P = prob->Mdtut_Pressure(t, QP, *c);
      VectorField Mdtut_V = prob->Mdtut_Velocity(t, QP, *c);

      double Aut_P = prob->Aut_Pressure(t, QP, *c);
      VectorField Aut_V = prob->Aut_Velocity(t, QP, *c);

      double P = Mdtut_P + Aut_P - prob->F_Pressure(t, *c, QP);
      VectorField V = Mdtut_V + Aut_V - prob->F_Velocity(t, *c, QP);

      normsq += w * (rho_inv * P * P + kappa * V * V);
    }
    res(c(), 0) = sqrt(normsq);
  }
  // printVTK_Eta(v, *this, 100+U.Level().space);
  return res.norm();
}

double STDGMatrixFreeViscoAcousticAssemble::ResidualErrorEstimateJumpCell(
    const cell &c, const Vector &u) const {
  double eta = 0.0;
  double rho = prob->Rho(*c, c());
  double rho_inv = 1.0 / prob->Rho(*c, c());
  double kappa = prob->Kappa(*c, c());
  double kappaInv = 1.0 / kappa;
  double Z = sqrt(rho * kappa);
  SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());

  for (int q = 0; q < elem.nQ(); ++q) {
    VectorField V = rho * elem.DtVelocity(q, u) - elem.GradPressure(q, u);
    double P = kappaInv * elem.DtPressure(q, u) - elem.DivVelocity(q, u);
    double w = elem.QWeight(q);
    Point QP = elem.QPoint(q);
    double t = QP.t();
    V -= prob->F_Velocity(t, *c, QP);
    P -= prob->F_Pressure(t, *c, QP);
    eta += w * (rho_inv * V * V + kappa * P * P);
  }
  eta *= u.GetMesh().MaxMeshWidth();
  for (int f = 0; f < c.Faces() - 2; ++f) {
    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, f, prob->nL());
    if (u.OnBoundary(*c, f)) {
      int bnd_id = prob->BndID(c.Face(f));
      for (int q = 0; q < felem.nQ(); ++q) {
        VectorField Nq = felem.QNormal(q);
        Point QP = felem.QPoint(q);
        double w = felem.QWeight(q);
        if (bnd_id == 1) {
          double P = felem.Pressure(q, u);
          P -= prob->ut_Pressure(QP.t(), QP, *c);
          eta += w * (1 / Z) * P * P;
        } else if (bnd_id == 2) {
          double Vn = felem.Velocity(q, u) * Nq;
          Vn -= prob->ut_Velocity(QP.t(), QP, *c) * Nq;
          eta += w * Z * Vn * Vn;
        }
      }
      continue;
    }
    cell cf = u.find_neighbour_cell(c, f);
    int f1 = u.find_neighbour_face_id(c.Face(f), cf);
    SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1, prob->nL());
    double Zf = sqrt(prob->Rho(*cf, cf()) * prob->Kappa(*cf, cf()));
    double Zinv = 0.5 / (Z + Zf);
    for (int q = 0; q < felem.nQ(); ++q) {
      Point QP = felem.QPoint(q);
      int q1 = find_q_id(QP, felem_1);
      VectorField Nq = felem.QNormal(q);
      double w = felem.QWeight(q);
      double Vn = (felem.Velocity(q, u) - felem_1.Velocity(q1, u)) * Nq;
      double P = felem.Pressure(q, u) - felem_1.Pressure(q1, u);
      eta += w * Zinv * (Z * Zf * Vn * Vn + P * P);
    }
  }
  SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 2,
                                             prob->nL(), "dg");
  if (c.min() == u.t(0)) {
    for (int q = 0; q < felem.nQ(); ++q) {
      double w = felem.QWeight(q);
      Point QP = felem.QPoint(q);
      VectorField V = felem.Velocity(q, u);
      double P = felem.Pressure(q, u);
      V -= prob->ut_Velocity(QP.t(), QP, *c);
      P -= prob->ut_Pressure(QP.t(), QP, *c);
      eta += w * (rho * V * V + kappaInv * P * P);
    }
  } else {
    cell c_prev = u.find_previous_cell(c);
    if (c() != c_prev()) {
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, c_prev,
                                                   c_prev.Faces() - 1,
                                                   prob->nL(), "dg", c);
      for (int q = 0; q < felem.nQ(); ++q) {
        Point QP = felem.QPoint(q);
        double w = felem.QWeight(q);
        int q1 = find_q_id(QP, felem_1);
        VectorField Nq = felem.QNormal(q);
        VectorField V = felem.Velocity(q, u) - felem_1.Velocity(q1, u);
        double P = felem.Pressure(q, u) - felem_1.Pressure(q1, u);
        eta += 0.5 * w * (rho * V * V + kappaInv * P * P);
      }
    }
  }
  if (c.max() == u.GetMesh().GetEndTime())
    return sqrt(eta);
  face f = u.find_face(c.Face(c.Faces() - 1));
  if (!u.has_cell_or_overlap_cell(f.Right())) {
    return sqrt(eta);
  }
  cell c_next = u.find_cell_or_overlap_cell(f.Right());

  SpaceTimeViscoAcousticDGTFaceElement felem_0(*disc, u, c, c.Faces() - 1,
                                               prob->nL(), "dg");
  SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, c_next,
                                               c_next.Faces() - 2, prob->nL(),
                                               "dg", c);
  for (int q = 0; q < felem_0.nQ(); ++q) {
    int q1 = find_q_id(felem_0.QPoint(q), felem_1);
    VectorField Nq = felem_0.QNormal(q);
    Point QP = felem_0.QPoint(q);
    double w = felem_0.QWeight(q);
    VectorField V = felem_0.Velocity(q, u) - felem_1.Velocity(q1, u);
    double P = felem_0.Pressure(q, u) - felem_1.Pressure(q1, u);
    eta += 0.5 * w * (rho * V * V + kappaInv * P * P);
  }

  return sqrt(eta);
}

double
STDGMatrixFreeViscoAcousticAssemble::L2ReliableResidualErrorEstimateJumpCell(
    const cell &c, const Vector &u, const Vector &conf) const {
  SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());

  double Mhalf_u_minus_u_conf = 0.0;
  double Mhalfinv_Lu_minus_f = 0.0;

  double rho = prob->Rho(*c, c());
  double rho_inv = 1.0 / prob->Rho(*c, c());
  double kappa = prob->Kappa(*c, c());
  double kappa_inv = 1.0 / prob->Kappa(*c, c());

  for (int q = 0; q < elem.nQ(); q++) {
    Point QP = elem.QPoint(q);
    double t = QP.t();
    double w = elem.QWeight(q);

    double dtP = elem.DtPressure(q, conf);
    VectorField gradP = elem.GradPressure(q, conf);
    VectorField dtV = elem.DtVelocity(q, conf);
    double divV = elem.DivVelocity(q, conf);

    double P = rho * dtP - divV - prob->F_Pressure(t, *c, QP);
    VectorField V = kappa_inv * dtV - gradP - prob->F_Velocity(t, *c, QP);

    Mhalfinv_Lu_minus_f += w * (rho_inv * V * V + kappa * P * P);
    double u_minus_u_conf_P = elem.Pressure(q, u) - elem.Pressure(q, conf);
    VectorField u_minus_u_conf_V = elem.Velocity(q, u) - elem.Velocity(q, conf);
    Mhalf_u_minus_u_conf +=
        w * (rho * u_minus_u_conf_P * u_minus_u_conf_P +
             kappa_inv * u_minus_u_conf_V * u_minus_u_conf_V);
  }
  double normsq_0 = 0;
  if (c.min() == 0.0) {
    double Mhalf_u_conf0_minus_u0 = 0.0;
    for (int q = 0; q < elem.space_quad.size(); q++) {
      Point QP_local = elem.SpaceQPoint(q).CopyWithT(0.0);
      Point QP_global = c.LocalToGlobal(QP_local);
      double w = elem.SpaceQWeight(q);
      double u_conf0_minus_u0_P =
          elem.PressureLocal(QP_local, conf) -
                                  prob->ut_Pressure(QP_global.t(), QP_global, *c);
      VectorField u_conf0_minus_u0_V =
          elem.VelocityLocal(QP_local, conf) -
                                       prob->ut_Velocity(QP_global.t(), QP_global, *c);
      Mhalf_u_conf0_minus_u0 +=
          w * (rho * u_conf0_minus_u0_P * u_conf0_minus_u0_P +
               kappa_inv * u_conf0_minus_u0_V * u_conf0_minus_u0_V);
    }
    normsq_0 += Mhalf_u_conf0_minus_u0;
  }

  double T = u.GetMesh().GetEndTime();
  double normsq =
      Mhalf_u_minus_u_conf + 4 * T * T * (Mhalfinv_Lu_minus_f + normsq_0);

  return sqrt(normsq);
};

double STDGMatrixFreeViscoAcousticAssemble::DGErrorEstimateJumpCell(
    const cell &c, const Vector &u, const Vector &conf, double eta_res_c,
    double eta_conf_c) const {
  STDGDGViscoAcousticElement elem(u, c, prob->nL());
  double kappa = prob->Kappa(*c, c());
  double kappa_inv = 1.0 / prob->Kappa(*c, c());
  double rho = prob->Rho(*c, c());
  double rho_inv = 1.0 / prob->Rho(*c, c());

  double Mhalfinv_Lu_conf_minus_f = 0.0;
  for (int q = 0; q < elem.nQ(); q++) {
    Point QP = elem.QPoint(q);
    double t = QP.t();
    double w = elem.QWeight(q);

    double dtP = elem.DtPressure(q, conf);
    VectorField gradP = elem.GradPressure(q, conf);
    VectorField dtV = elem.DtVelocity(q, conf);
    double divV = elem.DivVelocity(q, conf);

    double P = rho * dtP - divV - prob->F_Pressure(t, *c, QP);
    VectorField V = kappa_inv * dtV - gradP - prob->F_Velocity(t, *c, QP);

    Mhalfinv_Lu_conf_minus_f += w * (rho_inv * P * P + kappa * V * V);
  }

  double Mhalf_uT_minus_u_confT = 0.0;
  if (c.min() == 0) {
    for (int q = 0; q < elem.space_quad.size(); q++) {
      Point QP_local = elem.SpaceQPoint(q).CopyWithT(0.0);
      Point QP_global = c.LocalToGlobal(QP_local);
      double w = elem.SpaceQWeight(q);
      double u_confT_minus_uT_P =
          elem.PressureLocal(QP_local, conf) -
                                  prob->ut_Pressure(QP_global.t(), QP_global, *c);
      VectorField u_confT_minus_uT_V =
          elem.VelocityLocal(QP_local, conf) -
                                       prob->ut_Velocity(QP_global.t(), QP_global, *c);
      Mhalf_uT_minus_u_confT +=
          w * (rho * u_confT_minus_uT_P * u_confT_minus_uT_P +
               kappa_inv * u_confT_minus_uT_V * u_confT_minus_uT_V);
    }
  } else if (c.max() == u.GetMesh().GetEndTime()) {
    for (int q = 0; q < elem.space_quad.size(); q++) {
      Point QP_local = elem.SpaceQPoint(q).CopyWithT(1.0);
      Point QP_global = c.LocalToGlobal(QP_local);
      double w = elem.SpaceQWeight(q);
      double u_confT_minus_uT_P =
          elem.PressureLocal(QP_local, conf) - elem.PressureLocal(QP_local, u);
      VectorField u_confT_minus_uT_V =
          elem.VelocityLocal(QP_local, conf) - elem.VelocityLocal(QP_local, u);
      Mhalf_uT_minus_u_confT +=
          w * (rho * u_confT_minus_uT_P * u_confT_minus_uT_P +
               kappa_inv * u_confT_minus_uT_V * u_confT_minus_uT_V);
    }
  }
  double dg_normsq = eta_res_c * eta_res_c;
  dg_normsq +=
      2 * eta_conf_c * sqrt(Mhalfinv_Lu_conf_minus_f) + Mhalf_uT_minus_u_confT;
  return sqrt(dg_normsq);
};

void STDGMatrixFreeViscoAcousticAssemble::ConformingLagrangeInterpolation(
    const Vector &u, Vector &conf, int degree) const {
  mout.StartBlock("ConfLagInterpolation");

  // pout << u << endl;
  // pout << u.GetMesh() << endl;

  auto meshes = MeshesCreator()
                    .WithMeshName(u.GetMesh().Name())
                    .WithPLevel(u.GetDisc().GetMeshes().pLevel())
                    .WithLevel(u.Level().space)
                    .CreateUnique();
  auto ld = std::make_shared<LagrangeDiscretization>(*meshes, degree, prob->Dim());
  // ld().procSets.CheckConsistency();
  // ld().GetProcSets().CheckConsistency();

  Vector l_vec(0.0, ld, u.Level());
  Vector adds(0.0, ld, u.Level());
  l_vec.Clear();
  adds.Clear();
  conf.Clear();

  for (cell c = u.cells(); c != u.cells_end(); c++) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    vector<Point> nps = l_vec.GetDoF().GetNodalPoints(*c);
    for (int component = 0; component < prob->Dim(); component++) {
      for (Point np : nps) {
        l_vec(np, component) += elem.EvaluateComponentGlobal(np, u, component);
        adds(np, component) += 1;
      }
    }
  }

  l_vec.Accumulate();
  adds.Accumulate();

  for (row r = l_vec.rows(); r != l_vec.rows_end(); r++) {
    for (int component = 0; component < prob->Dim(); component++) {
      l_vec(r, component) = l_vec(r, component) / adds(r, component);
    }
  }
  for (cell c = conf.cells(); c != conf.cells_end(); c++) {
    SpaceTimeViscoAcousticDGTElement elem(conf, c, prob->nL());
    ScalarElementT l_elem(l_vec, *c);
    vector<Point> nps = conf.GetDoF().GetNodalPoints(*c);
    DegreePair deg = conf.GetDoF().GetDegree(*c);
    int cnt = 0;
    for (int ti = 0; ti <= deg.time; ti++) {
      for (int pi = 0; pi < prob->Dim(); pi++) {
        for (int s = 0; s < conf.GetDoF().NodalPointsLocal(deg.space, pi);
             s++) {
          conf(elem.r, cnt) +=
              l_elem.Value(c.GlobalToLocal(nps[cnt]), l_vec, pi);
          cnt++;
        }
      }
    }
  }
  conf.Accumulate();
  mout.EndBlock(verbose > 0);
}

void STDGMatrixFreeViscoAcousticAssemble::ConformingInterpolation(
    const Vector &u, Vector &conf) const {
  std::unordered_map<Point, vector<double>> conf_data;
  std::unordered_map<Point, int> number_of_additions;

  for (cell c = u.cells(); c != u.cells_end(); c++) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    vector<Point> nps = u.GetDoF().GetNodalPoints(*c);
    for (int i = 0; i < elem.i_dimension(); ++i) {
      int component = elem.GetComponent(i);
      double value = elem.EvaluateComponentGlobal(nps[i], u, component);
      if (auto entry = conf_data.find(nps[i]); entry == conf_data.end()) {
        conf_data[nps[i]] = vector<double>(prob->Dim(), 0.0);
        number_of_additions[nps[i]] = 0;
      }
      conf_data[nps[i]][component] += value;
      number_of_additions[nps[i]]++;
    }
  }

  ExchangeBuffer buffer(u.CommSplit());
  for (int q = 0; q < PPM->Size(u.CommSplit()); q++) {
    if (q != PPM->Proc(u.CommSplit())) {
      for (auto &[point, components] : conf_data) {
        buffer.Send(q) << point << number_of_additions[point]
                       << components.size();
        for (double component : components) {
          buffer.Send(q) << component;
        }
      }
    }
  }
  buffer.Communicate();
  for (int q = 0; q < PPM->Size(u.CommSplit()); q++) {
    if (q != PPM->Proc(u.CommSplit())) {
      while (buffer.Receive(q).size() < buffer.ReceiveSize(q)) {
        Point point;
        vector<double> components;
        size_t comp_size;
        int adds;
        buffer.Receive(q) >> point >> adds >> comp_size;
        components.resize(comp_size);
        for (int i = 0; i < comp_size; i++) {
          buffer.Receive(q) >> components[i];
        }
        if (auto entry = conf_data.find(point); entry == conf_data.end()) {
          conf_data[point] = vector<double>(prob->Dim(), 0.0);
          number_of_additions[point] = 0;
        }
        for (int i = 0; i < conf_data[point].size(); i++) {
          conf_data[point][i] += components[i];
        }
        number_of_additions[point] += adds;
      }
    }
  }
  buffer.ClearBuffers();
  for (auto &[point, components] : conf_data) {
    for (double &component : components) {
      component *= prob->Dim() * 1.0 / number_of_additions[point];
    }
  }

  for (cell c = u.cells(); c != u.cells_end(); c++) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    vector<Point> nps = u.GetDoF().GetNodalPoints(*c);
    for (int i = 0; i < nps.size(); ++i) {
      conf(elem.r, i) = conf_data[nps[i]][elem.GetComponent(i)];
    }
  }
}

void STDGMatrixFreeViscoAcousticAssemble::ConformingProjectionWithPenalty(
    const Vector &u, Vector &u_conf, Vector &r_conf) const {
  if (penalty == 0)
    return;
  Matrix A_conf(u);
  A_conf = 0;
  r_conf = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    row r_c = u.find_row(c());
    // const double *u_c = u(r_c);

    DGRowEntries A_cc(A_conf, *c, *c, false);

    const double rho = prob->Rho(*c, c());
    vector<double> kappaInv(2 + prob->nL());
    kappaInv[0] = 1.0 / prob->Kappa_i(*c, c(), 0);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      VectorField V = elem.Velocity(q, u);
      double P = elem.Pressure(q, u);
      for (int i = 0; i < elem.i_dimension(); ++i) {
        int var_i = elem.variable(i);

        VectorField V_i = elem.Velocity(q, i);
        double P_i = elem.Pressure(q, i);
        r_conf(r_c, i) += w * (rho * (V * V_i) + kappaInv[0] * (P * P_i));
        for (int j = 0; j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          VectorField V_j = elem.Velocity(q, j);
          double P_j = elem.Pressure(q, j);
          A_cc(i, j) += w * ((rho * (V_j * V_i) + kappaInv[0] * (P_j * P_i)));
        }
      }
    }

    // penalty for the jumps in time
    cell c_prev = A_conf.find_previous_cell(c);
    if (c != c_prev) {
      SpaceTimeViscoAcousticDGTFaceElement felem_t(*disc, A_conf, c,
                                                   c.Faces() - 2, prob->nL(),
                                                   "dg");
      SpaceTimeViscoAcousticDGTFaceElement felem_t1(*disc, A_conf, c_prev,
                                                    c_prev.Faces() - 1,
                                                    prob->nL(), "dg", c);
      DGRowEntries A_cp(A_conf, *c, *c_prev, false);
      DGRowEntries A_pc(A_conf, *c_prev, *c, false);
      DGRowEntries A_pp(A_conf, *c_prev, *c_prev, false);
      for (int q = 0; q < felem_t.nQ(); ++q) {
        double w = felem_t.QWeight(q);
        if (felem_t.QPoint(q) != felem_t1.QPoint(q)) {
          Exit("nono");
        }

        for (int i = 0; i < felem_t.i_dimension(); ++i) {
          if (felem_t.is_zero_test(q, i))
            continue;
          int var_i = felem_t.variable(i);
          VectorField V_i = felem_t.Velocity(q, i);
          double P_i = felem_t.Pressure(q, i);

          for (int j = 0; j < felem_t.j_dimension(); ++j) {
            int var_j = felem_t.variable(j);
            if (felem_t.is_zero(q, j) || var_i != var_j)
              continue;

            VectorField V_j = felem_t.Velocity(q, j);
            double P_j = felem_t.Pressure(q, j);
            A_cc(i, j) +=
                penalty * w * (rho * (V_j * V_i) + kappaInv[0] * P_j * P_i);
          }
          for (int j = 0; j < felem_t1.j_dimension(); ++j) {
            int var_j = felem_t1.variable(j);
            if (felem_t1.is_zero(q, j) || var_i != var_j)
              continue;
            VectorField V1_j = felem_t1.Velocity(q, j);
            double P1_j = felem_t1.Pressure(q, j);
            double a =
                penalty * w * (rho * (V1_j * V_i) + kappaInv[0] * P1_j * P_i);
            A_cp(i, j) -= a;
            A_pc(j, i) -= a;
          }
        }
        for (int i = 0; i < felem_t1.i_dimension(); ++i) {
          if (felem_t1.is_zero_test(q, i))
            continue;
          int var_i = felem_t1.variable(i);
          VectorField V1_i = felem_t1.Velocity(q, i);
          double P1_i = felem_t1.Pressure(q, i);
          for (int j = 0; j < felem_t1.j_dimension(); ++j) {
            int var_j = felem_t1.variable(j);
            if (felem_t1.is_zero(q, j) || var_i != var_j)
              continue;
            VectorField V1_j = felem_t1.Velocity(q, j);
            double P1_j = felem_t1.Pressure(q, j);
            A_pp(i, j) +=
                penalty * w * (rho * (V1_j * V1_i) + kappaInv[0] * P1_j * P1_i);
          }
        }
      }
    }

    // penalty for  the jumps in space
    for (int f = 0; f < c.Faces() - 2; ++f) {
      if (u.OnBoundary(*c, f))
        continue;
      cell cf = u.find_neighbour_cell(c, f);
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, f, prob->nL());
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1,
                                                   prob->nL());
      DGRowEntries A_cf(A_conf, *c, *cf, false);
      DGRowEntries A_fc(A_conf, *cf, *c, false);
      DGRowEntries A_ff(A_conf, *cf, *cf, false);

      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        VectorField Nq = felem.QNormal(q);
        int q1 = find_q_id(felem.QPoint(q), felem_1);

        for (int i = 0; i < felem.i_dimension(); ++i) {
          int var_i = felem.variable(i);
          if (felem.is_zero_test(q, i))
            continue;
          double VN_i = felem.Velocity(q, i) * Nq;
          double P_i = felem.Pressure(q, i);
          for (int j = 0; j < felem.i_dimension(); ++j) {
            int var_j = felem.variable(j);
            if (felem.is_zero_test(q, j) || var_i != var_j)
              continue;
            double VN_j = felem.Velocity(q, j) * Nq;
            double P_j = felem.Pressure(q, j);
            A_cc(i, j) += penalty * w *
                          ((rho * (VN_j * VN_i) + kappaInv[0] * (P_j * P_i)));
          }
          for (int j = 0; j < felem_1.i_dimension(); ++j) {
            int var_j = felem_1.variable(j);
            if (felem_1.is_zero_test(q1, j) || var_i != var_j)
              continue;
            double VN_j = felem_1.Velocity(q1, j) * Nq;
            double P_j = felem_1.Pressure(q1, j);
            A_cf(i, j) -= penalty * w *
                          ((rho * (VN_j * VN_i) + kappaInv[0] * (P_j * P_i)));
          }
        }
      }
    }
  }
  r_conf -= A_conf * u_conf;

  mout << "Solving ConformingProjection with penality=" +
              std::to_string(penalty)
       << endl;
  LinearSolver *S = GetLinearSolverByPrefix("ConformingProjection");

  u_conf += (*S)(A_conf)*r_conf;
}

double STDGMatrixFreeViscoAcousticAssemble::DualErrorEstimateJumpCell(
    const cell &c, const Vector &u, const Vector &u_star) const {
  double eta = 0.0;

  SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
  row r = u.find_row(c());

  const Scalar *u_r = u(r);
  const Scalar *u_star_r = u_star(r);

  double h = dist(c[0], c[c.Corners() - 1]);
  for (int cs = 0; cs < c.Corners() - 1; ++cs)
    h = max(h, dist(c[cs], c[cs + 1]));
  double I = c.max() - c.min();

  double rho_1_V = 0.0;
  double rho_1_P = 0.0;
  double omega_1_V = 0.0;
  double omega_1_P = 0.0;

  DegreePair deg = u.GetDoF().GetDegree(*c);

  double rho = prob->Rho(*c, c());
  double kappaInv = 1.0 / prob->Kappa(*c, c());

  for (int q = 0; q < elem.nQ(); ++q) {
    VectorField V = 0.0;
    double P = 0.0;

    double w = elem.QWeight(q);
    Point QP = elem.QPoint(q);
    double t = QP.t();
    for (int j = 0; j < elem.j_dimension(); ++j) {
      VectorField phi_V_j =
          -1.0 * rho * elem.DtVelocity(q, j) + elem.GradPressure(q, j);
      double phi_P_j =
          -1.0 * kappaInv * elem.DtPressure(q, j) + elem.DivVelocity(q, j);
      V += u_r[j] * phi_V_j;
      P += u_r[j] * phi_P_j;
    }

    VectorField V_rhs = prob->F_Velocity(t, *c, QP);
    ;
    Scalar P_rhs = prob->F_Pressure(t, *c, QP);

    rho_1_V += w * (V * V + V_rhs * V_rhs);
    rho_1_P += w * (P * P + P_rhs * P_rhs);
  }

  vector<VectorField> dual_V_projection(deg.time + 1, zero);
  vector<double> dual_P_projection(deg.time + 1, 0.0);

  double K_vol = 0;
  for (int q = 0; q < elem.nQ(); ++q) {
    K_vol += elem.QWeight(q);
    for (int i = 0; i < elem.i_dimension(); ++i) {
      int l = i / (elem.i_dimension() / (deg.time + 1));
      dual_V_projection[l] +=
          elem.QWeight(q) * u_star_r[i] * elem.Velocity(q, i);
      dual_P_projection[l] +=
          elem.QWeight(q) * u_star_r[i] * elem.Pressure(q, i);
    }
  }
  for (int l = 0; l <= deg.time; ++l) {
    dual_V_projection[l] *= 1. / K_vol;
    dual_P_projection[l] *= 1. / K_vol;
  }

  const Shape &time_shape(disc->GetTimeShape(deg.time));

  vector<double> V_Fluxes(c.Faces(), 0.0);
  vector<double> P_Fluxes(c.Faces(), 0.0);
  vector<double> V_dual_Fluxes(c.Faces(), 0.0);
  vector<double> P_dual_Fluxes(c.Faces(), 0.0);

  for (int f = 0; f < c.Faces() - 2; ++f) {
    if (u.OnBoundary(*c, f)) {
      continue;
    }
    cell cf = u.find_neighbour_cell(c, f);
    row rf = u.find_row(cf());

    const Scalar *u_rf = u(rf);
    const Scalar *u_star_rf = u_star(rf);
    DegreePair deg_cf = u.GetDoF().GetDegree(*cf);
    int cf_deg = deg_cf.space;

    int f1 = 0;
    Point f_c = u.find_face(c.Face(f))();
    Point f_cf = Origin;
    for (; f1 < cf.Faces() - 2; ++f1) {
      f_cf = u.find_face(cf.Face(f1))();
      Point pfcf = f_cf.CopyWithT(0);
      Point pfc = f_c.CopyWithT(0);
      if (pfc == pfcf)
        break;
    }

    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, f, prob->nL());
    SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1, prob->nL());
    SpaceTimeViscoAcousticDGTElement elem_1(u, cf, prob->nL());

    double eta_v = 0.;
    double fxI = 0;

    for (int q = 0; q < felem.nQ(); ++q) {
      int q1 = find_q_id(felem.QPoint(q), felem_1);

      VectorField V_Flux = 0.0;
      double P_Flux = 0.0;

      double w = felem.QWeight(q);
      VectorField Nq = felem.QNormal(q);

      fxI += w;

      for (int j = 0; j < felem.j_dimension(); ++j) {
        V_Flux -= u_r[j] * felem.Velocity(q, j);
        P_Flux -= u_r[j] * felem.Pressure(q, j);
      }
      for (int j = 0; j < felem_1.j_dimension(); ++j) {
        V_Flux += u_rf[j] * felem_1.Velocity(q1, j);
        P_Flux += u_rf[j] * felem_1.Pressure(q1, j);
      }
      V_Fluxes[f] += w * V_Flux * V_Flux;
      P_Fluxes[f] += w * P_Flux * P_Flux;
    }

    int time_deg_cf = deg_cf.time;

    vector<VectorField> dual_V_projection_cf(time_deg_cf + 1);
    vector<double> dual_P_projection_cf(time_deg_cf + 1);
    for (int l = 0; l <= time_deg_cf; ++l) {
      dual_V_projection_cf[l] = 0.0;
      dual_P_projection_cf[l] = 0.0;
    }

    double Kf_vol = 0;

    for (int q = 0; q < elem_1.nQ(); ++q) {
      Kf_vol += elem_1.QWeight(q);
      for (int i = 0; i < elem_1.i_dimension(); ++i) {
        int l = i / (elem_1.i_dimension() / (time_deg_cf + 1));
        double w_u_star = elem_1.QWeight(q) * u_star_rf[i];
        dual_V_projection_cf[l] += w_u_star * elem_1.Velocity(q, i);
        dual_P_projection_cf[l] += w_u_star * elem_1.Pressure(q, i);
      }
    }
    for (int l = 0; l <= time_deg_cf; ++l) {
      dual_V_projection_cf[l] *= 1. / Kf_vol;
      dual_P_projection_cf[l] *= 1. / Kf_vol;
    }

    const Shape &time_shape_cf(disc->GetTimeShape(deg_cf.time));

    double f_vol = fxI / I;

    auto &timeQuad = disc->GetTimeQuad(deg_cf.time);

    for (int q = 0; q < timeQuad.size(); ++q) {
      Scalar values_dual_V = 0.0;
      Scalar values_dual_P = 0.0;
      for (int l = 0; l < time_shape.size(); ++l) {
        double timeValue = time_shape(timeQuad.QPoint(q), l);
        values_dual_V += dual_V_projection[l] * felem.QNormal(0) * timeValue;
        values_dual_P += dual_P_projection[l] * timeValue;
      }

      Scalar values_dual_V_cf = 0.0;
      Scalar values_dual_P_cf = 0.0;
      for (int l = 0; l < time_shape_cf.size(); ++l) {
        double timeValue = time_shape_cf(timeQuad.QPoint(q), l);
        values_dual_V_cf +=
            dual_V_projection_cf[l] * felem.QNormal(0) * timeValue;
        values_dual_P_cf += dual_P_projection_cf[l] * timeValue;
      }

      V_dual_Fluxes[f] += f_vol * pow(values_dual_V_cf - values_dual_V, 2) *
                          timeQuad.Weight(q) * I;
      P_dual_Fluxes[f] += f_vol * pow(values_dual_P_cf - values_dual_P, 2) *
                          timeQuad.Weight(q) * I;
    }
  }
  for (int f = c.Faces() - 2; f < c.Faces() - 1; ++f) {
    face ff = u.find_face(c.Face(c.Faces() - 2));
    cell c_prev = u.find_cell_or_overlap_cell(ff.Left());
    if (c_prev == u.overlap_end())
      continue;

    cell cf = c_prev;

    V_Fluxes[f] = 0;
    V_dual_Fluxes[f] = 0;
    P_Fluxes[f] = 0;
    P_dual_Fluxes[f] = 0;

    row rf = u.find_row(cf());

    const Scalar *u_rf = u(rf);
    const Scalar *u_star_rf = u_star(rf);

    int f1 = cf.Faces() - 2;
    Point f_c = u.find_face(c.Face(f))();
    Point f_cf = Origin;
    for (; f1 < cf.Faces(); ++f1) {
      f_cf = u.find_face(cf.Face(f1))();
      Point pfcf = f_cf.CopyWithT(0.0);
      Point pfc = f_c.CopyWithT(0.0);
      if (pfc == pfcf)
        break;
    }

    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, f, prob->nL(),
                                               "dgdg");
    SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1, prob->nL(),
                                                 "dgdg", c);
    SpaceTimeViscoAcousticDGTElement elem_1(u, cf, prob->nL());

    double eta_v = 0.;
    double fxI = 0;

    for (int q = 0; q < felem.nQ(); ++q) {
      int q1 = q;
      // if (f<c.Faces()-2) find_q_id(felem.QPoint(q),*felem_1);

      VectorField V_Flux = 0.0;
      double P_Flux = 0.0;

      double w = felem.QWeight(q);
      // VectorField Nq = felem.QNormal(q);

      fxI += w;

      for (int j = 0; j < felem.j_dimension(); ++j) {
        V_Flux -= u_r[j] * felem.Velocity(q, j);
        P_Flux -= u_r[j] * felem.Pressure(q, j);
      }
      for (int j = 0; j < felem_1.j_dimension(); ++j) {
        VectorField Vf_j = u_rf[j] * felem_1.Velocity(q1, j);
        double Pf_j = u_rf[j] * felem_1.Pressure(q1, j);
        V_Flux += Vf_j;
        P_Flux += Pf_j;
      }
      V_Fluxes[f] += w * V_Flux * V_Flux;
      P_Fluxes[f] += w * P_Flux * P_Flux;
    }

    VectorField dual_V_projection_cf(0.0);
    Scalar dual_P_projection_cf(0.0);
    double Kf_vol = 0;

    for (int q = 0; q < elem_1.nQ(); ++q) {
      Kf_vol += elem_1.QWeight(q);
      for (int i = 0; i < elem_1.i_dimension(); ++i) {
        dual_V_projection_cf +=
            elem_1.QWeight(q) * u_star_rf[i] * elem_1.Velocity(q, i);
        dual_P_projection_cf +=
            elem_1.QWeight(q) * u_star_rf[i] * elem_1.Pressure(q, i);
      }
    }
    dual_V_projection_cf *= 1. / Kf_vol;
    dual_P_projection_cf *= 1. / Kf_vol;

    double f_vol = fxI;

    dual_V_projection_cf -= dual_V_projection[0];
    dual_P_projection_cf -= dual_P_projection[0];

    V_dual_Fluxes[f] += f_vol * (dual_V_projection_cf * dual_V_projection_cf);
    P_dual_Fluxes[f] += f_vol * (dual_P_projection_cf * dual_P_projection_cf);
  }

  double eta_f = 0;
  for (int f = 0; f < c.Faces(); ++f) {
    omega_1_V += V_dual_Fluxes[f];
    omega_1_P += P_dual_Fluxes[f];
    eta_f += sqrt(V_Fluxes[f]) * sqrt(V_dual_Fluxes[f]) +
             sqrt(P_Fluxes[f]) * sqrt(P_dual_Fluxes[f]);
  }
  eta = sqrt(rho_1_P) * sqrt(h) * sqrt(omega_1_P) +
        sqrt(rho_1_V) * sqrt(h) * sqrt(omega_1_V) + 0.5 * eta_f;

  return eta;
}

double STDGMatrixFreeViscoAcousticAssemble::Goal_Functional(const Vector &u,
                                                            Point a,
                                                            Point b) const {
  Scalar e = 0.0;

  if (a.t() != b.t()) {  // volume evaluation
    Exit("ToDo");
  } else {  // evaluation at some time point
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != a.t())
        continue;
      SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
      row r = u.find_row(c());
      const int dim = c.dim();
      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        Point QP = elem.QPoint(q).CopyWithT(0.0);
        if ((a[0] > QP[0]) || (QP[0] > b[0]))
          continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1]))
          continue;
        double w = elem.SpaceQWeight(q);

        VectorField j_V;
        for (int i = 0; i < c.dim(); ++i) {
          COMPONENT comp = prob->GetComponents()[i];
          j_V[i] = prob->ut_dual(QP, comp, a, b);
        }
        vector<Scalar> j_P(2 + prob->nL());
        for (int k = 1; k <= 1 + prob->nL(); ++k) {
          COMPONENT comp = prob->GetComponents()[c.dim() - 1 + k];
          j_P[k] = prob->ut_dual(QP, COMPONENT::P0, a, b);
        }

        int time_deg = u.GetDoF().GetDegree(*c).time;
        const int start_idx = time_deg * (elem.j_dimension() / (time_deg + 1));
        Scalar scale = 0.0;
        for (int j = start_idx; j < elem.j_dimension(); ++j)
          scale += elem.Velocity(q, j)[0];

        VectorField V = 0.0;
        vector<Scalar> P(2 + prob->nL());
        for (int j = start_idx; j < elem.j_dimension(); ++j) {
          VectorField V_j = elem.Velocity(q, j) / scale;
          Scalar P_j = elem.Pressure(q, j) / scale;
          V += u(r, j) * V_j;
          P[elem.variable(j)] += u(r, j) * P_j;
        }
        e += w * (j_V * V);
        for (int k = 1; k <= 1 + prob->nL(); ++k)
          e += w * j_P[k] * P[k];
      }
    }
  }
  return abs(PPM->SumOnCommSplit(e, u.CommSplit()));
}

double STDGMatrixFreeViscoAcousticAssemble::Energy(const Vector &u, Point a,
                                                   Point b) const {
  Exit("TODO") Exit("Not implemented");
}

double STDGMatrixFreeViscoAcousticAssemble::L2ScalarProduct(
    const Vector &a, const Vector &b) const {
  if (a.GetDisc().DiscName() != b.GetDisc().DiscName()) {
    THROW("Discretization not the same in Scalarproduct!");
  }
  if (a.Level() != b.Level()) {
    THROW("Level not the same in Scalarproduct!")
  }

  double nrm = 0;

  for (cell c = a.cells(); c != a.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem_a(a, c, prob->nL());
    row r = a.find_row(c());
    DegreePair deg_a = a.GetDoF().GetDegree(*c);
    row rb = b.find_row(c());
    DegreePair deg_b = b.GetDoF().GetDegree(*c);
    if (deg_a.time != deg_b.time || deg_a.space != deg_b.space) {
      THROW("Degs not the same in Scalarproduct");
    }

    const Scalar *u_a = a(r);
    const Scalar *u_b = b(rb);

    for (int q = 0; q < elem_a.nQ(); ++q) {
      VectorField V_a = 0.0;
      Scalar P_a = 0.0;
      VectorField V_b = 0.0;
      Scalar P_b = 0.0;
      double w = elem_a.QWeight(q);
      Point QP = elem_a.QPoint(q);
      for (int j = 0; j < elem_a.j_dimension(); ++j) {
        VectorField phi_V_j = elem_a.Velocity(q, j);
        double phi_P_j = elem_a.Pressure(q, j);
        P_a += u_a[j] * phi_P_j;
        V_a += u_a[j] * phi_V_j;
        P_b += u_b[j] * phi_P_j;
        V_b += u_b[j] * phi_V_j;
      }
      nrm += w * (V_a * V_b + P_a * P_b);
    }
  }
  return PPM->SumOnCommSplit(nrm, a.CommSplit());
}

void STDGMatrixFreeViscoAcousticAssemble::get_projected_exact_solution(
    Vector &Exact_Solution) const {
  Exact_Solution.Clear();
  std::unordered_map<DegreePair, RMatrix> InvertedMassMatrixCache_t;
  std::unordered_map<DegreePair, SpaceTimeCellQuadrature> SpaceTimeQuadCache;
  for (cell c_t = Exact_Solution.cells(); c_t != Exact_Solution.cells_end();
       c_t++) {
    row row_t = Exact_Solution.find_row(c_t());
    DegreePair deg_t = Exact_Solution.GetDoF().GetDegree(*c_t);
    auto &space_shape_t = disc->GetCellShape(deg_t.space);
    auto &time_shape_t = disc->GetTimeShape(deg_t.time);
    int N_t = time_shape_t.size() * prob->Dim() * space_shape_t.size();

    auto elem_q = SpaceTimeQuadCache.find(deg_t);
    if (elem_q == SpaceTimeQuadCache.end()) {
      const Quadrature &s_quad = disc->GetCellQuad(deg_t.space);
      const Quadrature &t_quad = disc->GetTimeQuad(deg_t.time);
      elem_q =
          SpaceTimeQuadCache
              .insert({deg_t, SpaceTimeCellQuadrature(c_t, s_quad, t_quad)})
              .first;
    }
    SpaceTimeCellQuadrature &quad_t = elem_q->second;

    // SpaceTimeCellQuadrature quad_t(c_t, disc->GetCellQuad(deg_t.space),
    //                                disc->GetTimeQuad(deg_t.time));

    auto elem = InvertedMassMatrixCache_t.find(deg_t);

    if (elem == InvertedMassMatrixCache_t.end()) {
      RMatrix InvertedMassMatrix(N_t, N_t);
      for (int q = 0; q < quad_t.size(); q++) {
        Point local = quad_t.qPointLocal[q];
        for (int tc = 0; tc < time_shape_t.size(); ++tc) {
          double ts_t = time_shape_t(Point(local.t()), tc);
          for (int k = 0; k < prob->Dim(); ++k) {
            for (int sc = 0; sc < space_shape_t.size(); ++sc) {
              double ss_t = space_shape_t(local, sc);
              double v_t = ss_t * ts_t;
              int row = sc + (k + tc * prob->Dim()) * space_shape_t.size();
              for (int tc2 = 0; tc2 < time_shape_t.size(); ++tc2) {
                double ts_t2 = time_shape_t(Point(local.t()), tc2);
                for (int sc2 = 0; sc2 < space_shape_t.size(); ++sc2) {
                  int col =
                      sc2 + (k + tc2 * prob->Dim()) * space_shape_t.size();
                  double ss_t2 = space_shape_t(local, sc2);
                  double v_t2 = ss_t2 * ts_t2;
                  InvertedMassMatrix(row, col) +=
                      quad_t.qWeight[q] * v_t * v_t2;
                }
              }
            }
          }
        }
      }
      InvertedMassMatrix.Invert();
      elem =
          InvertedMassMatrixCache_t.insert({deg_t, InvertedMassMatrix}).first;
    }
    RMatrix InvertedMassMatrix_t = elem->second;

    // RMatrix MassMatrix_ft(N_t, N_t);
    RVector integral_t(N_t);
    for (int q = 0; q < quad_t.size(); q++) {
      Point local = quad_t.qPointLocal[q];
      Point global = c_t.LocalToGlobal(local);
      for (int tt = 0; tt < time_shape_t.size(); ++tt) {
        double value_time_shape = time_shape_t(Point(local.t()), tt);
        for (int st = 0; st < space_shape_t.size(); ++st) {
          double value_space_shape = space_shape_t(local, st);
          double v_t = value_space_shape * value_time_shape;
          // for (int tf = 0; tf < time_shape_t.size(); ++tf) {
          // double value_time_shape_f = time_shape_t(Point(local.t()), tf);
          // for (int sf = 0; sf < time_shape_t.size(); ++sf) {
          // double value_space_shape_f = time_shape_t(local, sf);
          for (int k = 0; k < prob->Dim(); ++k) {
            COMPONENT comp = prob->GetComponents()[k];
            int t_idx = st + (k + tt * prob->Dim()) * space_shape_t.size();
            // int f_idx = sf + (k + tf * prob->SpaceDim()) * space_shape_t.size();
            double v_f =
                prob->ut(global.t(), global, *c_t,
                         comp);  // value_space_shape_f * value_time_shape_f;
            // MassMatrix_ft[f_idx][t_idx] += quad_t.qWeight[q] * v_t * v_f;
            integral_t[t_idx] += quad_t.qWeight[q] * v_t * v_f;
          }
          //}
          //}
        }
      }
    }

    // RVector value_f(N_f);
    // for (int i = 0; i < value_f.size(); i++) {
    //   value_f[i] = from(row_f)[i];
    // }
    // mout << MassMatrix_ft << endl;
    // mout << InvertedMassMatrix_t << endl;
    // return;
    // RVector integral_t(MassMatrix_ft * value_f);
    RVector value_t = InvertedMassMatrix_t * integral_t;

    for (int i = 0; i < value_t.size(); i++) {
      Exact_Solution(row_t)[i] += value_t[i];
    }
  }

  /*
  Exact_Solution.Clear();
  const DoFs &dofs = Exact_Solution.GetDoF();
  for (cell c = Exact_Solution.cells(); c != Exact_Solution.cells_end(); ++c) {
    vector<Point> z;
    Exact_Solution.GetDoF().NodalPoints(c, z);
    row r = Exact_Solution.find_row(c());
    int space_deg = dofs.get_space_deg(c);
    int time_deg = dofs.get_time_deg(c);
    int cnt = 0;
    for (int ti = 0; ti <= time_deg; ti++) {
      for (int pi = 0; pi < prob->SpaceDim(); pi++) {
        for (int s = 0; s < dofs.NodalPointsLocal(space_deg, pi); s++) {
          Exact_Solution(r, cnt) = prob->ut(z[cnt], c, pi);
          cnt++;
        }
      }
    }
  }*/
  Exact_Solution.Accumulate();
  // mout << Exact_Solution << endl;
}


void STDGMatrixFreeViscoAcousticAssemble::PlotSingleSolution(const Vector &u,
                                                             const string &filename) const {
  auto vtuDisc = std::make_shared<STDiscretizationT_DGDG<>>(u.GetDisc().GetMeshes(), DegreePair{0, 0}, prob->Dim());
  Vector U_vtu(0.0, vtuDisc);

  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    row r = U_vtu.find_row(c());
    for (size_t i = 0; i < elem.GetComponentCount(); i++){
      U_vtu(r, i) = elem.EvaluateComponentGlobal(c(), u, i);
    }
  }

  VtuPlot plot{filename};
  vector<int> velocity_indices(SpaceDimension);
  std::iota(begin(velocity_indices), end(velocity_indices), 0);
  plot.AddData("Velocity", U_vtu, velocity_indices);
  plot.AddData("Pressure", U_vtu, SpaceDimension);
  for (size_t i = 0; i < prob->nL(); i++) {
    plot.AddData("PDamping"+std::to_string(i), U_vtu, SpaceDimension + i);
  }
  plot.PlotFile();
}
