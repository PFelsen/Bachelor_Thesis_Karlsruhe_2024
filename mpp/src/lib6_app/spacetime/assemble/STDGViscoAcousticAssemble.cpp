#include "STDGViscoAcousticAssemble.hpp"
#include "STDGDGViscoAcousticElement.hpp"
#include "STDGDGViscoAcousticFaceElement.hpp"
#include "DebuggingTools.hpp"

#include "LagrangeDiscretization.hpp"
#include "ScalarElement.hpp"

#include "SymRMatrix.hpp"
#include "LinearSolver.hpp"

using namespace std;

void STDGViscoAcousticAssemble::Energy(const Vector &u, vector<double> &Energy) const {
  vector<double> ts = u.GetMesh().GetTimesteps();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    if (c.min() == 0.0) {
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 2,
                                                 prob->nL(), "dgdg");
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        VectorField V = felem.Velocity(q, u);
        Energy[0] += w * prob->Rho(*c, c()) * V * V;
        for (int i = 0; i < felem.numL + 1; i++){
          double P = felem.Pressure(q, u, 1 + i);
          Energy[0] += w / prob->Kappa(*c, c()) * P * P;
        }
      }
    }

    int index = -1;
    for (int i = 0; i < ts.size(); i++) {
      if (ts[i] == c.max()) {
        index = i;
        break;
      }
    }
    if (index < 0) {
      THROW("Not found.");
    }

    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 1,
                                               prob->nL(), "dgdg");
    for (int q = 0; q < felem.nQ(); ++q) {
      double w = felem.QWeight(q);
      VectorField V = felem.Velocity(q, u);
      Energy[index] += w * prob->Rho(*c, c()) * V * V;
      for (int i = 0; i < felem.numL + 1; i++){
        double P = felem.Pressure(q, u, 1 + i);
        Energy[index] += w / prob->Kappa(*c, c()) * P * P; // TODO Kappa_i oder?
      }
    }
  }
  for (int i = 0; i < Energy.size(); i++){
    Energy[i] = sqrt(PPM->SumOnCommSplit(Energy[i], u.CommSplit()));
  }
}

void STDGViscoAcousticAssemble::MassMatrix(Matrix &M) const {
  M = 0.0;
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(M, c, prob->nL());
    DGRowEntries M_c(M, *c, *c, false);
    double rho = prob->Rho(*c, c());
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      for (int i = 0; i < elem.i_dimension(); ++i) {
        int var_i = elem.variable(i);
        VectorField V_i = elem.Velocity(q, i);
        double P_i = elem.Pressure(q, i);
        for (int j = 0; j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          if (var_j != var_i) continue;
          VectorField V_j = elem.Velocity(q, j);
          Scalar P_j = elem.Pressure(q, j);

          double kappaInv = 1.0 / prob->Kappa_i(*c, c(), var_j);
          M_c(i, j) += w * (rho * (V_j * V_i) + kappaInv * (P_j * P_i));
        }
      }
    }
  }
}

double STDGViscoAcousticAssemble::MhalfNorm(const Vector &u) const {
  double normsq = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); c++) {

    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());

    double rho = prob->Rho(*c, c() );
    double rho_inv = 1.0 / prob->Rho(*c, c() );
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

double STDGViscoAcousticAssemble::MhalfInvLNorm(const Vector &u) const {
  double normsq = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); c++) {

    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());

    double rho = prob->Rho(*c, c() );
    double rho_inv = 1.0 / prob->Rho(*c, c() );
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

void STDGViscoAcousticAssemble::System(Matrix &M, Vector &RHS) const {

  Date Start;
  mout.PrintInfo("TDGViscoAcousticAssemble_DGT", verbose,
                 PrintInfoEntry<int>("Problem Size", M.pSize(), 2),
                 PrintInfoEntry<Date>("start assemble", Start, 2)
  );
  Time t_cell;
  Time t_face;
  const double T = M.GetMesh().GetEndTime();
  std::function<double(Point)> d_T = [T](Point QP) { return (T - QP.t()); };
  std::function<double(Point)> constOne_func = [](Point QP) { return 1.0; };
  bool WeightedAssemble = false;
  Config::Get("WeightedAssemble", WeightedAssemble);
  std::function<double(Point)> weight_func = WeightedAssemble ? d_T : constOne_func;

  M = 0;
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    Date Start_cell;
    SpaceTimeViscoAcousticDGTElement elem(M, c, prob->nL());
    row r_c = M.find_row(c());
    DGRowEntries M_c(M, *c, *c, false);
    const double rho = prob->Rho(*c, c());
    Point globalPointSourceLocation = Infty;
    vector<Point> nodalPoints = RHS.GetNodalPoints(*c);
    if ( prob->hasPointSource(globalPointSourceLocation)) {
      if (globalPointSourceLocation == Infty){
        THROW("globalPointSourceLocation not set in Problem!")
      }
      if (c.SpaceCell().PointInCell(globalPointSourceLocation.GetSpatialPoint())) {
        for (int i = 0; i < elem.i_dimension(); ++i) {
          int var_i = elem.variable(i);
          if (var_i == 1) {
            for (int q = 0; q < elem.nTimeQ(); ++q) {
              double w = elem.TimeQWeight(q);
              Point localTimeQ = elem.TimeQPoint(q);
              Point globalTimeQ = c.TimeCell().LocalToGlobal(localTimeQ);
              Point temp = globalPointSourceLocation.CopyWithT(globalTimeQ[0]);
                RHS(c(), i) += w * elem.PressureGlobal(temp, i) *
                             prob->F_PointSource(temp.t(),temp);
            }
          }
        }
      }
    }

    vector<double> kappaInv(2 + prob->nL());
    vector<double> kappaTauInv(2 + prob->nL());
    kappaInv[1] = 1.0 / prob->Kappa_i(*c, c(), 0);
    for (int j = 2; j < 2 + prob->nL(); j++) {
      kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j - 1);
      kappaTauInv[j] = kappaInv[j] / prob->Tau_i(*c, c(), j - 1);
    }
    // Initialize elements
    for (int q = 0; q < elem.nQ(); ++q) {
      Point QP = elem.QPoint(q);
      double t = QP.t();
      double w = elem.QWeight(q) * weight_func(QP);
      VectorField V_rhs = prob->F_Velocity(t, *c, QP);
      Scalar P_rhs = prob->F_Pressure(t, *c, QP);
      DampingVector DP_rhs = prob->F_DampingPressure(t, *c, QP);
      // Assemble RHS

      for (int i = 0; i < elem.i_dimension(); ++i) {
        int var_i = elem.variable(i);
        Velocity V_i = elem.Velocity(q, i);
        double P_i = elem.Pressure(q, i);
        if (prob->HasRHS()) {
          if (var_i == 0){
            RHS(c(), i) += w * (V_rhs * V_i);
          } else if (var_i == 1) {
            RHS(c(), i) += w * (P_rhs * P_i);
          } else {
            RHS(c(), i) += w * (DP_rhs[var_i - 2] * P_i);
          }
        }
        //Assemble System matrix?
        for (int j = 0; j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          VectorField DtV_j = elem.DtVelocity(q, j);
          Scalar DivV_j = elem.DivVelocity(q, j);
          Scalar P_j = elem.Pressure(q, j);
          Scalar DtP_j = elem.DtPressure(q, j);
          VectorField GradP_j = elem.GradPressure(q, j);
          if (var_i != var_j) DtP_j = 0.0;
          M_c(i, j) += w * ((rho * (DtV_j * V_i)
                             + kappaInv[var_j] * (DtP_j * P_i))
                            - ((GradP_j * V_i) + (DivV_j * P_i)));
          if (var_j > 1 && var_i == var_j) {
            M_c(i, j) += w * kappaTauInv[var_j] * (P_j * P_i);
          }
        }
      }
    }

    // Initialize faces
    t_cell += Date() - Start_cell;
    Date Start_face;

    const double z_K = sqrt(prob->Rho(*c, c() ) * prob->Kappa(*c, c()));
    for (int f = 0; f < c.Faces() - 2; ++f) {
      bool bnd = M.OnBoundary(*c, f);
      int bnd_id = bnd ? prob->BndID(c.Face(f)) : -1;
      cell cf = bnd ? c : M.find_neighbour_cell(c, f);
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, M, c, f, prob->nL());
      int f1 = M.find_neighbour_face_id(c.Face(f), cf);
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, M, cf, f1, prob->nL());

      DGRowEntries M_cf(M, *c, *cf, false);
      const double z_Kf = sqrt(prob->Rho(*cf, cf()) * prob->Kappa(*cf, cf()));
      const double alpha1 = 1.0 / (z_K + z_Kf);
      const double alpha2 = z_Kf * z_K * alpha1;
      const double alpha3 = z_Kf * alpha1;
      const double alpha4 = z_K * alpha1;
      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        Point QP = felem.QPoint(q);
        double w = felem.QWeight(q) * weight_func(QP);
        VectorField Nq = felem.QNormal(q);

        Scalar VN_bnd = Nq * prob->ut_Velocity(QP.t(), QP, *c);
        Scalar P_bnd = prob->ut_Pressure(QP.t(), QP, *c) +
                       prob->ut_DampingPressure(QP.t(), QP, *c).sum();

        for (int i = 0; i < felem.i_dimension(); ++i) {
          if (felem.is_zero_test(q, i)) continue;
          VectorField V_i = felem.Velocity(q, i);
          double P_i = felem.Pressure(q, i);
          double VFlux_c_Vi = VFlux_c(V_i, Nq);
          if (bnd_id == 1) {
            RHS(c(), i) += w * (1 / z_K * P_i + VFlux_c_Vi) * P_bnd;
          } else if (bnd_id == 2 || bnd_id == 0) {
            RHS(c(), i) += w * (z_K * VFlux_c_Vi + P_i) * VN_bnd;
          }
          for (int j = 0; j < felem.j_dimension(); ++j) {
            if (felem.is_zero(q, j)) continue;
            VectorField V_j = felem.Velocity(q, j);
            double P_j = felem.Pressure(q, j);
            double VFlux_c_Vj = VFlux_c(V_j, Nq);
            M_c(i, j) += w * (alpha1 * P_j * P_i
                              + alpha2 * VFlux_c_Vj * VFlux_c_Vi
                              + alpha3 * VFlux_c_Vj * P_i
                              + alpha4 * P_j * VFlux_c_Vi);
          }
          for (int j = 0; j < felem_1.j_dimension(); ++j) {
            if (felem_1.is_zero(q1, j)) continue;
            VectorField V1_j = felem_1.Velocity(q1, j);
            double P1_j = felem_1.Pressure(q1, j);
            double PFlux_cf_P1j = PFlux_cf(P1_j, bnd, bnd_id);
            double VFlux_cf_V1j = VFlux_cf(V1_j, Nq, bnd, bnd_id);

            M_cf(i, j) -= w * (alpha1 * PFlux_cf_P1j * P_i
                               + alpha2 * VFlux_cf_V1j * VFlux_c_Vi
                               + alpha3 * VFlux_cf_V1j * P_i
                               + alpha4 * PFlux_cf_P1j * VFlux_c_Vi);
          }
        }
      }
    }

    cell c_prev = M.find_previous_cell(c);
    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, M, c, c.Faces() - 2, prob->nL(), "dg");
    SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, M, c_prev, c_prev.Faces() - 1,
                                                 prob->nL(), "dg", c);

    DGRowEntries M_c_prev(M, *c, *c_prev, false);
    for (int q = 0; q < felem.nQ(); ++q) {
      Point QP = felem.QPoint(q);
      double w = felem.QWeight(q) * weight_func(QP);
      for (int i = 0; i < felem.i_dimension(); ++i) {
        if (felem.is_zero_test(q, i)) continue;
        int var_i = felem.variable(i);
        VectorField V_i = felem.Velocity(q, i);
        double P_i = felem.Pressure(q, i);

        for (int j = 0; j < felem.j_dimension(); ++j) {
          int var_j = felem.variable(j);
          if (felem.is_zero(q, j) || var_i != var_j) continue;

          VectorField V_j = felem.Velocity(q, j);
          double P_j = felem.Pressure(q, j);
          M_c(i, j) += w * (rho * (V_j * V_i) + kappaInv[var_j] * P_j * P_i);
        }
        // Initial condition
        if (c.min() == 0.0) {
          VectorField V1_j = prob->ut_Velocity(QP.t(), QP, *c);
          double P1_j = prob->ut_Pressure(QP.t(), QP, *c);
          DampingVector DP1_j = prob->ut_DampingPressure(QP.t(), QP, *c);
          if (var_i == 0){
            RHS(c(), i) += w * (rho * (V1_j * V_i));
          } else if (var_i == 1) {
            RHS(c(), i) += w * (kappaInv[var_i] * P1_j * P_i);
          } else {
            RHS(c(), i) += w * (kappaInv[var_i] * DP1_j[var_i - 2] * P_i);
          }
        } else if (c != c_prev)
          for (int j = 0; j < felem_1.j_dimension(); ++j) {
            int var_j = felem_1.variable(j);
            if (felem_1.is_zero(q, j) || var_i != var_j) continue;
            VectorField V1_j = felem_1.Velocity(q, j);
            double P1_j = felem_1.Pressure(q, j);
            M_c_prev(i, j) -= w * (rho * (V1_j * V_i) + kappaInv[var_j] * P1_j * P_i);
          }
      }
    }
    t_face += Date() - Start_face;
  }

  Time t_cell_max = PPM->Max(t_cell.t);
  Time t_cell_min = PPM->Min(t_cell.t);
  Time t_face_max = PPM->Max(t_face.t);
  Time t_face_min = PPM->Min(t_face.t);

  mout.PrintInfo("TDGViscoAcousticAssemble_DGT", verbose,
                 PrintInfoEntry<Time>("min cell assemble time", t_cell_min),
                 PrintInfoEntry<Time>("max cell assemble time", t_cell_max),
                 PrintInfoEntry<Time>("min face assemble time", t_face_min),
                 PrintInfoEntry<Time>("max face assemble time", t_face_max),
                 PrintInfoEntry<Time>("assemble time", Date() - Start),
                 PrintInfoEntry<Date>("finish assemble", Date()));
}

void STDGViscoAcousticAssemble::RHS(Vector &RHS) const {

  const double T = RHS.GetMesh().GetEndTime();
  std::function<double(Point)> d_T = [T](Point QP) { return (T - QP.t()); };
  std::function<double(Point)> constOne_func = [](Point QP) { return 1.0; };
  bool WeightedAssemble = false;
  Config::Get("WeightedAssemble", WeightedAssemble);
  std::function<double(Point)> weight_func = WeightedAssemble ? d_T : constOne_func;

  const std::shared_ptr<SeismogramData> sourceData = prob->GetSourceData();
  if (sourceData) {
    const ObservationSpecification &spec = sourceData->GetParaDat().spec;
    //std::map<cell, vector<size_t>> cellToPointSource;
    for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {
      SpaceTimeViscoAcousticDGTElement elem(RHS, c, prob->nL());
      //const Shape &spaceShape = elem.shape.GetSpaceShape();
      for (size_t spaceLocationIndex = 0; spaceLocationIndex <
                                          sourceData->GetParaDat().LocalReceivers.size(); spaceLocationIndex++) {
        Point spaceLocation = sourceData->GetParaDat().LocalReceivers[spaceLocationIndex];
        if (c.SpaceCell().PointInCell(spaceLocation)) {
          for (int i = 0; i < elem.i_dimension(); ++i) {
            int var_i = elem.variable(i);
            if (var_i == 1) {
              for (int q = 0; q < elem.nTimeQ(); ++q) {
                double w = elem.TimeQWeight(q);
                Point localTimeQ = elem.LocalTimeQPoint(q);
                Point globalTimeQ = c.TimeCell().LocalToGlobal(localTimeQ);
                Point globalSTPoint = spaceLocation.CopyWithT(globalTimeQ[0]);
                int globalTimeIndex = (int) (globalTimeQ[0] / spec.GetTimeDelta());
                double globalTime = spec.TimeAtStep(globalTimeIndex);
                double sourceAmplitude = (*sourceData)(spaceLocationIndex, globalTimeIndex, 0);

                // int spaceIdx = i % (spaceShape.size() * elem.GetFullDimensions());
                // spaceIdx = spaceIdx % spaceShape.size();
                // spaceShape(spaceLocation, spaceIdx)
                RHS(c(), i) += w * elem.PressureGlobal(globalSTPoint, i) * sourceAmplitude;
              }
            }
          }
        }
      }
    }

  }

  for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(RHS, c, prob->nL());

    if (Point global_src; prob->HasRHS() && prob->hasPointSource(global_src)) {
      //Exit("TODO pkt quelle");
      if (c.PointInCell(global_src)) {
        for (int i = 0; i < elem.i_dimension(); ++i) {
          RHS(c(), i) += elem.PressureGlobal(global_src, i);
        }
      }
    }
    double rho = prob->Rho(*c, c());
    vector<double> kappaInv(2 + prob->nL());
    vector<double> kappaTauInv(2 + prob->nL());
    kappaInv[1] = 1.0 / prob->Kappa_i(*c, c(), 0);
    for (int j = 2; j < 2 + prob->nL(); j++) {
      kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j - 1);
      kappaTauInv[j] = kappaInv[j] / prob->Tau_i(*c, c(), j - 1);
    }


    for (int q = 0; q < elem.nQ(); ++q) {
      Point QP = elem.QPoint(q);
      double w = elem.QWeight(q) * weight_func(QP);
      VectorField V_rhs = prob->F_Velocity(QP.t(), *c, QP);
      Scalar P_rhs = prob->F_Pressure(QP.t(), *c, QP);
      DampingVector DP_rhs = prob->F_DampingPressure(QP.t(), *c, QP);

      for (int i = 0; i < elem.i_dimension(); ++i) {
        int var_i = elem.variable(i);
        Velocity V_i = elem.Velocity(q, i);
        double P_i = elem.Pressure(q, i);
        if (prob->HasRHS()) {
          if (var_i == 0) {
            RHS(c(), i) += w * (V_rhs * V_i);
          } else if (var_i == 1) {
            RHS(c(), i) += w * (P_rhs * P_i);
          } else {
            RHS(c(), i) += w * (DP_rhs[var_i - 2] * P_i);
          }
        }

      }
    }

    const double z_K = sqrt(prob->Rho(*c, c()) * prob->Kappa(*c, c()));
    for (int f = 0; f < c.Faces() - 2; ++f) {
      bool bnd = RHS.OnBoundary(*c, f);
      int bnd_id = bnd ? prob->BndID(c.Face(f)) : -1;
      cell cf = bnd ? c : RHS.find_neighbour_cell(c, f);
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, RHS, c, f, prob->nL());
      int f1 = RHS.find_neighbour_face_id(c.Face(f), cf);
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, RHS, cf, f1, prob->nL());

      const double z_Kf = sqrt(prob->Rho(*cf, cf()) * prob->Kappa(*cf, cf()));

      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        Point QP = felem.QPoint(q);
        double w = felem.QWeight(q) * weight_func(QP);
        VectorField Nq = felem.QNormal(q);

        Scalar VN_bnd = Nq * prob->ut_Velocity(QP.t(), QP, *c);
        Scalar P_bnd = prob->ut_Pressure(QP.t(), QP, *c) + prob->ut_DampingPressure(QP.t(), QP, *c).sum();

        for (int i = 0; i < felem.i_dimension(); ++i) {
          if (felem.is_zero_test(q, i)) continue;
          VectorField V_i = felem.Velocity(q, i);
          double P_i = felem.Pressure(q, i);
          double VFlux_c_Vi = VFlux_c(V_i, Nq);
          if (bnd_id == 1) {
            RHS(c(), i) += w * (1 / z_K * P_i + VFlux_c_Vi) * P_bnd;
          } else if (bnd_id == 2 || bnd_id == 0) {
            RHS(c(), i) += w * (z_K * VFlux_c_Vi + P_i) * VN_bnd;
          }
        }
      }
    }

    cell c_prev = RHS.find_previous_cell(c);
    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, RHS, c, c.Faces() - 2, prob->nL(), "dg");
    SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, RHS, c_prev, c_prev.Faces() - 1, prob->nL(),
                                                 "dg", c);

    for (int q = 0; q < felem.nQ(); ++q) {
      Point QP = felem.QPoint(q);
      double w = felem.QWeight(q) * weight_func(QP);
      for (int i = 0; i < felem.i_dimension(); ++i) {
        if (felem.is_zero_test(q, i)) continue;
        int var_i = felem.variable(i);
        VectorField V_i = felem.Velocity(q, i);
        double P_i = felem.Pressure(q, i);

        if (c.min() == 0.0) {
          VectorField V1_j = prob->ut_Velocity(QP.t(), QP, *c);
          double P1_j = prob->ut_Pressure(QP.t(), QP, *c);
          DampingVector DP1_j = prob->ut_DampingPressure(QP.t(), QP, *c);
          if (var_i == 0) {
            RHS(c(), i) += w * (rho * (V1_j * V_i));
          } else if (var_i == 1) {
            RHS(c(), i) += w * (kappaInv[var_i] * P1_j * P_i);
          } else {
            RHS(c(), i) += w * (kappaInv[var_i] * DP1_j[var_i - 2] * P_i);
          }
        }
      }
    }
  }
}


SeismogramData STDGViscoAcousticAssemble::measure(const ObservationSpecification &specification,
                                                  const Vector &u) const {
  if (u.GetMesh().GetTimesteps().size() > specification.GetNumberOfTimeSteps()) {
    THROW("Need more measure timesteps for this mesh." \
           " MeasureTimeSteps: " + std::to_string(specification.GetNumberOfTimeSteps())
          + " MeshTimesteps: " + std::to_string(u.GetMesh().GetTimesteps().size()));
  }

  std::map<cell, vector<int>> cellsToMeasureOnLocalReceiver;
  SeismogramParallelizationData pData(specification);
  for (int i = 0; i < specification.size(); i++) {
    int localIndex = -1;
    for (cell c = u.GetMesh().cells(); c != u.GetMesh().cells_end(); c++) {
      if (c.SpaceCell().PointInCell(specification[i])) {
        if (localIndex < 0) {
          localIndex = pData.LocalReceivers.size();
          pData.LocalReceivers.push_back(specification[i]);

          pData.CoordsToLocalIndex.insert({specification[i], localIndex});
          pData.LocalIndexToGlobalIndex.push_back(i);
          pData.GlobalTimesOwned.emplace_back(std::vector<int>{std::numeric_limits<int>::max(), 0});
        }
        cellsToMeasureOnLocalReceiver[c].push_back(localIndex);
        double multiple = c.min() - std::round(c.min() / specification.GetTimeDelta()) *
                                    specification.GetTimeDelta();
        if (abs(multiple) > 1e-9) {Warning("measuregird in time does not align with cells"+ to_string(multiple))}
        multiple = c.max() - std::round(c.max() / specification.GetTimeDelta()) *
                             specification.GetTimeDelta();
        if (abs(multiple) > 1e-9) {Warning("measuregird in time does not align with cells")}

        int indexMin = (int) std::round(c.min() / specification.GetTimeDelta());
        int indexMax = (int) std::round(c.max() / specification.GetTimeDelta());

        vector<int> oldTimesOwned = pData.GlobalTimesOwned[localIndex];
        pData.GlobalTimesOwned[localIndex] = std::vector<int>{
            std::min(oldTimesOwned[0], indexMin),
            std::max(oldTimesOwned[1], indexMax)};


      }
    }
    if (localIndex >= 0) {
      auto &globTimes = pData.GlobalTimesOwned[localIndex];
      int sizeOfRange = (globTimes[1] - globTimes[0] + 1) * specification.GetNrOfMeasure();

      int start = (localIndex == 0) ? 0 : pData.LocalIndexToLocalRange[localIndex - 1][1];
      pData.LocalIndexToLocalRange.push_back({start, start + sizeOfRange});
      int iTimesSteps = i * specification.GetNumberOfTimeSteps();
      pData.LocalIndexToGlobalRange.push_back(
          {(iTimesSteps + globTimes[0]) * specification.GetNrOfMeasure(),
           (iTimesSteps + globTimes[1]) * specification.GetNrOfMeasure()});
    }
  }

  SeismogramData data(pData);

  std::unordered_map<int, std::map<int, int>> localReceiverToGlobalTimeIndexCount;
  for (auto &[c, localReceiverIDs]: cellsToMeasureOnLocalReceiver) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    int globalTimeIndex = (int)std::round(c.min() / specification.GetTimeDelta());
    int measureCount = (int) std::round((c.max() - c.min()) / specification.GetTimeDelta()) + 1;
    for (int i = 0;  i < measureCount; i++, globalTimeIndex++) {
      double globalTime = c.min() + i * specification.GetTimeDelta();
      for (int localReceiverID: localReceiverIDs) {
        Point measurePkt = pData.LocalReceivers[localReceiverID].CopyWithT(globalTime);
        double pressure = elem.PressureGlobal(measurePkt, u);

        data(localReceiverID, globalTimeIndex, 0) += pressure;
        localReceiverToGlobalTimeIndexCount[localReceiverID][globalTimeIndex]++;
      }
    }
  }
  for (auto &[localReceiverID, globalTimesIndexCount]: localReceiverToGlobalTimeIndexCount) {
    for (auto &[globalTimeIndex, count]: globalTimesIndexCount) {
      data(localReceiverID, globalTimeIndex, 0) /= count;
    }
  }


  return data;
}

void STDGViscoAcousticAssemble::SystemAddDoubleD(Matrix &M) const {
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(M, c, prob->nL());
    DGRowEntries M_c(M, *c, *c, false);
    const double rho = prob->Rho(*c, c() );
    vector<double> kappaTauInv(2 + prob->nL());
    kappaTauInv[0] = kappaTauInv[1] = 0.0;
    for (int j = 2; j < 2 + prob->nL(); j++) {
      kappaTauInv[j] = 1.0 / (prob->Kappa_i(*c, c(), j - 1) * prob->Tau_i(*c, c(), j - 1));
    }

    for (int q = 0; q < elem.nQ(); ++q) {
      Point QP = elem.QPoint(q);
      double w = elem.QWeight(q);

      VectorField V_rhs = prob->F_Velocity(QP.t(), *c, QP);
      Scalar P_rhs = prob->F_Pressure(QP.t(), *c, QP);

      for (int i = 0; i < elem.i_dimension(); ++i) {
        int var_i = elem.variable(i);
        if (var_i < 2) continue;
        double P_i = elem.Pressure(q, i);
        for (int j = 0; j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          if (var_i != var_j) continue;
          Scalar P_j = elem.Pressure(q, j);
          M_c(i, j) += 2 * w * kappaTauInv[var_j] * (P_j * P_i);
        }
      }
    }
  }
}

pair<double, double> STDGViscoAcousticAssemble::DiscNorm(const Matrix &L, const Vector &u) const {
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

  //mout << "discrete norm ||u-U||_Wh  = " << scientific << sqrt(discNormW) << endl;
  //mout << "discrete norm ||u-U||_Vh  = " << scientific << sqrt(discNormW + discNormV) << endl;

  //double theta = -1.0;
  //Config::Get("theta", theta);
  //Vector V(L.GetVector());
  //if (theta > 0) System(L, V);

  return pair<double, double>(std::sqrt(discNormW), std::sqrt(discNormW + discNormV));
}

double STDGViscoAcousticAssemble::L1Error(const Vector &u) const {
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

void STDGViscoAcousticAssemble::printConformReconstructionL2Error(Vector &u) const {
  if (u.dim() != 2 || prob->nL() != 0) Exit("check");

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

        Point QP = c.SpaceCell().LocalToGlobal(SpaceQ.QPoint(sq)).CopyWithT(time);

        VectorField V = 0.0;
        Scalar P = 0.0;
        for (int pi = 0; pi < prob->Dim(); pi++) {
          COMPONENT comp = prob->GetComponents()[pi];
          for (int si = 0; si < spaceShape.size(); si++) {
            if (c.min() == 0) {
              vector<Point> c_nodal = u.GetDoF().GetNodalPoints(*c);
              for (int i = 0; i < c_nodal.size(); ++i)
                c_nodal[i] = c_nodal[i].CopyWithT(0.0);

              double tmp = prob->ut(c_nodal[si].t(), c_nodal[si], *c, comp) * spaceShape(sq, si);;
              for (int tii = 1; tii <= deg.time + 1; tii++)
                tmp *= (1.0 - (TimeQ.QPoint(tq)[0] / ci[deg.time][tii]));
              if (pi < c.dim()) V[pi] += tmp;
              else P += tmp;
            }
            if (c.min() != 0.0) {
              double tmp = u_r_prev[shift + pi * spaceShape.size() + si] * spaceShape(sq, si);
              for (int tii = 1; tii <= deg.time + 1; tii++)
                tmp *= (1.0 - (TimeQ.QPoint(tq)[0] / ci[deg.time][tii]));
              if (pi < c.dim()) V[pi] += tmp;
              else P += tmp;
            }
            for (int ti = 1; ti < deg.time + 1; ti++) {
              double tmp = 0.0;
              for (int tii = 1; tii <= deg.time + 1; tii++) {
                tmp += u_r[(tii - 1) * dofPerTdeg + pi * spaceShape.size() + si]
                       * spaceShape(sq, si)
                       * timeShape(Point(ci[deg.time][ti]), tii - 1);
              }
              for (int tii = 0; tii <= deg.time + 1; tii++) {
                if (tii == ti) continue;
                tmp *= (TimeQ.QPoint(tq)[0] - ci[deg.time][tii])
                       / (ci[deg.time][ti] - ci[deg.time][tii]);
              }
              if (pi < c.dim()) V[pi] += tmp;
              else P += tmp;
            }
            for (int ti = deg.time + 1; ti <= deg.time + 1; ti++) {
              double tmp = u_r[deg.time * dofPerTdeg + pi * spaceShape.size() + si]
                           * spaceShape(sq, si);
              for (int tii = 0; tii < deg.time + 1; tii++)
                tmp *= (TimeQ.QPoint(tq)[0] - ci[deg.time][tii]) / (1.0 - ci[deg.time][tii]);
              if (pi < c.dim()) V[pi] += tmp;
              else P += tmp;
            }
          }
        }


        V -= prob->ut_Velocity(QP.t(), QP, *c);
        P -= prob->ut_Pressure(QP.t(), QP, *c);

        double w = timeWeight * space_qWeight;

        nrmConformReconstr += w * (V * V + P * P);
        nrmWConformReconstr += w * (V * V * prob->Rho(*c, QP) + P * P / prob->Kappa(*c, QP));
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
            X_conformReconstruction[p] += shape(localPT, n)
                                          * timeShapeUp(Point(localPT.t()), 0)
                                          * prob->ut(c_nodal[j++].t(), c_nodal[j++], *c, comp);
          }
          NX += u.GetDoF().NodalPointsLocal(deg.space, p);
        }
      } else {
        //row r_prev = u.find_row(c_prev());
        //int prev_time_deg = u.GetDoF().get_time_deg(c_prev);
        //int shift = prev_time_deg * r_prev.n() / (prev_time_deg + 1);

        int NX = 0;
        //int prev_deg = u.GetDoF().get_space_deg(c_prev);

        for (int p = 0; p < prob->Dim(); ++p) {
          int nploc = u.GetDoF().NodalPointsLocal(deg_c_prev.space, p);
          const Shape &shape(disc->GetCellShape(deg_c_prev.space));
          for (int n = 0; n < nploc; n++) {
            X_conformReconstruction[p] += shape(localPT, n)
                                          * timeShapeUp(Point(localPT.t()), 0)
                                          * u(r_prev)[shift + NX + n];
          }
          NX += nploc;
        }
      }

      int NX = 0;
      //int deg = u.GetDoF().get_space_deg(c);

      for (int p = 0; p < prob->Dim(); ++p) {
        const int nploc = u.GetDoF().NodalPointsLocal(deg.space, p);
        //const Shape &shape(disc->GetCellShape(deg.space));
        //const Shape &timeShape(disc->GetTimeShape(deg.time));
        for (int n = 0; n < nploc; n++) {
          double SSvalue = spaceShape(localPT, n);
          for (int s = 0; s <= deg.time; ++s) {
            //double TSvalue = 0.0;
            for (int ss = 0; ss <= deg.time; ++ss) {
              int shift = ss * r.n() / (deg.time + 1);
              double ts = (s + 1) / (double(deg.time + 1));
              double TSvalue = timeShapeUp(Point(localPT.t()), s + 1) * timeShape(Point(ts), ss);
              X_conformReconstruction[p] += SSvalue * TSvalue * u(r)[shift + NX + n];
            }
          }
        }
        NX += nploc;
      }

      VectorField V = VectorField{X_conformReconstruction[0], X_conformReconstruction[1]};
      V -= prob->ut_Velocity(QP.t(), QP, *c);
      Scalar P = X_conformReconstruction[c.dim()] - prob->ut_Pressure(QP.t(), QP, *c);
      if (prob->Dim() != 2 || prob->nL() != 0) Exit("check");
      nrmConformReconstrAlternative += w * (V * V + P * P);
      nrmWConformReconstrAlternative += w * (prob->Rho(*c, QP) * V * V + P * P / prob->Kappa(*c, QP));
    }
  }
  vout(10) << "NormError in W                      = "
           << sqrt(PPM->SumOnCommSplit(nrmW, u.CommSplit())) << endl;
  vout(10) << "NormConfReconstrErrorRadau in L2    = "
           << sqrt(PPM->SumOnCommSplit(nrmConformReconstr, u.CommSplit()))
           << endl;
  vout(10) << "NormConfReconstrErrorRadau in W     = "
           << sqrt(PPM->SumOnCommSplit(nrmWConformReconstr, u.CommSplit())) << endl;
  vout(10) << "NormConfReconstrErrorEquidist in L2 = "
           << sqrt(PPM->SumOnCommSplit(nrmConformReconstrAlternative, u.CommSplit())) << endl;
  vout(10) << "NormConfReconstrErrorEquidist in W  = "
           << sqrt(PPM->SumOnCommSplit(nrmWConformReconstrAlternative, u.CommSplit())) << endl;
}

double STDGViscoAcousticAssemble::DGError(const Vector &u) const {
//  if (prob->nL() != 0) Exit("check");
  double nrm = 0;
  double h = u.GetMesh().MaxMeshWidth();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    double rho = prob->Rho(*c, c() );
    double rho_inv = 1.0 / prob->Rho(*c, c() );
    double kappa = prob->Kappa(*c, c());
    double kappaInv = 1.0 / kappa;
    double Z = sqrt(rho * kappa);
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    for (int q = 0; q < elem.nQ(); ++q) {
      Point QP = elem.QPoint(q);
      double t = QP.t();
      VectorField V = rho * elem.DtVelocity(q, u) - elem.GradPressure(q, u);
      double P = kappaInv * elem.DtPressure(q, u) - elem.DivVelocity(q, u);
      // P_i = kappaInv_i*elem.DtPressure(q,u) - elem.divVelocity(q, u)+tauInv_i * elem.Pressure_i(q, u);
      double w = elem.QWeight(q);
      V -= prob->F_Velocity(t, *c, QP);
      P -= prob->F_Pressure(t, *c, QP);
      //P_i -= prob->F_Pressure_i(QP);
      nrm += h * w * (rho_inv * V * V + kappa * P * P);
      //nrm += h * w * (kappa * tau_p * P_i);
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
      if (c() < cf()) continue;
      double Zf = sqrt(prob->Rho(*cf, cf()) * prob->Kappa(*cf, cf()));
      double Zinv = 1 / (Z + Zf);
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1, prob->nL());
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
    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 2, prob->nL(), "dg");
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
      if (c_prev == u.cells_end() || c_prev == u.overlap_end()) {
        std::stringstream s;
        s << c();
        Warning("A)) Cell has no previous cell!: " + s.str());
      } else if (c_prev == c) {
        std::stringstream s;
        s << c() << c_prev();
        Warning("B)) Cell has no previous cell!: " + s.str());
      }
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, c_prev, c_prev.Faces() - 1, prob->nL(),
                                                   "dg", c);
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
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 1, prob->nL(), "dg");
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

double STDGViscoAcousticAssemble::L2Error(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    double rho = prob->Rho(*c, c() );
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

double STDGViscoAcousticAssemble::LInfError(const Vector &u) const {
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

double STDGViscoAcousticAssemble::GNError(const Vector &) const {
  return -1.0;
}

double STDGViscoAcousticAssemble::EnergyNorm(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    row r_c = u.find_row(c());
    const double *u_c = u(r_c);
    const double rho = prob->Rho(*c, c());
    vector<double> kappaInv(2 + prob->nL());
    kappaInv[0] = 1.0 / prob->Kappa(*c, c());

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

double STDGViscoAcousticAssemble::VNorm(const Matrix &L, const Vector &u) const {
  Vector Lu = L * u;
  double normLu = L1Norm(Lu);
  double normSpaceIntegral = 0;
  double normSpaceTimeIntegral = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    row r_c = u.find_row(c());
    int c_deg = u.GetDoF().GetDegree(*c).space;
    const double rho = prob->Rho(*c, c() );
    vector<double> kappaInv(2 + prob->nL());
    vector<double> kappaTauInv(2 + prob->nL());
    kappaInv[1] = 1.0 / prob->Kappa_i(*c, c(), 0);
    for (int j = 2; j < 2 + prob->nL(); j++) {
      kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j - 1);
      kappaTauInv[j] = kappaInv[j] / prob->Tau_i(*c, c(), j - 1);
    }
    const double z_K = sqrt(prob->Rho(*c, c() ) * prob->Kappa(*c, c()));
    for (int f = 0; f < c.Faces() - 2; ++f) {
      double face_norm = 0;
      bool bnd = u.OnBoundary(*c, f);
      if (bnd) continue; // TODO check this continue;
      int bnd_id = bnd ? prob->BndID(c.Face(f)) : -1;
      cell cf = bnd ? c : u.find_neighbour_cell(c, f);
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, f, prob->nL());
      int cf_deg = u.GetDoF().GetDegree(*cf).space;
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);
      row r_cf = u.find_row(cf());
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1, prob->nL());
      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        double w = felem.QWeight(q);
        VectorField Nq = felem.QNormal(q);
        Point QP = felem.QPoint(q);

        double face_jump_p = felem.Pressure(q, u(r_c)) - felem_1.Pressure(q1, u(r_cf));
        double face_jump_v = Nq * (felem.Velocity(q, u(r_c)) - felem_1.Velocity(q1, u(r_cf)));
        normSpaceIntegral += w * (abs(face_jump_p) + abs(face_jump_v));
      }
    }

    cell c_prev = u.find_previous_cell(c);
    row r_cp = u.find_row(c_prev());
    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 2, prob->nL(), "dg");
    int c_prev_deg = u.GetDoF().GetDegree(*c_prev).space;
    SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, c_prev, c_prev.Faces() - 1, prob->nL(),
                                                 "dg",
                                                 c);

    if (c == c_prev) {
      continue; // TODO check this continue
    }

    for (int q = 0; q < felem.nQ(); ++q) {
      double w = felem.QWeight(q);
      int q1 = find_q_id(felem.QPoint(q), felem_1);

      double face_jump_p = felem.Pressure(q, u(r_c)) - felem_1.Pressure(q1, u(r_cp));
      VectorField face_jump_v = felem.Velocity(q, u(r_c)) - felem_1.Velocity(q1, u(r_cp));
      double face_jump_v_value = abs(face_jump_v[0]) + abs(face_jump_v[1]);
      normSpaceIntegral += w * (kappaInv[0] * abs(face_jump_p) + rho * face_jump_v_value);
    }

  }


  return PPM->SumOnCommSplit(normLu + normSpaceIntegral + normSpaceTimeIntegral / 2.0,
                             u.CommSplit());
}

double STDGViscoAcousticAssemble::L1Norm(const Vector &u) const {
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

double STDGViscoAcousticAssemble::L2Norm(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    double rho = prob->Rho(*c, c() );
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

double STDGViscoAcousticAssemble::L2SpaceNormAtTime(const Vector &u, double time) const {
  double nrm = 0.0;
  if (time == u.GetMesh().GetEndTime()) {
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != time) continue;
      double rho = prob->Rho(*c, c() );
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
      if (!(c.min() == 0 && time == 0) || time <= c.min() || c.max() < time) continue;
      double rho = prob->Rho(*c, c() );
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

double STDGViscoAcousticAssemble::L2SpaceNormAtTimeError(const Vector &u, double time) const {
  double nrm = 0.0;
  const double T = u.GetMesh().GetEndTime();
  if (time == T) {
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != T) continue;
      double rho = prob->Rho(*c, c() );
      double kappaInv = 1.0 / prob->Kappa(*c, c());
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 1, prob->nL(), "dg");
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
      double rho = prob->Rho(*c, c() );
      double kappaInv = 1.0 / prob->Kappa(*c, c());
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 2, prob->nL(), "dg");
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

/*
double STDGViscoAcousticAssemble::L2SpaceNormAtTimeError(const Vector &u, double time) const {
  double nrm = 0.0;
  double normP = 0;
  double normV1 = 0;
  double normV2 = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    if (time == 0) {
      if (c.min() > 0) {
        continue;
      }
    } else {
      if (time <= c.min() || c.max() < time) {
        continue;
      }
    }
    double rho = prob->Rho(*c, c() );
    double kappaInv = 1.0 / prob->Kappa(*c, c());
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL(), true);
    double localtime = c.TimeCell().GlobalToLocal(time)[0];
    for (int q = 0; q < elem.nSpaceQ(); ++q) {
      double w = elem.SpaceQWeight(q);
      Point localQP = elem.SpaceQPoint(q).CopyWithT(localtime);
      Point globalQP = c.LocalToGlobal(localQP);
      VectorField V = elem.VelocityLocal(localQP, u) - prob->ut_Velocity(globalQP, c);
      double P = elem.PressureLocal(localQP, u) - prob->ut_Pressure(globalQP, c);

      //mout << "QP " << globalQP << " P " << P << " V " << V << endl;

      nrm += w * (kappaInv * P * P + rho * V * V);
      normP += w * P * P;
      normV1 += w * V[0] * V[0];
      normV2 += w * V[1] * V[1];
    }
  }
  //normP = PPM->Sum(normP);
  //normV1 = PPM->Sum(normV1);
  //normV2 = PPM->Sum(normV2);
  //mout << time << " " << DOUT(normP) << DOUT(normV1) <<DOUT(normV2) << endl;


  mout << " time  " << time << " nrm2 " << nrm << endl;



  return sqrt(PPM->Sum(nrm));
}
*/

double STDGViscoAcousticAssemble::DGSemiNorm(const Vector &u) const {
  double nrm = 0.0;

  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    row r_c = u.find_row(c());
    int c_deg = u.GetDoF().GetDegree(*c).space;
    const double rho = prob->Rho(*c, c() );
    vector<double> kappaInv(2 + prob->nL());
    vector<double> kappaTauInv(2 + prob->nL());
    kappaInv[0] = 1.0 / prob->Kappa(*c, c());
    for (int j = 1; j < 2 + prob->nL(); j++) {
      kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j);
      kappaTauInv[j] = kappaInv[j] / prob->Tau_i(*c, c(), j);
    }
    const double z_K = sqrt(prob->Rho(*c, c() ) * prob->Kappa(*c, c()));
    for (int f = 0; f < c.Faces() - 2; ++f) {
      double face_norm = 0;
      bool bnd = u.OnBoundary(*c, f);
      if (bnd) continue; // TODO check this continue;
      int bnd_id = bnd ? prob->BndID(c.Face(f)) : -1;
      cell cf = bnd ? c : u.find_neighbour_cell(c, f);
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, f, prob->nL());
      int cf_deg = u.GetDoF().GetDegree(*cf).space;
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);
      row r_cf = u.find_row(cf());
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1, prob->nL());

      const double z_Kf = sqrt(prob->Rho(*cf, cf()) * prob->Kappa(*cf, cf()));
      const double alpha1 = 1.0 / (z_K + z_Kf);
      const double alpha2 = z_Kf * z_K * alpha1;
      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        double w = felem.QWeight(q);
        VectorField Nq = felem.QNormal(q);
        Point QP = felem.QPoint(q);

        double face_jump_p = felem.Pressure(q, u(r_c)) - felem_1.Pressure(q1, u(r_cf));
        double face_jump_v = Nq * (felem.Velocity(q, u(r_c)) - felem_1.Velocity(q1, u(r_cf)));
        face_jump_p = alpha1 * face_jump_p * face_jump_p;
        face_jump_v = alpha2 * face_jump_v * face_jump_v;
        nrm += w * (face_jump_p + face_jump_v);
      }
    }
    cell c_prev = u.find_previous_cell(c);
    row r_cp = u.find_row(c_prev());
    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 2, prob->nL(), "dg");
    int c_prev_deg = u.GetDoF().GetDegree(*c_prev).space;
    SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, c_prev, c_prev.Faces() - 1, prob->nL(),
                                                 "dg", c);

    if (c == c_prev) {
      continue; // TODO check this continue
    }

    for (int q = 0; q < felem.nQ(); ++q) {
      double w = felem.QWeight(q);
      int q1 = find_q_id(felem.QPoint(q), felem_1);

      double face_jump_p = felem.Pressure(q, u(r_c)) - felem_1.Pressure(q1, u(r_cp));
      VectorField face_jump_v = felem.Velocity(q, u(r_c)) - felem_1.Velocity(q1, u(r_cp));
      face_jump_p = kappaInv[0] * face_jump_p * face_jump_p;
      double face_jump_v_value = rho * face_jump_v * face_jump_v;
      nrm += w * (face_jump_p + face_jump_v_value);
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()) / 2);
};

double STDGViscoAcousticAssemble::DGNormError(const Vector &u) {
//
  THROW("needs to be reimplemented")
  return 1.0;
//  std::shared_ptr<ueberclass> temp_prob = prob;
//  prob = CreateViscoAcousticProblemShared("Zero");
//  Matrix M_zero(u);
//  Vector rhs_zero(0.0, u);
//  System(M_zero, rhs_zero);
//  prob = temp_prob;
//
//  double normsq = 0;
//  Vector Lu = M_zero * u;
//  double buu = Lu * u;
//
//  double h = u.GetMesh().MaxMeshWidth();
//
//  for (cell c = u.cells(); c != u.cells_end(); c++) {
//
//    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
//
//    double rho = prob->Rho(*c, c() );
//    double rho_inv = 1.0 / prob->Rho(*c, c() );
//    double kappa = prob->Kappa(*c, c());
//    double kappa_inv = 1.0 / prob->Kappa(*c, c());
//    for (int q = 0; q < elem.nQ(); q++) {
//      Point QP = elem.QPoint(q);
//      double w = elem.QWeight(q);
//
//      double dtP = elem.DtPressure(q, u);
//      VectorField gradP = elem.GradPressure(q, u);
//      VectorField dtV = elem.DtVelocity(q, u);
//      double divV = elem.DivVelocity(q, u);
//
//      double P = rho * dtP - divV - prob->F_Pressure(t, *c, QP);
//      VectorField V = kappa_inv * dtV - gradP - prob->F_Velocity(t, *c, QP);
//
//      normsq += w * (rho_inv * P * P + kappa * V * V);
//    }
//  }
//  return sqrt(buu + PPM->SumOnCommSplit(h * normsq, u.CommSplit()));
}


double STDGViscoAcousticAssemble::DGNorm(const Vector &u, const Matrix &L) const {
  double normsq = 0;
  Vector Lu = L * u;
  double buu = Lu * u;

  for (cell c = u.cells(); c != u.cells_end(); c++) {
    double h = max(0.0, dist(c[0], c[c.Corners() - 1]));
    for (int cs = 0; cs < c.Corners() - 1; ++cs)
      h = max(h, dist(c[cs], c[cs + 1]));

    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());

    double rho = prob->Rho(*c, c() );
    double rho_inv = 1.0 / prob->Rho(*c, c() );
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

void STDGViscoAcousticAssemble::DualRHS_linear(Vector &RHS,
                                               const Vector &U,
                                               Point a,
                                               Point b) const {

  RHS = 0;

  if (a.t() != b.t()) { // volume evaluation
    Exit("TODO")
  } else { // evaluation at some time point
    for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {
      if (c.max() != a.t()) continue;

      SpaceTimeViscoAcousticDGTElement elem(RHS, c, prob->nL());

      row r_RHS = RHS.find_row(c());

      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        double w = elem.SpaceQWeight(q);
        Point QP = elem.QPoint(q).CopyWithT(0.0);

        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;

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

void STDGViscoAcousticAssemble::DualRHS_quadratic(Vector &RHS,
                                                  const Vector &U,
                                                  Point a,
                                                  Point b) const {
  Exit("TODO")
  Exit("Not implemented");
}

double STDGViscoAcousticAssemble::MminushalfLuMinusF(const cell &c,
                                                     const Vector &u) const {

  double normsq = 0;
  SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
  double rho = prob->Rho(*c, c() );
  double rho_inv = 1.0 / prob->Rho(*c, c() );
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

double STDGViscoAcousticAssemble::MhalfLu_exactMinusF(LevelPair levels) const {
  Vector res(disc, levels);
  res.Clear();
  for (cell c = res.cells(); c != res.cells_end(); c++) {
    double normsq = 0.0;
    SpaceTimeViscoAcousticDGTElement elem(res, c, prob->nL());
    double rho = prob->Rho(*c, c() );
    double rho_inv = 1.0 / prob->Rho(*c, c() );
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
  //printVTK_Eta(v, *this, 100+U.Level().space);
  return res.norm();
}


double STDGViscoAcousticAssemble::ResidualErrorEstimateJumpCell(const cell &c,
                                                                const Vector &u) const {
  double eta = 0.0;
  double rho = prob->Rho(*c, c() );
  double rho_inv = 1.0 / prob->Rho(*c, c() );
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
  SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, c.Faces() - 2, prob->nL(), "dg");
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
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, c_prev, c_prev.Faces() - 1, prob->nL(),
                                                   "dg", c);
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
  if (c.max() == u.GetMesh().GetEndTime()) return sqrt(eta);
  face f = u.find_face(c.Face(c.Faces() - 1));
  if (!u.has_cell_or_overlap_cell(f.Right())) {
    return sqrt(eta);
  }
  cell c_next = u.find_cell_or_overlap_cell(f.Right());

  SpaceTimeViscoAcousticDGTFaceElement felem_0(*disc, u, c, c.Faces() - 1, prob->nL(), "dg");
  SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, c_next, c_next.Faces() - 2, prob->nL(),
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

double STDGViscoAcousticAssemble::L2ReliableResidualErrorEstimateJumpCell(const cell &c,
                                                                          const Vector &u,
                                                                          const Vector &conf) const {
  SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());

  double Mhalf_u_minus_u_conf = 0.0;
  double Mhalfinv_Lu_minus_f = 0.0;

  double rho = prob->Rho(*c, c() );
  double rho_inv = 1.0 / prob->Rho(*c, c() );
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
    Mhalf_u_minus_u_conf += w * (rho * u_minus_u_conf_P * u_minus_u_conf_P +
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
      Mhalf_u_conf0_minus_u0 += w * (rho * u_conf0_minus_u0_P * u_conf0_minus_u0_P +
                                     kappa_inv * u_conf0_minus_u0_V * u_conf0_minus_u0_V);
    }
    normsq_0 += Mhalf_u_conf0_minus_u0;
  }

  double T = u.GetMesh().GetEndTime();
  double normsq = Mhalf_u_minus_u_conf + 4 * T * T * (Mhalfinv_Lu_minus_f + normsq_0);

  return sqrt(normsq);
};

double STDGViscoAcousticAssemble::DGErrorEstimateJumpCell(const cell &c,
                                                          const Vector &u,
                                                          const Vector &conf,
                                                          double eta_res_c,
                                                          double eta_conf_c) const {
  STDGDGViscoAcousticElement elem(u, c, prob->nL());
  double kappa = prob->Kappa(*c, c());
  double kappa_inv = 1.0 / prob->Kappa(*c, c());
  double rho = prob->Rho(*c, c() );
  double rho_inv = 1.0 / prob->Rho(*c, c() );

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
      Mhalf_uT_minus_u_confT += w * (rho * u_confT_minus_uT_P * u_confT_minus_uT_P +
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
      Mhalf_uT_minus_u_confT += w * (rho * u_confT_minus_uT_P * u_confT_minus_uT_P +
                                     kappa_inv * u_confT_minus_uT_V * u_confT_minus_uT_V);
    }
  }
  double dg_normsq = eta_res_c * eta_res_c;
  dg_normsq += 2 * eta_conf_c * sqrt(Mhalfinv_Lu_conf_minus_f) + Mhalf_uT_minus_u_confT;
  return sqrt(dg_normsq);
};


void
STDGViscoAcousticAssemble::ConformingLagrangeInterpolation(const Vector &u, Vector &conf,
                                                           int degree) const {
  mout.StartBlock("ConfLagInterpolation");

  //pout << u << endl;
  //pout << u.GetMesh() << endl;

  auto meshes = MeshesCreator().
      WithMeshName(u.GetMesh().Name()).
      WithPLevel(u.GetDisc().GetMeshes().pLevel()).
      WithLevel(u.Level().space).
      CreateUnique();
  auto ld = std::make_shared<LagrangeDiscretization>(*meshes, degree, prob->Dim());
  //ld().procSets.CheckConsistency();
  //ld().GetProcSets().CheckConsistency();

  Vector l_vec(0.0, ld, u.Level());
  Vector adds(0.0, ld, u.Level());
  l_vec.Clear();
  adds.Clear();
  conf.Clear();

  for (cell c = u.cells(); c != u.cells_end(); c++) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    vector<Point> nps = l_vec.GetDoF().GetNodalPoints(*c);
    for (int component = 0; component < prob->Dim(); component++) {
      for (Point np: nps) {
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
        for (int s = 0; s < conf.GetDoF().NodalPointsLocal(deg.space, pi); s++) {
          conf(elem.r, cnt) += l_elem.Value(c.GlobalToLocal(nps[cnt]), l_vec, pi);
          cnt++;
        }
      }
    }
  }
  conf.Accumulate();
  mout.EndBlock(verbose > 0);
}

void STDGViscoAcousticAssemble::ConformingInterpolation(const Vector &u, Vector &conf) const {

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
      for (auto &[point, components]: conf_data) {
        buffer.Send(q) << point << number_of_additions[point] << components.size();
        for (double component: components) {
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
  for (auto &[point, components]: conf_data) {
    for (double &component: components) {
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

void STDGViscoAcousticAssemble::ConformingProjectionWithPenalty(const Vector &u,
                                                                Vector &u_conf,
                                                                Vector &r_conf) const {


  if (penalty == 0) return;
  Matrix A_conf(u);
  A_conf = 0;
  r_conf = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    row r_c = u.find_row(c());
    //const double *u_c = u(r_c);

    DGRowEntries A_cc(A_conf, *c, *c, false);

    const double rho = prob->Rho(*c, c() );
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
        r_conf(r_c, i) += w * (rho * (V * V_i)
                               + kappaInv[0] * (P * P_i));
        for (int j = 0; j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          VectorField V_j = elem.Velocity(q, j);
          double P_j = elem.Pressure(q, j);
          A_cc(i, j) += w * ((rho * (V_j * V_i)
                              + kappaInv[0] * (P_j * P_i)));
        }
      }
    }

    // penalty for the jumps in time
    cell c_prev = A_conf.find_previous_cell(c);
    if (c != c_prev) {
      SpaceTimeViscoAcousticDGTFaceElement felem_t(*disc, A_conf, c, c.Faces() - 2, prob->nL(),
                                                   "dg");
      SpaceTimeViscoAcousticDGTFaceElement felem_t1(*disc, A_conf, c_prev, c_prev.Faces() - 1,
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
          if (felem_t.is_zero_test(q, i)) continue;
          int var_i = felem_t.variable(i);
          VectorField V_i = felem_t.Velocity(q, i);
          double P_i = felem_t.Pressure(q, i);

          for (int j = 0; j < felem_t.j_dimension(); ++j) {
            int var_j = felem_t.variable(j);
            if (felem_t.is_zero(q, j) || var_i != var_j) continue;

            VectorField V_j = felem_t.Velocity(q, j);
            double P_j = felem_t.Pressure(q, j);
            A_cc(i, j) += penalty * w * (rho * (V_j * V_i)
                                         + kappaInv[0] * P_j * P_i
            );
          }
          for (int j = 0; j < felem_t1.j_dimension(); ++j) {
            int var_j = felem_t1.variable(j);
            if (felem_t1.is_zero(q, j) || var_i != var_j) continue;
            VectorField V1_j = felem_t1.Velocity(q, j);
            double P1_j = felem_t1.Pressure(q, j);
            double a = penalty * w * (rho * (V1_j * V_i)
                                      + kappaInv[0] * P1_j * P_i
            );
            A_cp(i, j) -= a;
            A_pc(j, i) -= a;
          }
        }
        for (int i = 0; i < felem_t1.i_dimension(); ++i) {
          if (felem_t1.is_zero_test(q, i)) continue;
          int var_i = felem_t1.variable(i);
          VectorField V1_i = felem_t1.Velocity(q, i);
          double P1_i = felem_t1.Pressure(q, i);
          for (int j = 0; j < felem_t1.j_dimension(); ++j) {
            int var_j = felem_t1.variable(j);
            if (felem_t1.is_zero(q, j) || var_i != var_j) continue;
            VectorField V1_j = felem_t1.Velocity(q, j);
            double P1_j = felem_t1.Pressure(q, j);
            A_pp(i, j) += penalty * w * (rho * (V1_j * V1_i)
                                         + kappaInv[0] * P1_j * P1_i
            );
          }
        }
      }
    }

    // penalty for  the jumps in space
    for (int f = 0; f < c.Faces() - 2; ++f) {
      if (u.OnBoundary(*c, f)) continue;
      cell cf = u.find_neighbour_cell(c, f);
      SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, f, prob->nL());
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);
      SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1, prob->nL());
      DGRowEntries A_cf(A_conf, *c, *cf, false);
      DGRowEntries A_fc(A_conf, *cf, *c, false);
      DGRowEntries A_ff(A_conf, *cf, *cf, false);

      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        VectorField Nq = felem.QNormal(q);
        int q1 = find_q_id(felem.QPoint(q), felem_1);

        for (int i = 0; i < felem.i_dimension(); ++i) {
          int var_i = felem.variable(i);
          if (felem.is_zero_test(q, i)) continue;
          double VN_i = felem.Velocity(q, i) * Nq;
          double P_i = felem.Pressure(q, i);
          for (int j = 0; j < felem.i_dimension(); ++j) {
            int var_j = felem.variable(j);
            if (felem.is_zero_test(q, j) || var_i != var_j) continue;
            double VN_j = felem.Velocity(q, j) * Nq;
            double P_j = felem.Pressure(q, j);
            A_cc(i, j) += penalty * w * ((rho * (VN_j * VN_i)
                                          + kappaInv[0] * (P_j * P_i)));
          }
          for (int j = 0; j < felem_1.i_dimension(); ++j) {
            int var_j = felem_1.variable(j);
            if (felem_1.is_zero_test(q1, j) || var_i != var_j) continue;
            double VN_j = felem_1.Velocity(q1, j) * Nq;
            double P_j = felem_1.Pressure(q1, j);
            A_cf(i, j) -= penalty * w * ((rho * (VN_j * VN_i)
                                          + kappaInv[0] * (P_j * P_i)));
          }
        }
      }
    }
  }
  r_conf -= A_conf * u_conf;

  mout << "Solving ConformingProjection with penality=" + std::to_string(penalty) << endl;
  LinearSolver *S = GetLinearSolverByPrefix("ConformingProjection");

  u_conf += (*S)(A_conf) * r_conf;
}

double STDGViscoAcousticAssemble::DualErrorEstimateJumpCell(const cell &c,
                                                            const Vector &u,
                                                            const Vector &u_star) const {
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

  double rho = prob->Rho(*c, c() );
  double kappaInv = 1.0 / prob->Kappa(*c, c());

  for (int q = 0; q < elem.nQ(); ++q) {
    VectorField V = 0.0;
    double P = 0.0;
    double w = elem.QWeight(q);
    Point QP = elem.QPoint(q);
    double t = QP.t();

    for (int j = 0; j < elem.j_dimension(); ++j) {
      VectorField phi_V_j = -1.0 * rho * elem.DtVelocity(q, j) + elem.GradPressure(q, j);
      double phi_P_j = -1.0 * kappaInv * elem.DtPressure(q, j) + elem.DivVelocity(q, j);
      V += u_r[j] * phi_V_j;
      P += u_r[j] * phi_P_j;
    }

    VectorField V_rhs = prob->F_Velocity(t, *c, QP);;
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
      dual_V_projection[l] += elem.QWeight(q) * u_star_r[i] * elem.Velocity(q, i);
      dual_P_projection[l] += elem.QWeight(q) * u_star_r[i] * elem.Pressure(q, i);
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
      if (pfc == pfcf) break;
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
        values_dual_V_cf += dual_V_projection_cf[l] * felem.QNormal(0) * timeValue;
        values_dual_P_cf += dual_P_projection_cf[l] * timeValue;
      }

      V_dual_Fluxes[f] += f_vol * pow(values_dual_V_cf - values_dual_V, 2)
                          * timeQuad.Weight(q) * I;
      P_dual_Fluxes[f] += f_vol * pow(values_dual_P_cf - values_dual_P, 2)
                          * timeQuad.Weight(q) * I;
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
      if (pfc == pfcf) break;
    }

    SpaceTimeViscoAcousticDGTFaceElement felem(*disc, u, c, f, prob->nL(), "dgdg");
    SpaceTimeViscoAcousticDGTFaceElement felem_1(*disc, u, cf, f1, prob->nL(), "dgdg", c);
    SpaceTimeViscoAcousticDGTElement elem_1(u, cf, prob->nL());

    double eta_v = 0.;
    double fxI = 0;

    for (int q = 0; q < felem.nQ(); ++q) {
      int q1 = q;
      //if (f<c.Faces()-2) find_q_id(felem.QPoint(q),*felem_1);

      VectorField V_Flux = 0.0;
      double P_Flux = 0.0;

      double w = felem.QWeight(q);
      //VectorField Nq = felem.QNormal(q);

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
        dual_V_projection_cf += elem_1.QWeight(q) * u_star_rf[i] * elem_1.Velocity(q, i);
        dual_P_projection_cf += elem_1.QWeight(q) * u_star_rf[i] * elem_1.Pressure(q, i);
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
    eta_f += sqrt(V_Fluxes[f]) * sqrt(V_dual_Fluxes[f])
             + sqrt(P_Fluxes[f]) * sqrt(P_dual_Fluxes[f]);
  }
  eta = sqrt(rho_1_P) * sqrt(h) * sqrt(omega_1_P)
        + sqrt(rho_1_V) * sqrt(h) * sqrt(omega_1_V) + 0.5 * eta_f;

  return eta;
}

double STDGViscoAcousticAssemble::Goal_Functional(const Vector &u,
                                                  Point a,
                                                  Point b) const {
  Scalar e = 0.0;

  if (a.t() != b.t()) { // volume evaluation
    Exit("ToDo");
  } else { // evaluation at some time point
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != a.t()) continue;
      SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
      row r = u.find_row(c());
      const int dim = c.dim();
      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        Point QP = elem.QPoint(q).CopyWithT(0.0);
        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;
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

double STDGViscoAcousticAssemble::Energy(const Vector &u, Point a, Point b) const {
  double e = 0;
  if (a.t() != b.t()) { // volume evaluation
    Exit("TODO")
  } else { // evaluation at some time point
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != a.t()) continue;
      SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
      row r = u.find_row(c());
      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        Point QP = elem.QPoint(q).CopyWithT(0.0);
        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;
        const double w = elem.SpaceQWeight(q);

        VectorField V = elem.Velocity(q, u);
        double P = elem.Pressure(q, u);
        DampingVector D = elem.DampingPressure(q, u);
        e += 0.5 * w * V * V;
        e += 0.5 * w * P * P;
        for (int k = 0; k <= D.size(); ++k)
          e += 0.5 * w * D[k] * D[k];
      }
    }
  }
  return abs(PPM->SumOnCommSplit(e, u.CommSplit()));
}

double STDGViscoAcousticAssemble::L2ScalarProduct(const Vector &a, const Vector &b) const {
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

void STDGViscoAcousticAssemble::get_projected_exact_solution(Vector &Exact_Solution) const {
  Exact_Solution.Clear();
  std::unordered_map<DegreePair, RMatrix> InvertedMassMatrixCache_t;
  std::unordered_map<DegreePair, SpaceTimeCellQuadrature> SpaceTimeQuadCache;
  for (cell c_t = Exact_Solution.cells(); c_t != Exact_Solution.cells_end(); c_t++) {
    row row_t = Exact_Solution.find_row(c_t());
    DegreePair deg_t = Exact_Solution.GetDoF().GetDegree(*c_t);
    auto &space_shape_t = disc->GetCellShape(deg_t.space);
    auto &time_shape_t = disc->GetTimeShape(deg_t.time);
    int N_t = time_shape_t.size() * prob->Dim() * space_shape_t.size();


    auto elem_q = SpaceTimeQuadCache.find(deg_t);
    if (elem_q == SpaceTimeQuadCache.end()) {
      const Quadrature &s_quad = disc->GetCellQuad(deg_t.space);
      const Quadrature &t_quad = disc->GetTimeQuad(deg_t.time);
      elem_q = SpaceTimeQuadCache.insert(
          {deg_t, SpaceTimeCellQuadrature(c_t, s_quad, t_quad)}).first;
    }
    SpaceTimeCellQuadrature &quad_t = elem_q->second;


    //SpaceTimeCellQuadrature quad_t(c_t, disc->GetCellQuad(deg_t.space),
    //                               disc->GetTimeQuad(deg_t.time));


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
                  int col = sc2 + (k + tc2 * prob->Dim()) * space_shape_t.size();
                  double ss_t2 = space_shape_t(local, sc2);
                  double v_t2 = ss_t2 * ts_t2;
                  InvertedMassMatrix(row, col) += quad_t.qWeight[q] * v_t * v_t2;
                }
              }
            }
          }
        }
      }
      InvertedMassMatrix.Invert();
      elem = InvertedMassMatrixCache_t.insert({deg_t, InvertedMassMatrix}).first;
    }
    RMatrix InvertedMassMatrix_t = elem->second;

    //RMatrix MassMatrix_ft(N_t, N_t);
    RVector integral_t(N_t);
    for (int q = 0; q < quad_t.size(); q++) {
      Point local = quad_t.qPointLocal[q];
      Point global = c_t.LocalToGlobal(local);
      for (int tt = 0; tt < time_shape_t.size(); ++tt) {
        double value_time_shape = time_shape_t(Point(local.t()), tt);
        for (int st = 0; st < space_shape_t.size(); ++st) {
          double value_space_shape = space_shape_t(local, st);
          double v_t = value_space_shape * value_time_shape;
          //for (int tf = 0; tf < time_shape_t.size(); ++tf) {
          //double value_time_shape_f = time_shape_t(Point(local.t()), tf);
          //for (int sf = 0; sf < time_shape_t.size(); ++sf) {
          //double value_space_shape_f = time_shape_t(local, sf);
          for (int k = 0; k < prob->Dim(); ++k) {
            COMPONENT comp = prob->GetComponents()[k];
            int t_idx = st + (k + tt * prob->Dim()) * space_shape_t.size();
            //int f_idx = sf + (k + tf * prob->SpaceDim()) * space_shape_t.size();
            double v_f = prob->ut(global.t(), global, *c_t, comp);//value_space_shape_f * value_time_shape_f;
            //MassMatrix_ft[f_idx][t_idx] += quad_t.qWeight[q] * v_t * v_f;
            integral_t[t_idx] += quad_t.qWeight[q] * v_t * v_f;
          }
          //}
          //}
        }
      }
    }

    //RVector value_f(N_f);
    //for (int i = 0; i < value_f.size(); i++) {
    //  value_f[i] = from(row_f)[i];
    //}
    //mout << MassMatrix_ft << endl;
    //mout << InvertedMassMatrix_t << endl;
    //return;
    //RVector integral_t(MassMatrix_ft * value_f);
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
  //mout << Exact_Solution << endl;
}

void STDGViscoAcousticAssemble::PlotSingleSolution(const Vector &u, const string &filename) const {
  auto vtuDisc = std::make_shared<STDiscretizationT_DGDG<>>(u.GetDisc().GetMeshes(), DegreePair{0, 0}, prob->Dim());
  Vector U_vtu(0.0, vtuDisc, u.Level());
  Vector PSumComponent(0.0, vtuDisc, u.Level());
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(u, c, prob->nL());
    row r = U_vtu.find_row(c());
    for (size_t i = 0; i < elem.GetComponentCount(); i++) {
      U_vtu(r, i) = elem.EvaluateComponentGlobal(c(), u, i);
    }
    for (size_t i = 0; i < prob->nL() + 1; i++) {
      PSumComponent(r, 0) += U_vtu(r, SpaceDimension + i);
    }
  }

  VtuPlot plot{filename};
  vector<int> velocity_indices(SpaceDimension);
  std::iota(begin(velocity_indices), end(velocity_indices), 0);
  plot.AddData("Velocity", U_vtu, velocity_indices);
  plot.AddData("p_sum", PSumComponent, 0);
  for (size_t i = 0; i < prob->nL() + 1; i++) {
    plot.AddData("p"+std::to_string(i), U_vtu, SpaceDimension + i);
  }
  plot.PlotFile();
}
void STDGViscoAcousticAssemble::PlotParameters(const Vector &u, std::string suffix) const {
  auto paramDisc = std::make_shared<const STDiscretizationT_DGDG<>>(u.GetDisc().GetMeshes(), DegreePair{0, 0}, prob-> nL() + 2);
  Vector rhoAndKappa(0.0, paramDisc, u.Level());

  for (cell c = rhoAndKappa.cells(); c != rhoAndKappa.cells_end(); ++c) {
    rhoAndKappa(c())[0] = GetProblem().Rho(*c, c());
    rhoAndKappa(c())[1] = GetProblem().Kappa(*c, c());
    for (size_t i = 2; i < prob->nL() + 2; i++) {
      rhoAndKappa(c())[i] = GetProblem().Kappa_i(*c, c(), i - 2);
    }
  }
  VtuPlot plot("Parameters" + suffix);
  plot.AddData("rho", rhoAndKappa, 0);
  plot.AddData("kappa", rhoAndKappa, 1);
  for (size_t i = 2; i < prob->nL() + 2; i++) {
    plot.AddData("kappa"+std::to_string(i-2), rhoAndKappa,  i);
  }
  plot.PlotFile();
}

void STDGViscoAcousticAssemble::PlotDualSolution(const Vector &U_dual, const Vector &Eta, int run) {
  const Mesh &STMesh = U_dual.GetMesh();

  int probDim = prob->Dim();
  Vector U_vtk_dual(0.0, U_dual);

  for (cell c = U_dual.cells(); c != U_dual.cells_end(); ++c) {
    row r = U_dual.find_row(c());
    DegreePair deg = U_dual.GetDoF().GetDegree(*c);
    int shift = r.n() * deg.time / (deg.time + 1);
    Point localPT;
    transformPointGlobalToLocal(c(), localPT, c);
    int NX = 0;
    double tmp = 0.0;
    for (int p = 0; p < probDim; ++p) {
      const Shape &shape(disc->GetCellShape(deg.space));
      for (int n = 0; n < U_dual.GetDoF().NodalPointsLocal(deg.space, p); n++) {
        tmp += shape(localPT, n) * U_dual(r)[shift + NX + n];
      }
      U_vtk_dual(r)[p] = tmp;
      NX += U_dual.GetDoF().NodalPointsLocal(deg.space, p);
    }
  }
  VtuPlot plot(STMesh);
  for (int i = 0; i < probDim; ++i) {
    COMPONENT comp = prob->GetComponents()[i];
    string name = "Dual_" + to_string(comp) + "_" + to_string(run);
    plot.AddData(name, U_vtk_dual, i);
  }
  plot.PlotFile("Dual");
  VtuPlot plot_eta(STMesh);
  string name = "Est_Error_" + to_string(run);
  plot_eta.AddData(name, Eta, 0);
  plot_eta.PlotFile(name);
}

Vector STDGViscoAcousticAssemble::CalculateMaterialUpdate(const Vector &material,
                                                          const Vector &forwardSolution,
                                                          const Vector &backwardSolution) {
  Vector gradient(0.0, material);
  for (cell c = material.cells(); c != material.cells_end(); c++) {
    SpaceTimeViscoAcousticDGTElement elem(forwardSolution, c, prob->nL());
    double spaceArea = elem.Area() / (c.max() - c.min());
    for (int q = 0; q < elem.nQ(); q++) {
      double w = elem.QWeight(q) / spaceArea;
      double kappa = material(c(), 0);
      gradient(c(), 0) -= w / (kappa * kappa) * elem.DtPressure(q, forwardSolution) *
                          elem.Pressure(q, backwardSolution);
      gradient(c(), 1) +=
          w * elem.DtVelocity(q, forwardSolution) * elem.Velocity(q, backwardSolution);
    }
  }
  std::map<Point, double> gradientValuesKappa;
  std::map<Point, double> gradientValuesRho;
  for (cell c = material.cells(); c != material.cells_end(); c++) {
    gradientValuesKappa[c.SpaceCell().Center()] += gradient(c(), 0);
    gradientValuesRho[c.SpaceCell().Center()] += gradient(c(), 1);
  }

  ExchangeBuffer buffer;
  ExchangeBuffer bufferReceiver;
  if (!PPM->Master()) {
    buffer.Send(0) << gradientValuesKappa << gradientValuesRho;
  }
  buffer.Communicate();

  if (PPM->Master()) {
    for (int i = 1; i < PPM->Size(); i++) {
      std::map<Point, double> gradientValuesKappaLocal;
      std::map<Point, double> gradientValuesRhoLocal;
      buffer.Receive(i) >> gradientValuesKappaLocal >> gradientValuesRhoLocal;
      for (auto &[spacePoint, value]: gradientValuesKappaLocal) {
        gradientValuesKappa[spacePoint] += value;
      }
      for (auto &[spacePoint, value]: gradientValuesRhoLocal) {
        gradientValuesRho[spacePoint] += value;
      }

    }
    for (int i = 1; i < PPM->Size(); i++) {
      bufferReceiver.Send(i) << gradientValuesKappa << gradientValuesRho;
    }
  }
  bufferReceiver.Communicate();

  if (!PPM->Master()) {
    gradientValuesKappa.clear();
    gradientValuesRho.clear();
    bufferReceiver.Receive(0) >> gradientValuesKappa >> gradientValuesRho;
  }

  for (cell c = material.cells(); c != material.cells_end(); c++) {
    gradient(c(), 0) = gradientValuesKappa[c.SpaceCell().Center()];
    gradient(c(), 1) = gradientValuesRho[c.SpaceCell().Center()];
  }

  return gradient;
}