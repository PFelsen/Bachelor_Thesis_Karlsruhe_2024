#include "STPGViscoAcousticAssemble.hpp"
#include "Results.hpp"
#include "LinearSolver.hpp"

using namespace std;

void STPGViscoAcousticAssemble::System(Matrix &M, Vector &RHS) const {
  Date Start;
  mout.PrintInfo("TDGViscoAcousticAssemble", verbose,
                 PrintInfoEntry<int>("Problem Size", M.pSize(), 2),
                 PrintInfoEntry<Date>("start assemble", Start, 2)
  );
  t_cell = Time();
  t_face = Time();

  RHS.Clear(); // RHS = 0;
  if (prob->hasInitialData()) {
    M = 0;
    Vector u_0(0.0, RHS);
    u_0.SetAccumulateFlag(false);
    for (cell c = M.cells(); c != M.cells_end(); ++c) {
      cell c_prev = M.find_previous_cell(c);
      if (c_prev != c) continue;
      Date Start_cell;
      SpaceTimeViscoAcousticElement elem(*disc, M, c, prob->nL());
      row r_c = M.find_row(c());
      DegreePair deg_c = RHS.GetDoF().GetDegree(*c);
      DGRowEntries M_c(M, *c, *c, false);
      const double rho = prob->Rho(*c, c());
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
          VectorField V_i = elem.VelocityTestSpace(q, i);
          double P_i = elem.PressureTestSpace(q, i);
          for (int j = 0; j < elem.j_dimension() - r_c.n(); ++j) {
            int var_j = elem.variable(j);
            VectorField DtV_j = elem.DtVelocity(q, j);
            Scalar DivV_j = elem.DivVelocity(q, j);
            Scalar P_j = elem.Pressure(q, j);
            Scalar DtP_j = elem.DtPressure(q, j);
            VectorField GradP_j = elem.GradPressure(q, j);

//            double kappaInv = 1.0 / prob->Kappa_i(*c, c(), var_j);
//            double kappaTauInv = kappaInv / prob->Tau_i(c(), var_j);
            if (var_i != var_j) DtP_j = 0.0;

            M_c(i, j) += w * (rho * (DtV_j * V_i) + kappaInv[var_j] * (DtP_j * P_i)
                              - ((GradP_j * V_i) + (DivV_j * P_i)));

            if (var_j > 1 && var_i == var_j)
              M_c(i, j) += w * kappaTauInv[var_j] * (P_j * P_i);
          }
        }
      }
      vector<Point> c_nodal = RHS.GetDoF().GetNodalPoints(*c);
      int probDim = prob->Dim();
      int cnt = 0;
      for (int ti = 0; ti < deg_c.time; ++ti)
        for (int pi = 0; pi < probDim; ++pi)
          for (int i = 0; i < RHS.GetDoF().NodalPointsLocal(deg_c.space); ++i) {
            COMPONENT comp = prob->GetComponents()[pi];
            u_0(r_c, cnt) = prob->ut(c_nodal[cnt].CopyWithT(0.0).t(), c_nodal[cnt].CopyWithT(0.0), *c, comp);
            cnt++;
          }

      t_cell += Date() - Start_cell;
      Date Start_face;

      double z_K = sqrt(rho * prob->Kappa(*c, c()));
      for (int f = 0; f < c.Faces() - 2; ++f) {

        bool bnd = M.OnBoundary(*c, f);
        int bnd_id = bnd ? prob->BndID(c.Face(f)) : -1;
        cell cf = bnd ? c : M.find_neighbour_cell(c, f);

        SpaceTimeViscoAcousticFaceElement felem(*disc, M, c, f, prob->nL());
        DegreePair deg_cf = RHS.GetDoF().GetDegree(*cf);
        int f1 = M.find_neighbour_face_id(c.Face(f), cf);
        SpaceTimeViscoAcousticFaceElement felem_1(*disc, M, cf, f1, prob->nL());
        DGRowEntries M_cf(M, *c, *cf, false);
        row r_cf = M.find_row(cf());

        double z_Kf = sqrt(prob->Rho(*cf, cf()) * prob->Kappa(*cf, cf()));
        double alpha1 = 1.0 / (z_K + z_Kf);
        double alpha2 = z_Kf * z_K * alpha1;
        double alpha3 = z_Kf * alpha1;
        double alpha4 = z_K * alpha1;

        for (int q = 0; q < felem.nQ(); ++q) {
          int q1 = find_q_id(felem.QPoint(q), felem_1);
          double w = felem.QWeight(q);
          VectorField Nq = felem.QNormal(q);
          for (int i = 0; i < felem.i_dimension(); ++i) {
            if (felem.is_zero_test(q, i)) continue;
            VectorField V_i = felem.VelocityTestSpace(q, i);
            Scalar P_i = felem.PressureTestSpace(q, i);
            double VFlux_c_Vi = VFlux_c(V_i, Nq);
            for (int j = 0; j < felem.j_dimension() - r_c.n(); ++j) {
              if (felem.is_zero(q, j)) continue;
              VectorField V_j = felem.Velocity(q, j);
              Scalar P_j = felem.Pressure(q, j);
              double VFlux_c_Vj = VFlux_c(V_j, Nq);
              M_c(i, j) += w * (alpha1 * P_j * P_i
                                + alpha2 * VFlux_c_Vj * VFlux_c_Vi
                                + alpha3 * VFlux_c_Vj * P_i
                                + alpha4 * P_j * VFlux_c_Vi);
            }
            for (int j = 0; j < felem_1.j_dimension() - r_cf.n(); ++j) {
              if (felem_1.is_zero(q1, j)) continue;
              VectorField V1_j = felem_1.Velocity(q1, j);
              Scalar P1_j = felem_1.Pressure(q1, j);
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
      t_face += Date() - Start_face;
    }
    u_0.Accumulate();
    vout(3) << "||u_0|| = " << u_0.norm() << endl;
    Vector Bu_0 = M * u_0;
    RHS -= Bu_0;
    vout(3) << "||Mu_0|| = " << RHS.norm() << endl;


    t_cell.Max();
    t_face.Max();
    vout(1) << "      cell init assemble time " << t_cell << endl;
    vout(1) << "      face init assemble time " << t_face << endl;
  }

  M = 0;
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    SystemCell(RHS, M, c);
  }

  t_cell.Max();
  t_face.Max();
  mout.PrintInfo("TDGViscoAcousticAssemble", verbose,
                 PrintInfoEntry<Time>("cell assemble time", t_cell),
                 PrintInfoEntry<Time>("face assemble time", t_face),
                 PrintInfoEntry<Time>("assemble time", Date() - Start),
                 PrintInfoEntry<Date>("finish assemble", Date())
  );
}

void STPGViscoAcousticAssemble::SystemCell(Vector &RHS, Matrix &M, cell &c) const {
  Date Start_cell;
  SpaceTimeViscoAcousticElement elem_c(*disc, M, c, prob->nL());
  cell c_prev = M.find_previous_cell(c);
  row r_c = M.find_row(c());
  //int c_deg = elem_c.sdeg;
  DegreePair deg_c_prev = RHS.GetDoF().GetDegree(*c_prev);
  vector<vector<Scalar>> T_Matrix{calc_TMatrix(RHS.GetDoF().get_m(c()),
                                               RHS.GetDoF().get_m(c_prev()),
                                               elem_c.deg.space, deg_c_prev.space)};

  DGRowEntries M_c(M, *c, *c, false);
  DGRowEntries M_c_prev(M, *c, *c_prev, true);
  const double rho = prob->Rho(*c, c());

  if (Point temp; prob->hasPointSource(temp)) {
    Point Source = temp.CopyWithT(c().t());
    Point localPT;
    const int tmp = elem_c.SS.size() * prob->Dim();
    if (Point temp2; transformPointGlobalToLocal(Source, temp2, c)) {
      for (int q = 0; q < elem_c.nTimeQ(); ++q) {
        localPT = temp2.CopyWithT(elem_c.TimeQ.QPoint(q)[0]);
        const double time = c.min() + elem_c.TimeQ.QPoint(q)[0] * (c.max() - c.min());
        for (int i = 0; i < elem_c.i_dimension(); ++i) {
          int comp = (i % tmp) / elem_c.SS.size();
          double timeFactor = prob->F_Pkt(time, *c, comp);
          RHS(c(), i) += elem_c.time_qWeight[q] * timeFactor * elem_c.Evaluate(localPT, i);
        }
      }
    }
  }

//  double kappaInv[2 + prob->nL()];
//  double tau[2 + prob->nL()];
//  for (int j = 0; j < 2 + prob->nL(); j++) {
//    kappaInv[j] = 1.0 / prob->Kappa(*c, c(), j);
//    tau[j] = prob->Tau_i(c(), j);
//  }
  vector<double> kappaInv(2 + prob->nL());
  vector<double> tau(2 + prob->nL());
  kappaInv[1] = 1.0 / prob->Kappa_i(*c, c(), 0);
  for (int j = 2; j < 2 + prob->nL(); j++) {
    kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j - 1);
    tau[j] = prob->Tau_i(*c, c(), j - 1);
  }
  for (int i = 0; i < elem_c.i_dimension(); ++i) {
    int var_i = elem_c.variable(i);
    for (int q = 0; q < elem_c.nQ(); ++q) {
      double w = elem_c.QWeight(q);
      Point QP = elem_c.QPoint(q);
      double t = QP.t();
      VectorField V_i = elem_c.VelocityTestSpace(q, i);
      double P_i = elem_c.PressureTestSpace(q, i);

      if (prob->HasRHS()) {
        VectorField V_rhs = prob->F_Velocity(t, *c, QP);
        if (var_i == 0) {
          RHS(c(), i) += w * (V_rhs * V_i);
        } else if (var_i == 1) {
          Scalar P_rhs = prob->F_Pressure(t, *c, QP);
          RHS(c(), i) += w * (P_rhs * P_i);
        } else {
          DampingVector DV = prob->F_DampingPressure(t, *c, QP);
          RHS(c(), i) += w * (DV[var_i - 2] * P_i);
        }
      }

      if (c_prev != c) {
        for (int j = 0; j < elem_c.j_dimension() - r_c.n(); ++j) {
          int var_j = elem_c.variable(j);
          VectorField DtV_j = elem_c.DtVelocity(q, j);
          Scalar DivV_j = elem_c.DivVelocity(q, j);
          Scalar DtP_j = 0.0;
          if (var_i == var_j) DtP_j = elem_c.DtPressure(q, j);
          VectorField GradP_j = elem_c.GradPressure(q, j);

          double value = w * ((rho * (DtV_j * V_i)
                               + kappaInv[var_j] * (DtP_j * P_i))
                              - ((GradP_j * V_i) + (DivV_j * P_i)));

          if (elem_c.deg.space == deg_c_prev.space) {
            M_c_prev(i, j) += value;
          } else {
            for (int k = 0; k < T_Matrix[j].size(); ++k) {
              M_c_prev(i, k) += value * T_Matrix[j][k];
            }
          }

          if (var_j > 1 && var_i == var_j) {
            Scalar P_j = elem_c.Pressure(q, j);
            double kappaTauInv = kappaInv[var_j] / tau[var_j];
            double valueInner = w * kappaTauInv * (P_j * P_i);
            if (elem_c.deg.space == deg_c_prev.space) {
              M_c_prev(i, j) += valueInner;
            } else {
              for (int k = 0; k < T_Matrix[j].size(); ++k) {
                M_c_prev(i, k) += valueInner * T_Matrix[j][k];
              }
            }
          }
        }
      }
      for (int j = elem_c.j_dimension() - r_c.n(); j < elem_c.j_dimension(); ++j) {
        int var_j = elem_c.variable(j);
        VectorField DtV_j = elem_c.DtVelocity(q, j);
        Scalar DivV_j = elem_c.DivVelocity(q, j);
        Scalar DtP_j = 0.0;
        if (var_i == var_j) DtP_j = elem_c.DtPressure(q, j);
        VectorField GradP_j = elem_c.GradPressure(q, j);
        int j_idx = j - elem_c.j_dimension() + r_c.n();
        M_c(i, j_idx) += w * ((rho * (DtV_j * V_i)
                               + kappaInv[var_j] * (DtP_j * P_i))
                              - ((GradP_j * V_i) + (DivV_j * P_i)));

        if (var_j > 1 && var_i == var_j) {
          Scalar P_j = elem_c.Pressure(q, j);
          double kappaTauInv = kappaInv[var_j] / tau[var_j];
          M_c(i, j_idx) += w * kappaTauInv * (P_j * P_i);
        }
      }
    }
  }

  t_cell += Date() - Start_cell;
  Date Start_face;

  const double z_K = sqrt(prob->Rho(*c, c()) * prob->Kappa(*c, c()));
  for (int f = 0; f < c.Faces() - 2; ++f) {
    bool bnd = M.OnBoundary(*c, f);
    int bnd_id = bnd ? prob->BndID(c.Face(f)) : -1;
    cell cf = bnd ? c : M.find_neighbour_cell(c, f);
    SpaceTimeViscoAcousticFaceElement felem_c(*disc, M, c, f, prob->nL());
    int f1 = M.find_neighbour_face_id(c.Face(f), cf);
    SpaceTimeViscoAcousticFaceElement felem_cf(*disc, M, cf, f1, prob->nL());
    cell cf_prev = M.find_previous_cell(cf);

    if (cf_prev == M.cells_end()) {
      THROW("cell not found")
    }

    DegreePair deg_cf_prev = RHS.GetDoF().GetDegree(*cf_prev);

    vector<vector<Scalar>> T_Matrix_f_prev{calc_TMatrix(RHS.GetDoF().get_m(cf()),
                                                        RHS.GetDoF().get_m(cf_prev()),
                                                        felem_cf.deg.space,
                                                        deg_cf_prev.space)};

    DGRowEntries M_cf(M, *c, *cf, false);
    row r_cf = M.find_row(cf());

    const double z_Kf = sqrt(prob->Rho(*cf, cf()) * prob->Kappa(*cf, cf()));
    const double alpha1 = 1.0 / (z_K + z_Kf);
    const double alpha2 = z_Kf * z_K * alpha1;
    const double alpha3 = z_Kf * alpha1;
    const double alpha4 = z_K * alpha1;
    for (int q = 0; q < felem_c.nQ(); ++q) {
      int q1 = find_q_id(felem_c.QPoint(q), felem_cf);
      double w = felem_c.QWeight(q);
      VectorField Nq = felem_c.QNormal(q);
      Point QP = felem_c.QPoint(q);

      Scalar VN_bnd = Nq * prob->ut_Velocity(QP.t(), QP, *c);
      Scalar P_bnd = prob->ut_Pressure(QP.t(), QP, *c) + prob->ut_DampingPressure(QP.t(), QP, *c).sum();

      for (int i = 0; i < felem_c.i_dimension(); ++i) {
        if (felem_c.is_zero_test(q, i)) continue;
        VectorField V_i = felem_c.VelocityTestSpace(q, i);
        double P_i = felem_c.PressureTestSpace(q, i);
        double VFlux_c_Vi = VFlux_c(V_i, Nq);

        if (bnd_id == 1) {
          RHS(c(), i) += w * (1 / z_K * P_i + VFlux_c_Vi) * P_bnd;
        } else if (bnd_id == 2 || bnd_id == 0) {
          RHS(c(), i) += w * (z_K * VFlux_c_Vi + P_i) * VN_bnd;
        }

        for (int j = 0; j < felem_c.j_dimension() - r_c.n(); ++j) {
          if (c_prev == c) break;
          if (felem_c.is_zero(q, j)) continue;
          VectorField V_j = felem_c.Velocity(q, j);
          double P_j = felem_c.Pressure(q, j);
          double VFlux_c_Vj = VFlux_c(V_j, Nq);
          double value = w * (alpha1 * P_j * P_i
                              + alpha2 * VFlux_c_Vj * VFlux_c_Vi
                              + alpha3 * VFlux_c_Vj * P_i
                              + alpha4 * P_j * VFlux_c_Vi);

          if (elem_c.deg.space == deg_c_prev.space) {
            M_c_prev(i, j) += value;
          } else {
            for (int k = 0; k < T_Matrix[j].size(); ++k) {
              M_c_prev(i, k) += value * T_Matrix[j][k];
            }
          }
        }
        for (int j = felem_c.j_dimension() - r_c.n(); j < felem_c.j_dimension(); ++j) {
          if (felem_c.is_zero(q, j)) continue;
          int idx_j = j - (felem_c.j_dimension() - r_c.n());
          VectorField V_j = felem_c.Velocity(q, j);
          double P_j = felem_c.Pressure(q, j);
          double VFlux_c_Vj = VFlux_c(V_j, Nq);
          M_c(i, idx_j) += w * (alpha1 * P_j * P_i
                                + alpha2 * VFlux_c_Vj * VFlux_c_Vi
                                + alpha3 * VFlux_c_Vj * P_i
                                + alpha4 * P_j * VFlux_c_Vi);
        }
        for (int j = 0; j < felem_cf.j_dimension() - r_cf.n(); ++j) {
          if (cf_prev == cf || c == c_prev) break;
          DGRowEntries M_cf_prev(M, *c, *cf_prev, true);

          if (felem_cf.is_zero(q1, j)) continue;
          VectorField V1_j = felem_cf.Velocity(q1, j);
          double P1_j = felem_cf.Pressure(q1, j);

          double PFlux_cf_P1j = PFlux_cf(P1_j, bnd, bnd_id);
          double VFlux_cf_V1j = VFlux_cf(V1_j, Nq, bnd, bnd_id);
          double value = w * (alpha1 * PFlux_cf_P1j * P_i
                              + alpha2 * VFlux_cf_V1j * VFlux_c_Vi
                              + alpha3 * VFlux_cf_V1j * P_i
                              + alpha4 * PFlux_cf_P1j * VFlux_c_Vi);

          if (felem_cf.deg.space == deg_cf_prev.space) {
            M_cf_prev(i, j) -= value;
          } else {
            for (int k = 0; k < T_Matrix_f_prev[j].size(); ++k)
              M_cf_prev(i, k) -= value * T_Matrix_f_prev[j][k];
          }
        }
        for (int j = felem_cf.j_dimension() - r_cf.n(); j < felem_cf.j_dimension(); ++j) {
          if (felem_cf.is_zero(q1, j)) continue;
          VectorField V1_j = felem_cf.Velocity(q1, j);
          double P1_j = felem_cf.Pressure(q1, j);

          int idx_j = j - (felem_cf.j_dimension() - r_cf.n());
          double PFlux_cf_P1j = PFlux_cf(P1_j, bnd, bnd_id);
          double VFlux_cf_V1j = VFlux_cf(V1_j, Nq, bnd, bnd_id);

          M_cf(i, idx_j) -= w * (alpha1 * PFlux_cf_P1j * P_i
                                 + alpha2 * VFlux_cf_V1j * VFlux_c_Vi
                                 + alpha3 * VFlux_cf_V1j * P_i
                                 + alpha4 * PFlux_cf_P1j * VFlux_c_Vi);
        }
      }
    }
  }
  t_face += Date() - Start_face;
}

void STPGViscoAcousticAssemble::SystemAddDoubleD(Matrix &M) const {
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    Date Start_cell;
    SpaceTimeViscoAcousticElement elem_c(*disc, M, c, prob->nL());
    cell c_prev = M.find_previous_cell(c);
    row r_c = M.find_row(c());
    DegreePair deg_c_prev = M.GetDoF().GetDegree(*c_prev);
    vector<vector<Scalar> > T_Matrix(M.GetDoF().get_m(c()));
    vector<Scalar> vc_diag(M.GetDoF().get_m(c()));
    vector<Scalar> vc(M.GetDoF().get_m(c_prev()));
    if (elem_c.deg.space != deg_c_prev.space) {
      for (int i = 0; i < vc_diag.size(); ++i)
        T_Matrix[i].resize(vc.size());
      for (int i = 0; i < vc.size(); ++i) {
        for (int k = 0; k < vc.size(); ++k)
          vc[k] = 0.0;
        vc[i] = 1.;
        PX_P0(vc, vc_diag, deg_c_prev.space, elem_c.deg.space);
        for (int k = 0; k < vc_diag.size(); ++k)
          T_Matrix[k][i] = vc_diag[k];
      }
    }
    DGRowEntries M_c(M, *c, *c, false);
    DGRowEntries M_c_prev(M, *c, *c_prev, true);
    const double rho = prob->Rho(*c, c());
    vector<double> kappaInv(2 + prob->nL());
    vector<double> kappaTauInv(2 + prob->nL());
    kappaInv[1] = 1.0 / prob->Kappa_i(*c, c(), 0);
    for (int j = 2; j < 2 + prob->nL(); j++) {
      kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j - 1);
      kappaTauInv[j] = kappaInv[j] / prob->Tau_i(*c, c(), j - 1);
    }

    for (int q = 0; q < elem_c.nQ(); ++q) {
      double w = elem_c.QWeight(q);
      Point QP = elem_c.QPoint(q);

      for (int i = 0; i < elem_c.i_dimension(); ++i) {
        int var_i = elem_c.variable(i);
        if (var_i < 2) continue;
        double P_i = elem_c.PressureTestSpace(q, i);
        for (int j = 0; j < elem_c.j_dimension() - r_c.n(); ++j) {
          if (c_prev == c) break;
          int var_j = elem_c.variable(j);
          if (var_i != var_j) continue;
          Scalar P_j = elem_c.Pressure(q, j);
//          double kappaTauInv =
//              1.0 / (prob->Kappa(c(), var_j) * prob->Tau_i(c(), var_j));
          if (elem_c.deg.space == deg_c_prev.space) {
            M_c_prev(i, j) -= 2 * w * kappaTauInv[var_j] * (P_j * P_i);
          } else {
            for (int k = 0; k < vc.size(); ++k) {
              M_c_prev(i, k) += 2 * w * kappaTauInv[var_j] * (P_j * P_i) * T_Matrix[j][k];
            }
          }
        }
        for (int j = elem_c.j_dimension() - r_c.n(); j < elem_c.j_dimension(); ++j) {
          int var_j = elem_c.variable(j);
          if (var_i != var_j) continue;
          Scalar P_j = elem_c.Pressure(q, j);
//          double kappaTauInv = 1.0 / (prob->Kappa(c(), var_j) * prob->Tau_i(c(), var_j));
          M_c(i, j - elem_c.j_dimension() + r_c.n()) += 2 * w * kappaTauInv[var_j] * (P_j * P_i);
        }
      }
    }
  }
}

std::pair<double, double>
STPGViscoAcousticAssemble::DiscNorm(const Matrix &L, const Vector &u) const {
  Matrix M(u);
  M = 0;
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    SpaceTimeViscoAcousticElement elem_c(*disc, M, c, prob->nL());
    cell c_prev = M.find_previous_cell(c);
    row r_c = M.find_row(c());
    DegreePair deg_c_prev = M.GetDoF().GetDegree(*c_prev);
    vector<vector<Scalar> > T_Matrix(M.GetDoF().get_m(c()));
    vector<Scalar> vc_diag(M.GetDoF().get_m(c()));
    vector<Scalar> vc(M.GetDoF().get_m(c_prev()));
    if (elem_c.deg.space != deg_c_prev.space) {
      for (int i = 0; i < vc_diag.size(); ++i)
        T_Matrix[i].resize(vc.size());
      for (int i = 0; i < vc.size(); ++i) {
        for (int k = 0; k < vc.size(); ++k)
          vc[k] = 0.0;
        vc[i] = 1.;
        PX_P0(vc, vc_diag, deg_c_prev.space, elem_c.deg.space);
        for (int k = 0; k < vc_diag.size(); ++k)
          T_Matrix[k][i] = vc_diag[k];
      }
    }
    DGRowEntries M_c(M, *c, *c, false);
    DGRowEntries M_c_prev(M, *c, *c_prev, true);
    double rho = prob->Rho(*c, c());
    vector<double> kappaInv(2 + prob->nL());
    vector<double> kappaTauInv(2 + prob->nL());
    kappaInv[1] = 1.0 / prob->Kappa_i(*c, c(), 0);
    for (int j = 2; j < 2 + prob->nL(); j++) {
      kappaInv[j] = 1.0 / prob->Kappa_i(*c, c(), j - 1);
      kappaTauInv[j] = kappaInv[j] / prob->Tau_i(*c, c(), j - 1);
    }
    for (int q = 0; q < elem_c.nQ(); ++q) {
      double w = elem_c.QWeight(q);
      Point QP = elem_c.QPoint(q);
      for (int i = 0; i < elem_c.i_dimension(); ++i) {
        int var_i = elem_c.variable(i);
        VectorField V_i = elem_c.VelocityTestSpace(q, i);
        double P_i = elem_c.PressureTestSpace(q, i);
        for (int j = 0; j < elem_c.j_dimension() - r_c.n(); ++j) {
          if (c_prev == c) break;
          int var_j = elem_c.variable(j);
          if (var_j != var_i) continue;
          VectorField V_j = elem_c.Velocity(q, j);
          Scalar P_j = elem_c.Pressure(q, j);
//          double kappaInv = 1.0 / prob->Kappa(c(), var_j);
          if (elem_c.deg.space == deg_c_prev.space) {
            M_c_prev(i, j) += w * (rho * (V_j * V_i) + kappaInv[var_j] * (P_j * P_i));
          } else {
            for (int k = 0; k < vc.size(); ++k) {
              M_c_prev(i, k) += w * (rho * (V_j * V_i) + kappaInv[var_j] * (P_j * P_i)) * T_Matrix[j][k];
            }
          }
        }
        for (int j = elem_c.j_dimension() - r_c.n(); j < elem_c.j_dimension(); ++j) {
          int var_j = elem_c.variable(j);
          if (var_j != var_i) continue;
          VectorField V_j = elem_c.Velocity(q, j);
          Scalar P_j = elem_c.Pressure(q, j);
//          double kappaInv = 1.0 / prob->Kappa(c(), var_j);
          int idx_j = j - elem_c.j_dimension() + r_c.n();
          M_c(i, idx_j) += w * (rho * (V_j * V_i) + kappaInv[var_j] * (P_j * P_i));
        }
      }
    }
  }

  Vector Mu(M * u);
  Vector Lu(L * u);

  Preconditioner *PC = GetPC("PointBlockGaussSeidel");
  LinearSolver *S = GetLinearSolver(PC);
  (*S)(M);

  Vector mLu(u);
  mLu = (*S) * Lu;

  double discNormW = u * Mu;
  double discNormV = Lu * mLu;

  mout << "discrete norm ||u-U||_W  = " << scientific //<< setprecision(10)
       << sqrt(discNormW) << endl;
  mout << "discrete norm ||u-U||_V  = " << scientific //<< setprecision(10)
       << sqrt(discNormW + discNormV) << endl;

  return pair<double, double>(std::sqrt(discNormW), std::sqrt(discNormW + discNormV));
}

double STPGViscoAcousticAssemble::L2Error(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    SpaceTimeViscoAcousticElement elem(*disc, u, c, prob->nL());
    row r = u.find_row(c());

    cell c_prev = c;
    vector<Scalar> c_nodal_value;

    DegreePair deg = u.GetDoF().GetDegree(*c);

    if (c.min() != 0) {
      face ff = u.find_face(c.Face(c.Faces() - 2));
      c_prev = u.find_cell_or_overlap_cell(ff.Left());
    } else {
      fillNodalValuesAtStart(deg.space, c, c_nodal_value, u);
    }

    row r_prev = u.find_row(c_prev());

    DegreePair deg_prev = u.GetDoF().GetDegree(*c_prev);
    int shift = (deg_prev.time - 1) * (r_prev.n() / deg_prev.time);

    const Scalar *u_r_prev = u(r_prev);
    const Scalar *u_r = u(r);

    vector<vector<Scalar>> T_Matrix{calc_TMatrix(u.GetDoF().get_m(c()),
                                                 u.GetDoF().get_m(c_prev()),
                                                 deg.space, deg_prev.space)};
    const int dim = c.dim();
    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField V = zero;
      vector<double> P(1 + prob->nL(), 0.0);

      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      double t = QP.t();
      for (int j = 0; j < elem.j_dimension() - r.n(); ++j) {
        int var_j = elem.variable(j);
        VectorField phi_V_j = elem.Velocity(q, j);
        double phi_P_j = elem.Pressure(q, j);

        if (deg_prev.space != deg.space) {
          for (int k = 0; k < T_Matrix[j].size(); ++k) {
            double value = u_r_prev[k + shift] * T_Matrix[j][k];
            V += value * phi_V_j;
            if (var_j > 0) {
              P[var_j - 1] += value * phi_P_j;
            }
          }
        } else {
          if (c.min() != 0.0) {
            V += u_r_prev[j + shift] * phi_V_j;
            if (var_j > 0)
              P[var_j - 1] += u_r_prev[j + shift] * phi_P_j;
          } else {
            V += c_nodal_value[j] * phi_V_j;
            if (var_j > 0)
              P[var_j - 1] += c_nodal_value[j] * phi_P_j;
          }
        }
      }
      for (int j = elem.j_dimension() - r.n(); j < elem.j_dimension(); ++j) {
        int var_j = elem.variable(j);
        VectorField phi_V_j = elem.Velocity(q, j);
        double phi_P_j = elem.Pressure(q, j);
        V += u_r[j - elem.j_dimension() + r.n()] * phi_V_j;
        if (var_j > 0)
          P[var_j - 1] += u_r[j - elem.j_dimension() + r.n()] * phi_P_j;
      }

      for (int k = 0; k < 1 + prob->nL(); ++k) {
        COMPONENT comp = prob->GetComponents()[dim + k];
        nrm += w * std::pow(P[k] - prob->ut(t, QP, *c, comp), 2);
      }
      V -= prob->ut_Velocity(QP.t(), QP, *c);
      nrm += w * V * V;
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
};

double STPGViscoAcousticAssemble::L1Error(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    SpaceTimeViscoAcousticElement elem(*disc, u, c, prob->nL());
    row r = u.find_row(c());

    cell c_prev = c;
    vector<Scalar> c_nodal_value;

    DegreePair deg = u.GetDoF().GetDegree(*c);

    if (c.min() != 0) {
      face ff = u.find_face(c.Face(c.Faces() - 2));
      c_prev = u.find_cell_or_overlap_cell(ff.Left());
    } else {
      fillNodalValuesAtStart(deg.space, c, c_nodal_value, u);
    }

    row r_prev = u.find_row(c_prev());

    DegreePair deg_prev = u.GetDoF().GetDegree(*c_prev);
    int shift = (deg_prev.time - 1) * (r_prev.n() / deg_prev.time);

    const Scalar *u_r_prev = u(r_prev);
    const Scalar *u_r = u(r);

    vector<vector<Scalar>> T_Matrix{calc_TMatrix(u.GetDoF().get_m(c()),
                                                 u.GetDoF().get_m(c_prev()),
                                                 deg.space, deg_prev.space)};
    const int dim = c.dim();
    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField V = zero;
      vector<double> P(1 + prob->nL(), 0.0);

      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);

      for (int j = 0; j < elem.j_dimension() - r.n(); ++j) {
        int var_j = elem.variable(j);
        VectorField phi_V_j = elem.Velocity(q, j);
        double phi_P_j = elem.Pressure(q, j);

        if (deg_prev.space != deg.space) {
          for (int k = 0; k < T_Matrix[j].size(); ++k) {
            double value = u_r_prev[k + shift] * T_Matrix[j][k];
            V += value * phi_V_j;
            if (var_j > 0) {
              P[var_j - 1] += value * phi_P_j;
            }
          }
        } else {
          if (c.min() != 0.0) {
            V += u_r_prev[j + shift] * phi_V_j;
            if (var_j > 0)
              P[var_j - 1] += u_r_prev[j + shift] * phi_P_j;
          } else {
            V += c_nodal_value[j] * phi_V_j;
            if (var_j > 0)
              P[var_j - 1] += c_nodal_value[j] * phi_P_j;
          }
        }
      }
      for (int j = elem.j_dimension() - r.n(); j < elem.j_dimension(); ++j) {
        int var_j = elem.variable(j);
        VectorField phi_V_j = elem.Velocity(q, j);
        double phi_P_j = elem.Pressure(q, j);
        V += u_r[j - elem.j_dimension() + r.n()] * phi_V_j;
        if (var_j > 0)
          P[var_j - 1] += u_r[j - elem.j_dimension() + r.n()] * phi_P_j;
      }

      for (int k = 0; k < 1 + prob->nL(); ++k) {
        COMPONENT comp = prob->GetComponents()[dim + k];
        nrm += w * std::abs(P[k] - prob->ut(QP.t(), QP, *c, comp));
      }
      nrm += w * absNorm(V - prob->ut_Velocity(QP.t(), QP, *c));

    }
  }
  return PPM->SumOnCommSplit(nrm, u.CommSplit());
}

double STPGViscoAcousticAssemble::LInfError(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    SpaceTimeViscoAcousticElement elem(*disc, u, c, prob->nL());
    row r = u.find_row(c());

    cell c_prev = c;
    vector<Scalar> c_nodal_value;

    DegreePair deg = u.GetDoF().GetDegree(*c);

    if (c.min() != 0) {
      face ff = u.find_face(c.Face(c.Faces() - 2));
      c_prev = u.find_cell_or_overlap_cell(ff.Left());
    } else {
      fillNodalValuesAtStart(deg.space, c, c_nodal_value, u);
    }

    row r_prev = u.find_row(c_prev());

    DegreePair deg_prev = u.GetDoF().GetDegree(*c_prev);
    int shift = (deg_prev.time - 1) * (r_prev.n() / deg_prev.time);

    const Scalar *u_r_prev = u(r_prev);
    const Scalar *u_r = u(r);

    vector<vector<Scalar>> T_Matrix{calc_TMatrix(u.GetDoF().get_m(c()),
                                                 u.GetDoF().get_m(c_prev()),
                                                 deg.space, deg_prev.space)};
    const int dim = c.dim();
    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField V = zero;
      vector<double> P(1 + prob->nL(), 0.0);

      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      double t = QP.t();
      for (int j = 0; j < elem.j_dimension() - r.n(); ++j) {
        int var_j = elem.variable(j);
        VectorField phi_V_j = elem.Velocity(q, j);
        double phi_P_j = elem.Pressure(q, j);

        if (deg_prev.space != deg.space) {
          for (int k = 0; k < T_Matrix[j].size(); ++k) {
            double value = u_r_prev[k + shift] * T_Matrix[j][k];
            V += value * phi_V_j;
            if (var_j > 0) {
              P[var_j - 1] += value * phi_P_j;
            }
          }
        } else {
          if (c.min() != 0.0) {
            V += u_r_prev[j + shift] * phi_V_j;
            if (var_j > 0)
              P[var_j - 1] += u_r_prev[j + shift] * phi_P_j;
          } else {
            V += c_nodal_value[j] * phi_V_j;
            if (var_j > 0)
              P[var_j - 1] += c_nodal_value[j] * phi_P_j;
          }
        }
      }
      for (int j = elem.j_dimension() - r.n(); j < elem.j_dimension(); ++j) {
        int var_j = elem.variable(j);
        VectorField phi_V_j = elem.Velocity(q, j);
        double phi_P_j = elem.Pressure(q, j);
        V += u_r[j - elem.j_dimension() + r.n()] * phi_V_j;
        if (var_j > 0)
          P[var_j - 1] += u_r[j - elem.j_dimension() + r.n()] * phi_P_j;
      }

      for (int k = 0; k < 1 + prob->nL(); ++k) {
        COMPONENT comp = prob->GetComponents()[dim + k];
        nrm = std::max(nrm, std::abs(P[k] - prob->ut(t, QP, *c, comp)));
      }
      nrm = maxNorm(V - prob->ut_Velocity(QP.t(), QP, *c));
    }
  }
  return PPM->Max(nrm, u.CommSplit());
};

void STPGViscoAcousticAssemble::fillNodalValuesAtStart(int deg, cell &c,
                                                       vector<Scalar> &c_nodal_value,
                                                       const VectorMatrixBase &m) const {
  vector<Point> c_nodal = m.GetDoF().GetNodalPoints(*c);
  c_nodal_value.resize(m.GetDoF().get_m(c()));
  int cnt = 0;
  for (int pi = 0; pi < prob->Dim(); ++pi)
    for (int i = 0; i < m.GetDoF().NodalPointsLocal(deg); ++i) {
      COMPONENT comp = prob->GetComponents()[pi];
      c_nodal_value[cnt] = prob->ut(c_nodal[cnt].CopyWithT(0.0).t(), c_nodal[cnt].CopyWithT(0.0), *c, comp);
      cnt++;
    }
}

double STPGViscoAcousticAssemble::GNError(const Vector &) const {
  mout << "WARNING ... GNError not implemented!" << endl;
  return 0.0;
}

double STPGViscoAcousticAssemble::L2Norm(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    SpaceTimeViscoAcousticElement elem(*disc, u, c, prob->nL());
    row r = u.find_row(c());

    cell c_prev = c;
    vector<Scalar> c_nodal_value;

    DegreePair deg = u.GetDoF().GetDegree(*c);

    if (c.min() != 0) {
      face ff = u.find_face(c.Face(c.Faces() - 2));
      c_prev = u.find_cell_or_overlap_cell(ff.Left());
    } else {
      fillNodalValuesAtStart(deg.space, c, c_nodal_value, u);
    }

    row r_prev = u.find_row(c_prev());

    DegreePair deg_prev = u.GetDoF().GetDegree(*c_prev);
    int shift = (deg_prev.time - 1) * (r_prev.n() / deg_prev.time);

    const Scalar *u_r_prev = u(r_prev);
    const Scalar *u_r = u(r);

    vector<vector<Scalar>> T_Matrix{calc_TMatrix(u.GetDoF().get_m(c()),
                                                 u.GetDoF().get_m(c_prev()),
                                                 deg.space, deg_prev.space)};
    const int dim = c.dim();
    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField V = zero;
      vector<double> P(1 + prob->nL(), 0.0);

      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);

      for (int j = 0; j < elem.j_dimension() - r.n(); ++j) {
        int var_j = elem.variable(j);
        VectorField phi_V_j = elem.Velocity(q, j);
        double phi_P_j = elem.Pressure(q, j);

        if (deg_prev.space != deg.space) {
          for (int k = 0; k < T_Matrix[j].size(); ++k) {
            double value = u_r_prev[k + shift] * T_Matrix[j][k];
            V += value * phi_V_j;
            if (var_j > 0) {
              P[var_j - 1] += value * phi_P_j;
            }
          }
        } else {
          if (c.min() != 0.0) {
            V += u_r_prev[j + shift] * phi_V_j;
            if (var_j > 0)
              P[var_j - 1] += u_r_prev[j + shift] * phi_P_j;
          } else {
            V += c_nodal_value[j] * phi_V_j;
            if (var_j > 0)
              P[var_j - 1] += c_nodal_value[j] * phi_P_j;
          }
        }
      }
      for (int j = elem.j_dimension() - r.n(); j < elem.j_dimension(); ++j) {
        int var_j = elem.variable(j);
        VectorField phi_V_j = elem.Velocity(q, j);
        double phi_P_j = elem.Pressure(q, j);
        V += u_r[j - elem.j_dimension() + r.n()] * phi_V_j;
        if (var_j > 0)
          P[var_j - 1] += u_r[j - elem.j_dimension() + r.n()] * phi_P_j;
      }

      for (int k = 0; k < 1 + prob->nL(); ++k)
        nrm += w * std::pow(P[k], 2);
      nrm += w * V * V;
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

void STPGViscoAcousticAssemble::DualRHS_linear(Vector &RHS,
                                               const Vector &U,
                                               Point a,
                                               Point b) const {

  RHS = 0;

  if (a.t() != b.t()) { // volume evaluation
    Exit("TODO")
  } else { // evaluation at some time point
    for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {
      if (c.max() != a.t()) continue;

      SpaceTimeViscoAcousticElement elem(*disc, RHS, c, prob->nL());

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
        Scalar scale = 0.0;
        const int start_idx = elem.deg.time * (elem.j_dimension() / (elem.deg.time + 1));
        for (int j = start_idx; j < elem.j_dimension(); ++j)
          scale += elem.Velocity(q, j)[0];

        for (int j = start_idx; j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          VectorField V_j = elem.Velocity(q, j) / scale;
          Scalar P_j = elem.Pressure(q, j) / scale;
          RHS(r_RHS)[j - elem.j_dimension() + r_RHS.n()] -= w * (j_P[var_j] * P_j + j_V * V_j);
        }
      }
    }
  }
}

void STPGViscoAcousticAssemble::DualRHS_quadratic(Vector &RHS,
                                                  const Vector &U,
                                                  Point a,
                                                  Point b) const {
  RHS = 0;

  if (a.t() != b.t()) { // volume evaluation
    Exit("CHECK");
    for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {
      face f_diag = RHS.find_face(c.Face(c.Faces() - 2));
      cell c_prev = RHS.find_cell(f_diag.Left());
      if (c_prev == RHS.cells_end())
        c_prev = RHS.find_overlap_cell(f_diag.Left());
      if (c_prev == RHS.overlap_end())
        c_prev = c;
      SpaceTimeViscoAcousticElement elem(*disc, RHS, c, prob->nL());
      row r_RHS = RHS.find_row(c());
      row r_RHS_prev = RHS.find_row(c_prev());
      for (int q = 0; q < elem.nQ(); ++q) {
        double w = elem.QWeight(q);
        Point QP = elem.QPoint(q);

        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;
        if (a.t() > QP.t() || QP.t() > b.t()) continue;

        for (int i = 0; i < elem.j_dimension(); ++i) {
          //int var_i = elem.variable(i);
          Scalar P_i = elem.Pressure(q, i);
          VectorField V_i = elem.Velocity(q, i);
          for (int j = 0; j < elem.j_dimension(); ++j) {
            //int var_j = elem.variable(j);
            VectorField V_j = elem.Velocity(q, j);
            Scalar P_j = elem.Pressure(q, j);
            if (j < elem.j_dimension() - r_RHS.n()) {
              if (c_prev != c) {
                Exit("what if space-deg are different?");
                RHS(r_RHS_prev)[2 * r_RHS.n() - elem.j_dimension() + j]
                    -=
                    w * U(r_RHS_prev)[i + 2 * r_RHS.n()
                                      - elem.j_dimension()] * (P_i * P_j + V_i * V_j);
              }
            } else {
              RHS(r_RHS)[j - elem.j_dimension() + r_RHS.n()]
                  -= w * U(r_RHS)[i - elem.j_dimension() + r_RHS.n()]
                     * (P_i * P_j + V_i * V_j);
            }
          }
        }
      }
    }
  } else { // evaluation at some time point
    for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {
      if (c.max() != a.t()) continue;
      SpaceTimeViscoAcousticElement elem(*disc, RHS, c, prob->nL());
      row r_RHS = RHS.find_row(c());
      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        double w = elem.SpaceQWeight(q);
        Point QP = elem.QPoint(q).CopyWithT(0.0);

        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;

        const int start_idx = elem.deg.time * (elem.j_dimension() / (elem.deg.time + 1));

        Scalar scale = 0.0;
        for (int j = start_idx; j < elem.j_dimension(); ++j)
          scale += elem.Velocity(q, j)[0];

        for (int i = start_idx; i < elem.j_dimension(); ++i) {
          VectorField V_i = elem.Velocity(q, i) / scale;
          Scalar P_i = elem.Pressure(q, i) / scale;

          for (int j = start_idx; j < elem.j_dimension(); ++j) {
            Scalar P_j = elem.Pressure(q, j) / scale;
            VectorField V_j = elem.Velocity(q, j) / scale;
            const int idx_j = j - elem.j_dimension() + r_RHS.n();
            const int idx_i = i - elem.j_dimension() + r_RHS.n();
            RHS(r_RHS)[idx_j] -= w * U(r_RHS)[idx_i] * (P_i * P_j + V_i * V_j);
          }
        }
      }
    }
  }
}

double STPGViscoAcousticAssemble::DualErrorEstimateCell(const cell &c,
                                                        const Vector &u,
                                                        const Vector &u_star) const {
  Exit("TODO")
  double eta = 0.0;

  SpaceTimeViscoAcousticElement elem(*disc, u, c, prob->nL());
  row r = u.find_row(c());

  cell c_prev = c;
  vector<Scalar> c_nodal_value;

  face f = u.find_face(c.Face(c.Faces() - 2));
  c_prev = u.find_cell(f.Left());
  if (c_prev == u.cells_end())
    c_prev = u.find_overlap_cell(f.Left());
  row r_prev = u.find_row(c_prev());

  const Scalar *u_r_prev = u(r_prev);
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

  double rho_2_V = 0.0;
  double rho_2_P = 0.0;
  double omega_2_V = 0.0;
  double omega_2_P = 0.0;

  DegreePair deg = u.GetDoF().GetDegree(*c);

  if (c.min() != 0) {
    c_prev = u.find_cell_or_overlap_cell(f.Left());
  } else {
    vector<Point> c_nodal = u.GetDoF().GetNodalPoints(*c);
    c_nodal_value.resize(c_nodal.size());
    int cnt = 0;
    for (int ti = 0; ti < deg.time; ++ti) {
      for (int pi = 0; pi < prob->Dim(); ++pi) {
        COMPONENT comp = prob->GetComponents()[pi];
        for (int i = 0; i < u.GetDoF().NodalPointsLocal(deg.space); ++i) {
          c_nodal_value[cnt] = prob->ut(c_nodal[cnt].CopyWithT(0.0).t(), c_nodal[cnt].CopyWithT(0.0), *c, comp);
          cnt++;
        }
      }
    }
  }

  DegreePair deg_prev = u.GetDoF().GetDegree(*c_prev);
  int shift = (deg_prev.time - 1) * (r_prev.n() / deg_prev.time);

  vector<Scalar> vc(u.GetDoF().get_m(c()));
  vector<Scalar> vc_prev(u.GetDoF().get_m(c_prev()));
  vector<vector<Scalar> > T_Matrix(vc.size());

  if (deg.space != deg_prev.space) {
    for (int i = 0; i < vc.size(); ++i)
      T_Matrix[i].resize(vc_prev.size());
    for (int i = 0; i < vc_prev.size(); ++i) {
      for (int k = 0; k < vc_prev.size(); ++k)
        vc_prev[k] = 0.0;
      vc_prev[i] = 1.;
      PX_P0(vc_prev, vc, deg_prev.space, deg.space);
      for (int k = 0; k < vc.size(); ++k)
        T_Matrix[k][i] = vc[k];
    }
  }

  double rho = prob->Rho(*c, c());
  double kappaInv = 1.0 / prob->Kappa(*c, c());

  for (int q = 0; q < elem.nQ(); ++q) {
    VectorField V = 0.0;
    double P = 0.0;

    double w = elem.QWeight(q);
    Point QP = elem.QPoint(q);

    for (int j = 0; j < elem.j_dimension() - r.n(); ++j) {
      VectorField phi_V_j = -1.0 * rho * elem.DtVelocity(q, j) + elem.GradPressure(q, j);
      double phi_P_j = -1.0 * kappaInv * elem.DtPressure(q, j) + elem.DivVelocity(q, j);

      if (deg_prev.space != deg.space) {
        for (int k = 0; k < vc_prev.size(); ++k) {
          V += u_r_prev[k + shift] * phi_V_j * T_Matrix[j][k];
          P += u_r_prev[k + shift] * phi_P_j * T_Matrix[j][k];
        }
      } else {
        if (c.min() != 0.0) {
          V += u_r_prev[j + shift] * phi_V_j;
          P += u_r_prev[j + shift] * phi_P_j;
        } else {
          V += c_nodal_value[j] * phi_V_j;
          P += c_nodal_value[j] * phi_P_j;
        }
      }
    }
    for (int j = elem.j_dimension() - r.n(); j < elem.j_dimension(); ++j) {
      VectorField
          phi_V_j = -1.0 * rho * elem.DtVelocity(q, j) + elem.GradPressure(q, j);
      double phi_P_j =
          -1.0 * kappaInv * elem.DtPressure(q, j) + elem.DivVelocity(q, j);
      V += u_r[j - elem.j_dimension() + r.n()] * phi_V_j;
      P += u_r[j - elem.j_dimension() + r.n()] * phi_P_j;
    }
    rho_1_V += w * V * V;
    rho_1_P += w * P * P;

    VectorField dual_V = 0.0;
    double dual_P = 0.0;
    for (int i = 0; i < elem.i_dimension(); ++i) {
      dual_V += u_star_r[i] * elem.VelocityTestSpace(q, i);
      dual_P += u_star_r[i] * elem.PressureTestSpace(q, i);
    }
    omega_1_V += w * dual_V * dual_V;
    omega_1_P += w * dual_P * dual_P;
  }

  //est_cell_time += Date() - Start_cell;
  //Date Start_face;

  for (int f = 0; f < c.Faces() - 2; ++f) {
    cell cf = c;
    bool bnd = false;
    if (u.OnBoundary(*c, f)) bnd = false;
    else cf = u.find_neighbour_cell(c, f);

    if (bnd) continue;

    row rf = u.find_row(cf());
    cell cf_prev = cf;
    vector<Scalar> cf_nodal_value;

    face ff = u.find_face(cf.Face(cf.Faces() - 2));
    cf_prev = u.find_cell(ff.Left());
    if (cf_prev == u.cells_end())
      cf_prev = u.find_overlap_cell(ff.Left());
    row rf_prev = u.find_row(cf_prev());

    const Scalar *u_rf = u(rf);
    const Scalar *u_rf_prev = u(rf_prev);

    DegreePair deg_cf_prev = u.GetDoF().GetDegree(*cf_prev);
    DegreePair deg_cf = u.GetDoF().GetDegree(*cf);
    int f1 = 0;
    Point f_c = u.find_face(c.Face(f))();
    Point f_cf = Origin;
    for (; f1 < cf.Faces() - 2; ++f1) {
      f_cf = u.find_face(cf.Face(f1))();
      Point pfcf = f_cf.CopyWithT(0.0);
      Point pfc = f_c.CopyWithT(0.0);
      if (pfc == pfcf) break;
    }

    SpaceTimeViscoAcousticFaceElement felem_c(*disc, u, c, f, prob->nL());
    SpaceTimeViscoAcousticFaceElement felem_cf(*disc, u, cf, f1, prob->nL());

    if (c.min() != 0) {
      cf_prev = u.find_cell_or_overlap_cell(ff.Left());
    } else {
      vector<Point> cf_nodal = u.GetDoF().GetNodalPoints(*cf);
      cf_nodal_value.resize(u.GetDoF().get_m(cf()));

      int cnt = 0;
      for (int pi = 0; pi < prob->Dim(); ++pi) {
        COMPONENT comp = prob->GetComponents()[pi];
        for (int i = 0; i < u.GetDoF().NodalPointsLocal(deg_cf.space); ++i) {
          cf_nodal_value[cnt] = prob->ut(cf_nodal[cnt].t(), cf_nodal[cnt], *c, comp);
          cnt++;
        }
      }
    }

    vector<Scalar> vcf(u.GetDoF().get_m(cf()));
    vector<Scalar> vcf_prev(u.GetDoF().get_m(cf_prev()));
    vector<vector<Scalar> > T_Matrix_f_prev(vcf.size());

    if (deg_cf.space != deg_cf_prev.space) {
      for (int i = 0; i < vcf.size(); ++i)
        T_Matrix_f_prev[i].resize(vcf_prev.size());
      for (int i = 0; i < vcf_prev.size(); ++i) {
        for (int k = 0; k < vcf_prev.size(); ++k)
          vcf_prev[k] = 0.0;
        vcf_prev[i] = 1.;
        PX_P0(vcf_prev, vcf, deg_cf_prev.space, deg_cf.space);
        for (int k = 0; k < vcf.size(); ++k)
          T_Matrix_f_prev[k][i] = vcf[k];
      }
    }

    double eta_v = 0.;
    int shift_1 = (deg_cf_prev.time - 1) * (rf_prev.n() / deg_cf_prev.time);

    for (int q = 0; q < felem_c.nQ(); ++q) {
      int q1 = find_q_id(felem_c.QPoint(q), felem_cf);
      VectorField V_Flux = 0.0;
      double P_Flux = 0.0;

      double w = felem_c.QWeight(q);
      VectorField Nq = felem_c.QNormal(q);

      double z_K = sqrt(prob->Rho(*c, c()) * prob->Kappa(*c, c()));
      double z_Kf = sqrt(prob->Rho(*cf, cf()) * prob->Kappa(*cf, cf()));
      double alpha1 = -1 / (z_K + z_Kf);
      double alpha2 = -z_Kf * z_K / (z_K + z_Kf);
      double alpha3 = z_Kf / (z_K + z_Kf);
      double alpha4 = z_K / (z_K + z_Kf);

      for (int j = 0; j < felem_c.j_dimension(); ++j) {
        VectorField V_j = 0.0;
        double P_j = 0.0;
        if (j < felem_c.j_dimension() - r.n()) {
          if (deg.space != deg_prev.space) {
            for (int k = 0; k < vc_prev.size(); ++k) {
              V_j += T_Matrix[j][k] * u_r_prev[k + shift] * felem_c.Velocity(q, j);
              P_j += T_Matrix[j][k] * u_r_prev[k + shift] * felem_c.Pressure(q, j);
            }
          } else {
            if (cf.min() != 0) {
              V_j += u_r_prev[j + shift] * felem_c.Velocity(q, j);
              P_j += u_r_prev[j + shift] * felem_c.Pressure(q, j);
            } else {
              V_j += c_nodal_value[j] * felem_c.Velocity(q, j);
              P_j += c_nodal_value[j] * felem_c.Pressure(q, j);
            }
          }
        } else {
          V_j += u_r[j - felem_c.j_dimension() + r.n()] * felem_c.Velocity(q, j);
          P_j += u_r[j - felem_c.j_dimension() + r.n()] * felem_c.Pressure(q, j);
        }
        V_Flux -= double(alpha1 * VFlux_c(V_j, Nq) + alpha3 * P_j) * Nq;
        P_Flux -= alpha2 * P_j + alpha4 * VFlux_c(V_j, Nq);
      }
      for (int j = 0; j < felem_cf.j_dimension(); ++j) {
        VectorField Vf_j = 0.0;
        double Pf_j = 0.0;

        if (j < felem_cf.j_dimension() - rf.n()) {
          if (deg_cf.space != deg_cf_prev.space) {
            for (int k = 0; k < vcf_prev.size(); ++k) {
              Vf_j += u_rf_prev[k + shift_1] * felem_cf.Velocity(q1, j) * T_Matrix_f_prev[j][k];
              Pf_j += u_rf_prev[k + shift_1] * felem_cf.Pressure(q1, j) * T_Matrix_f_prev[j][k];
            }
          } else {
            if (cf.min() != 0) {
              Vf_j += u_rf_prev[j + shift_1] * felem_cf.Velocity(q1, j);
              Pf_j += u_rf_prev[j + shift_1] * felem_cf.Pressure(q1, j);
            } else {
              Vf_j += cf_nodal_value[j] * felem_cf.Velocity(q1, j);
              Pf_j += cf_nodal_value[j] * felem_cf.Pressure(q1, j);
            }
          }
        } else {
          Vf_j += u_rf[j - felem_cf.j_dimension() + rf.n()] * felem_cf.Velocity(q1, j);
          Pf_j += u_rf[j - felem_cf.j_dimension() + rf.n()] * felem_cf.Pressure(q1, j);
        }
        V_Flux += double(alpha1 * VFlux_c(Vf_j, Nq) + alpha3 * Pf_j) * Nq;
        P_Flux += alpha2 * Pf_j + alpha4 * VFlux_c(Vf_j, Nq);
      }
      rho_2_V += w * V_Flux * V_Flux;
      rho_2_P += w * P_Flux * P_Flux;

      VectorField dual_V_Flux = 0.0;
      double dual_P_Flux = 0.0;
      for (int i = 0; i < felem_c.i_dimension(); ++i) {
        dual_V_Flux += u_star_r[i] * felem_c.VelocityTestSpace(q, i);
        dual_P_Flux += u_star_r[i] * felem_c.PressureTestSpace(q, i);
      }
      omega_2_V += w * dual_V_Flux * dual_V_Flux;
      omega_2_P += w * dual_P_Flux * dual_P_Flux;
    }
  }

  eta = sqrt(rho_1_V + rho_1_P) * sqrt(omega_1_V + omega_1_P)
        + sqrt(rho_2_V + rho_2_P) * sqrt(omega_2_V + omega_2_P);
  return eta;
}


vector<vector<Scalar>> STPGViscoAcousticAssemble::calc_TMatrix(int mc,
                                                               int prev_m,
                                                               int deg,
                                                               int prev_deg) const {
  vector<vector<Scalar>> T_Matrix(mc, vector<Scalar>(prev_m, 0.0));
  if (deg != prev_deg) {
    vector<Scalar> vc(mc, 0.0);
    vector<Scalar> vc_prev(prev_m, 0.0);

    for (int i = 0; i < vc_prev.size(); ++i) {
      ft::make_unit_vec(vc_prev, i);
      PX_P0(vc_prev, vc, prev_deg, deg);
      for (int k = 0; k < vc.size(); ++k)
        T_Matrix[k][i] = vc[k];
    }
  }
  return T_Matrix;
}


double STPGViscoAcousticAssemble::DualErrorEstimateJumpCell(const cell &c,
                                                            const Vector &u,
                                                            const Vector &u_star) const {
  double eta = 0.0;

  SpaceTimeViscoAcousticElement elem(*disc, u, c, prob->nL());
  row r = u.find_row(c());

  cell c_prev = c;
  vector<Scalar> c_nodal_value;

  face f = u.find_face(c.Face(c.Faces() - 2));
  c_prev = u.find_cell_or_overlap_cell(f.Left());
  row r_prev = u.find_row(c_prev());

  const Scalar *u_r_prev = u(r_prev);
  const Scalar *u_r = u(r);
  const Scalar *u_star_r = u_star(r);

  double h = dist(c.SpaceCell()[0], c.SpaceCell()[c.SpaceCell().Corners() - 1]);
  for (int cs = 0; cs < c.SpaceCell().Corners() - 1; ++cs)
    h = max(h, dist(c.SpaceCell()[cs], c.SpaceCell()[cs + 1]));
  double I = c.max() - c.min();

  double rho_1_V = 0.0;
  double rho_1_P = 0.0;
  double omega_1_V = 0.0;
  double omega_1_P = 0.0;

  DegreePair deg = u.GetDoF().GetDegree(*c);

  if (c.min() != 0) {
    c_prev = u.find_cell_or_overlap_cell(f.Left());
  } else {
    vector<Point> c_nodal = u.GetDoF().GetNodalPoints(*c);
    c_nodal_value.resize(c_nodal.size());
    int cnt = 0;
    for (int ti = 0; ti < deg.time; ++ti) {
      for (int pi = 0; pi < prob->Dim(); ++pi) {
        COMPONENT comp = prob->GetComponents()[pi];
        for (int i = 0; i < u.GetDoF().NodalPointsLocal(deg.space); ++i) {
          c_nodal_value[cnt] = prob->ut(c_nodal[cnt].CopyWithT(0.0).t(), c_nodal[cnt].CopyWithT(0.0), *c, comp);
          cnt++;
        }
      }
    }
  }

  DegreePair deg_prev = u.GetDoF().GetDegree(*c_prev);
  int shift = (deg_prev.time - 1) * (r_prev.n() / deg_prev.time);

  vector<vector<Scalar>> T_Matrix{calc_TMatrix(u.GetDoF().get_m(c()),
                                               u.GetDoF().get_m(c_prev()),
                                               deg.space,
                                               deg_prev.space)};


  double rho = prob->Rho(*c, c());
  double kappaInv = 1.0 / prob->Kappa(*c, c());

  for (int q = 0; q < elem.nQ(); ++q) {
    VectorField V = 0.0;
    double P = 0.0;

    double w = elem.QWeight(q);
    Point QP = elem.QPoint(q);
    double t = QP.t();
    for (int j = 0; j < elem.j_dimension() - r.n(); ++j) {
      VectorField phi_V_j = -1.0 * rho * elem.DtVelocity(q, j) + elem.GradPressure(q, j);
      double phi_P_j = -1.0 * kappaInv * elem.DtPressure(q, j) + elem.DivVelocity(q, j);

      if (deg_prev.space != deg.space) {
        for (int k = 0; k < T_Matrix[j].size(); ++k) {
          V += u_r_prev[k + shift] * phi_V_j * T_Matrix[j][k];
          P += u_r_prev[k + shift] * phi_P_j * T_Matrix[j][k];
        }
      } else {
        if (c.min() != 0.0) {
          V += u_r_prev[j + shift] * phi_V_j;
          P += u_r_prev[j + shift] * phi_P_j;
        } else {
          V += c_nodal_value[j] * phi_V_j;
          P += c_nodal_value[j] * phi_P_j;
        }
      }
    }
    for (int j = elem.j_dimension() - r.n(); j < elem.j_dimension(); ++j) {
      VectorField
          phi_V_j = -1.0 * rho * elem.DtVelocity(q, j) + elem.GradPressure(q, j);
      double phi_P_j =
          -1.0 * kappaInv * elem.DtPressure(q, j) + elem.DivVelocity(q, j);
      V += u_r[j - elem.j_dimension() + r.n()] * phi_V_j;
      P += u_r[j - elem.j_dimension() + r.n()] * phi_P_j;
    }

    VectorField V_rhs = prob->F_Velocity(t, *c, QP);
    Scalar P_rhs = prob->F_Pressure(t, *c, QP);

    rho_1_V += w * (V * V + V_rhs * V_rhs);
    rho_1_P += w * (P * P + P_rhs * P_rhs);
  }

  vector<VectorField> dual_V_projection(deg.time, zero);
  vector<double> dual_P_projection(deg.time, 0.0);

  double K_vol = 0;
  for (int q = 0; q < elem.nQ(); ++q) {
    K_vol += elem.QWeight(q);
    for (int i = 0; i < elem.i_dimension(); ++i) {
      for (int l = 0; l < deg.time; ++l) {
        dual_V_projection[l] += elem.QWeight(q) * u_star_r[i % (elem.i_dimension() / deg.time)
                                                           + l * elem.i_dimension() / deg.time] *
                                elem.VelocityTestSpace(q, i);
        dual_P_projection[l] += elem.QWeight(q) * u_star_r[i % (elem.i_dimension() / deg.time)
                                                           + l * elem.i_dimension() / deg.time] *
                                elem.PressureTestSpace(q, i);
      }
    }
  }
  for (int l = 0; l < deg.time; ++l) {
    dual_V_projection[l] *= 1. / K_vol;
    dual_P_projection[l] *= 1. / K_vol;
  }

  const Shape &time_shape_testspace(disc->GetTimeShape(deg.time - 1));

  vector<double> V_Fluxes(c.Faces() - 2);
  vector<double> P_Fluxes(c.Faces() - 2);
  vector<double> V_dual_Fluxes(c.Faces() - 2);
  vector<double> P_dual_Fluxes(c.Faces() - 2);

  for (int f = 0; f < c.Faces() - 2; ++f) {
    V_Fluxes[f] = 0;
    V_dual_Fluxes[f] = 0;
    P_Fluxes[f] = 0;
    P_dual_Fluxes[f] = 0;

    cell cf = c;
    bool bnd = false;
    if (u.OnBoundary(*c, f)) bnd = true;
    else cf = u.find_neighbour_cell(c, f);

    if (bnd) continue;

    double hf = dist(c.SpaceCell()[f], c.SpaceCell()[(f + 1) % c.SpaceCell().Corners()]);
    hf = sqrt(hf * hf + (c.max() - c.min()) * (c.max() - c.min()));

    row rf = u.find_row(cf());
    cell cf_prev = cf;
    vector<Scalar> cf_nodal_value;

    face ff = u.find_face(cf.Face(cf.Faces() - 2));
    cf_prev = u.find_cell_or_overlap_cell(ff.Left());
    row rf_prev = u.find_row(cf_prev());

    const Scalar *u_rf = u(rf);
    const Scalar *u_rf_prev = u(rf_prev);
    const Scalar *u_star_rf = u_star(rf);

    DegreePair deg_cf_prev = u.GetDoF().GetDegree(*cf_prev);
    DegreePair deg_cf = u.GetDoF().GetDegree(*cf);

    int f1 = 0;
    Point f_c = u.find_face(c.Face(f))();
    Point f_cf = Origin;
    for (; f1 < cf.Faces() - 2; ++f1) {
      f_cf = u.find_face(cf.Face(f1))();
      Point pfcf = f_cf.WithT(0.0);
      Point pfc = f_c.WithT(0.0);
      if (pfc == pfcf) break;
    }

    SpaceTimeViscoAcousticFaceElement felem_c(*disc, u, c, f, prob->nL());
    SpaceTimeViscoAcousticFaceElement felem_cf(*disc, u, cf, f1, prob->nL());
    SpaceTimeViscoAcousticElement elem_cf(*disc, u, cf, prob->nL());

    if (c.min() != 0) {
      //face ff = u.find_face(cf.Face(cf.Faces() - 2));
      cf_prev = u.find_cell_or_overlap_cell(ff.Left());
    } else {
      vector<Point> cf_nodal = u.GetDoF().GetNodalPoints(*cf);
      cf_nodal_value.resize(u.GetDoF().get_m(cf()));

      int cnt = 0;
      for (int pi = 0; pi < prob->Dim(); ++pi) {
        COMPONENT comp = prob->GetComponents()[pi];
        for (int i = 0; i < u.GetDoF().NodalPointsLocal(elem_cf.deg.space); ++i) {
          cf_nodal_value[cnt] = prob->ut(cf_nodal[cnt].WithT(0.0).t(), cf_nodal[cnt].WithT(0.0), *c, comp);
          cnt++;
        }
      }
    }

    vector<vector<Scalar>> T_Matrix_f_prev{calc_TMatrix(u.GetDoF().get_m(cf()),
                                                        u.GetDoF().get_m(cf_prev()),
                                                        elem_cf.deg.space,
                                                        deg_cf_prev.space)};

    double eta_v = 0.;
    int shift_1 = (deg_cf_prev.time - 1) * (rf_prev.n() / deg_cf_prev.time);

    double fxI = 0;

    for (int q = 0; q < felem_c.nQ(); ++q) {
      int q1 = find_q_id(felem_c.QPoint(q), felem_cf);
      VectorField V_Flux = 0.0;
      double P_Flux = 0.0;

      double w = felem_c.QWeight(q);
      VectorField Nq = felem_c.QNormal(q);

      fxI += w;

      for (int j = 0; j < felem_c.j_dimension(); ++j) {
        VectorField V_j = 0.0;
        double P_j = 0.0;
        if (j < felem_c.j_dimension() - r.n()) {
          if (deg.space != deg_prev.space) {
            for (int k = 0; k < T_Matrix[j].size(); ++k) {
              V_j += T_Matrix[j][k] * u_r_prev[k + shift] * felem_c.Velocity(q, j);
              P_j += T_Matrix[j][k] * u_r_prev[k + shift] * felem_c.Pressure(q, j);
            }
          } else {
            if (cf.min() != 0) {
              V_j += u_r_prev[j + shift] * felem_c.Velocity(q, j);
              P_j += u_r_prev[j + shift] * felem_c.Pressure(q, j);
            } else {
              V_j += c_nodal_value[j] * felem_c.Velocity(q, j);
              P_j += c_nodal_value[j] * felem_c.Pressure(q, j);
            }
          }
        } else {
          V_j += u_r[j - felem_c.j_dimension() + r.n()] * felem_c.Velocity(q, j);
          P_j += u_r[j - felem_c.j_dimension() + r.n()] * felem_c.Pressure(q, j);
        }
        V_Flux -= V_j;
        P_Flux -= P_j;
      }
      for (int j = 0; j < felem_cf.j_dimension(); ++j) {
        VectorField Vf_j = 0.0;
        double Pf_j = 0.0;

        if (j < felem_cf.j_dimension() - rf.n()) {
          if (deg_cf.space != deg_cf_prev.space) {
            for (int k = 0; k < T_Matrix_f_prev[j].size(); ++k) {
              double value = u_rf_prev[k + shift_1] * T_Matrix_f_prev[j][k];
              Vf_j += value * felem_cf.Velocity(q1, j);
              Pf_j += value * felem_cf.Pressure(q1, j);
            }
          } else {
            if (cf.min() != 0) {
              Vf_j += u_rf_prev[j + shift_1] * felem_cf.Velocity(q1, j);
              Pf_j += u_rf_prev[j + shift_1] * felem_cf.Pressure(q1, j);
            } else {
              Vf_j += cf_nodal_value[j] * felem_cf.Velocity(q1, j);
              Pf_j += cf_nodal_value[j] * felem_cf.Pressure(q1, j);
            }
          }
        } else {
          Vf_j += u_rf[j - felem_cf.j_dimension() + rf.n()] * felem_cf.Velocity(q1, j);
          Pf_j += u_rf[j - felem_cf.j_dimension() + rf.n()] * felem_cf.Pressure(q1, j);
        }
        V_Flux += Vf_j;
        P_Flux += Pf_j;
      }
      V_Fluxes[f] += w * V_Flux * V_Flux;
      P_Fluxes[f] += w * P_Flux * P_Flux;
    }

    //int time_deg_cf = u.GetDoF().get_time_deg(cf);

    vector<VectorField> dual_V_projection_cf(deg_cf.time, zero);
    vector<double> dual_P_projection_cf(deg_cf.time, 0.0);


    double Kf_vol = 0;

    for (int q = 0; q < elem_cf.nQ(); ++q) {
      double w = elem_cf.QWeight(q);
      Kf_vol += w;
      for (int i = 0; i < elem_cf.i_dimension(); ++i) {
        for (int l = 0; l < deg_cf.time; ++l) {
          int u_star_idx = elem_cf.i_dimension() / deg_cf.time;
          double u_star_rf_elem = u_star_rf[i % u_star_idx + l * u_star_idx];
          double value = w * u_star_rf_elem;
          dual_V_projection_cf[l] += value * elem_cf.VelocityTestSpace(q, i);
          dual_P_projection_cf[l] += value * elem_cf.PressureTestSpace(q, i);
        }
      }
    }

    for (int l = 0; l < deg_cf.time; ++l) {
      dual_V_projection_cf[l] *= 1. / Kf_vol;
      dual_P_projection_cf[l] *= 1. / Kf_vol;
    }

    const Shape &time_shape_testspace_cf(disc->GetTimeShape(deg_cf.time - 1));

    auto &timeQuad = disc->GetTimeQuad(deg_cf.time - 1);

    for (int q = 0; q < timeQuad.size(); ++q) {
      Point QP = timeQuad.QPoint(q);
      double w = timeQuad.Weight(q);

      Scalar values_dual_V = 0.0;
      Scalar values_dual_P = 0.0;
      for (int l = 0; l < time_shape_testspace.size(); ++l) {
        double timeValue = time_shape_testspace(QP, l);
        values_dual_V += dual_V_projection[l] * felem_c.QNormal(0) * timeValue;
        values_dual_P += dual_P_projection[l] * timeValue;
      }

      Scalar values_dual_V_cf = 0.0;
      Scalar values_dual_P_cf = 0.0;
      for (int l = 0; l < time_shape_testspace_cf.size(); ++l) {
        double timeValue = time_shape_testspace_cf(QP, l);
        values_dual_V_cf += dual_V_projection_cf[l] * felem_c.QNormal(0) * timeValue;
        values_dual_P_cf += dual_P_projection_cf[l] * timeValue;
      }
      V_dual_Fluxes[f] += w * fxI * pow(values_dual_V_cf - values_dual_V, 2);
      P_dual_Fluxes[f] += w * fxI * pow(values_dual_P_cf - values_dual_P, 2);
    }
  }

  double eta_f = 0;
  for (int f = 0; f < c.Faces() - 2; ++f) {
    omega_1_V += V_dual_Fluxes[f];
    omega_1_P += P_dual_Fluxes[f];
    eta_f += sqrt(V_Fluxes[f]) * sqrt(V_dual_Fluxes[f])
             + sqrt(P_Fluxes[f]) * sqrt(P_dual_Fluxes[f]);
  }
  eta = sqrt(rho_1_P) * sqrt(h) * sqrt(omega_1_P)
        + sqrt(rho_1_V) * sqrt(h) * sqrt(omega_1_V) + 0.5 * eta_f;
  return eta;
}

double STPGViscoAcousticAssemble::ErrorEstimateCell(const cell &c,
                                                    const Vector &u) const {
  Exit("Not implemented");
}

double STPGViscoAcousticAssemble::Goal_Functional(const Vector &u,
                                                  Point a,
                                                  Point b) const {
  Scalar e = 0.0;
  //mout << u.norm() << " " << prob->nL() << endl;
  if (a.t() != b.t()) { // volume evaluation
    Exit("ToDo");
  } else { // evaluation at some time point
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != a.t()) continue;
      SpaceTimeViscoAcousticElement elem(*disc, u, c, prob->nL());
      row r = u.find_row(c());
      const int dim = c.dim();
      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        Point QP = elem.QPoint(q).CopyWithT(0.0);
        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;
        const double w = elem.SpaceQWeight(q);

        VectorField j_V;
        for (int i = 0; i < dim; ++i) {
          COMPONENT comp = prob->GetComponents()[i];
          j_V[i] = prob->ut_dual(QP, comp, a, b);
        }

        vector<Scalar> j_P(2 + prob->nL());
        for (int k = 1; k <= 1 + prob->nL(); ++k) {
          COMPONENT comp = prob->GetComponents()[dim - 1 + k];
          j_P[k] = prob->ut_dual(QP, COMPONENT::P0, a, b);

          //mout << j_P[k] << endl;
        }
        const int start_idx = elem.deg.time * (elem.j_dimension() / (elem.deg.time + 1));

        Scalar scale = 0.0;
        for (int j = start_idx; j < elem.j_dimension(); ++j) {
          scale += elem.Velocity(q, j)[0];
        }
        VectorField V = 0.0;
        vector<Scalar> P(2 + prob->nL());
        for (int j = start_idx; j < elem.j_dimension(); ++j) {
          VectorField V_j = elem.Velocity(q, j) / scale;
          Scalar P_j = elem.Pressure(q, j) / scale;
          V += u(r, j - elem.j_dimension() + r.n()) * V_j;
          P[elem.variable(j)] += u(r, j - elem.j_dimension() + r.n()) * P_j;
        }
        e += w * (j_V * V);
        for (int k = 1; k <= 1 + prob->nL(); ++k) {
          e += w * j_P[k] * P[k];
        }
      }
    }
  }
  return abs(PPM->SumOnCommSplit(e, u.CommSplit()));
}

double STPGViscoAcousticAssemble::Energy(const Vector &u, Point a, Point b) const {
  Scalar e = 0;
  if (a.t() != b.t()) { // volume evaluation
    Exit("TODO")
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      SpaceTimeViscoAcousticElement elem(*disc, u, c, prob->nL());
      row r = u.find_row(c());

      cell c_prev = c;
      vector<Point> c_nodal;

      if (c.min() != 0) {
        face f = u.find_face(c.Face(c.Faces() - 2));
        c_prev = u.find_cell(f.Left());
        if (c_prev == u.cells_end())
          c_prev = u.find_overlap_cell(f.Left());
      } else {
        c_nodal = u.GetDoF().GetNodalPoints(*c);
        for (int i = 0; i < c_nodal.size(); ++i) c_nodal[i] = c_nodal[i].CopyWithT(0.0);
      }

      row r_prev = u.find_row(c_prev());

      DegreePair deg_c_prev = u.GetDoF().GetDegree(*c_prev);

      vector<Scalar> vc(u.GetDoF().get_m(c()));
      vector<Scalar> vc_prev(u.GetDoF().get_m(c_prev()));

      vector<vector<Scalar> > T_Matrix(vc.size());
      for (int i = 0; i < vc.size(); ++i)
        T_Matrix[i].resize(vc_prev.size());

      for (int i = 0; i < vc_prev.size(); ++i) {
        for (int k = 0; k < vc_prev.size(); ++k)
          vc_prev[k] = 0.0;
        vc_prev[i] = 1.;
        PX_P0(vc_prev, vc, deg_c_prev.space, elem.deg.space);

        for (int k = 0; k < vc.size(); ++k)
          T_Matrix[k][i] = vc[k];
      }

      for (int q = 0; q < elem.nQ(); ++q) {
        double w = elem.QWeight(q);
        Point QP = elem.QPoint(q);
        Scalar P = 0.0;
        VectorField V = 0.0;
        if ((a[0] <= QP[0]) && (QP[0] <= b[0])
            && (a[1] <= QP[1]) && (QP[1] <= b[1])
            && a.t() <= QP.t() && QP.t() <= b.t()) {
          for (int j = 0; j < elem.j_dimension(); ++j) {
            Scalar P_j = elem.Pressure(q, j);
            VectorField V_j = elem.Velocity(q, j);

            if (j < elem.j_dimension() - r.n()) {
              if (c.min() == 0) {
                P += prob->ut_Pressure(0, c_nodal[j], *c) * P_j;
                V += prob->ut_Velocity(0, c_nodal[j], *c) * V_j;
              } else
                for (int k = 0; k < vc_prev.size(); ++k) {
                  int index = j + (elem.deg.time - 1) * (r_prev.n() / elem.deg.time);
                  V += u(r_prev, index) * V_j * T_Matrix[j][k];
                  P += u(r_prev, index) * P_j * T_Matrix[j][k];
                }
            } else {
              V += u(r, j - elem.j_dimension() + r.n()) * V_j;
              P += u(r, j - elem.j_dimension() + r.n()) * P_j;
            }
          }
        }
        e += 0.5 * w * (V * V + P * P);
      }
    }
  } else { // evaluation at some time point
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != a.t()) continue;
      SpaceTimeViscoAcousticElement elem(*disc, u, c, prob->nL());
      row r = u.find_row(c());
      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        Point QP = elem.QPoint(q).CopyWithT(0.0);
        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;
        const double w = elem.SpaceQWeight(q);
        Scalar scale = 0.0;
        const int start_idx = elem.deg.time * (elem.j_dimension() / (elem.deg.time + 1));
        for (int j = start_idx; j < elem.j_dimension(); ++j)
          scale += elem.Velocity(q, j)[0];

        VectorField V = 0.0;
        vector<Scalar> P(2 + prob->nL());
        for (int j = start_idx; j < elem.j_dimension(); ++j) {
          VectorField V_j = elem.Velocity(q, j) / scale;
          Scalar P_j = elem.Pressure(q, j) / scale;

          V += u(r, j - elem.j_dimension() + r.n()) * V_j;
          P[elem.variable(j)] += u(r, j - elem.j_dimension() + r.n()) * P_j;
        }
        e += 0.5 * w * V * V;
        for (int k = 1; k <= 1 + prob->nL(); ++k)
          e += 0.5 * w * P[k] * P[k];
      }
    }
  }
  return abs(PPM->SumOnCommSplit(e, u.CommSplit()));
}

