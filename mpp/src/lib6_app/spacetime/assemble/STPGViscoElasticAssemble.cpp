#include "STPGViscoElasticAssemble.hpp"

using namespace std;


void STPGViscoElasticAssemble::System(Matrix &M, Vector &RHS) const {
  Date Start;

  vout(2) << "   Problem Size: " << M.pSize() << " start assemble " << Start;

  Time t_cell;
  Time t_face;
  M = 0;
  Vector u_0(0.0, RHS);
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    cell c_prev = M.find_previous_cell(c);
    if (c_prev != c) continue;
    Date Start_cell;
    SpaceTimeViscoElasticElement elem(*disc, M, c, prob->nL());
    row r_c = M.find_row(c());
    int c_deg = RHS.GetDoF().get_space_deg(*c);
    DGRowEntries M_c(M, *c, *c, false);
    const double rho = prob->Rho(c());
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      for (int i = 0; i < elem.i_dimension(); ++i) {
        int var_i = elem.variable(i);
        VectorField V_i = elem.VelocityTestSpace(q, i);
        Tensor S_i = elem.StressTestSpace(q, i);
        for (int j = 0; j < elem.j_dimension() - r_c.n(); ++j) {
          int var_j = elem.variable(j);
          double mu = prob->Mu(c(), var_j);
          //double kappa = problem->Kappa( c(),var_j );
          double tau = prob->Tau(c(), var_j);
          double s = 1 / (2.0 * mu);
          //double ss = 3/kappa - 1/(6.0 * mu);
          double lambda = prob->Lambda(c(), var_j);
          double ss = -lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));

          VectorField DtV_j = elem.DtVelocity(q, j);
          Tensor DtS_j = elem.DtStress(q, j);
          Tensor S_j = elem.Stress(q, j);
          VectorField DivS_j = elem.DivStress(q, j);
          Tensor E_j = elem.Strain(q, j);

          if (var_i != var_j) DtS_j = Zero;

          M_c(i, j) += w * ((rho * (DtV_j * V_i) + s * Frobenius(DtS_j, S_i)
                             + ss * trace(DtS_j) * trace(S_i))
                            - ((DivS_j * V_i) + Frobenius(E_j, S_i)));
          if (var_j > 1 && var_i == var_j)
            M_c(i, j) += w / tau * (
                s * Frobenius(S_j, S_i) + ss * trace(S_j) * trace(S_i));
        }
      }
    }

    int time_deg = RHS.GetDoF().get_time_deg(*c);
    vector<Point> c_nodal = RHS.GetDoF().GetNodalPoints(*c);
    int probDim = prob->Dim();
    int cnt = 0;
    for (int ti = 0; ti < time_deg; ++ti) {
      for (int pi = 0; pi < probDim; ++pi) {
        COMPONENT comp = prob->GetComponents()[pi];
        for (int i = 0; i < RHS.GetDoF().NodalPointsLocal(c_deg, pi); ++i) {
          u_0(r_c, cnt) = prob->ut(c_nodal[cnt].CopyWithT(0.0), comp);
          cnt++;
        }
      }
    }

    t_cell += Date() - Start_cell;
    Date Start_face;

    const double c_s = prob->C_s(c());
    const double c_p = prob->C_p(c());
    for (int f = 0; f < c.Faces() - 2; ++f) {
      cell cf = c;
      bool bnd = false;
      //int bnd_id = M.BoundaryFaces::Part(c.Face(f));
      int bnd_id = -1;
      if (M.OnBoundary(*c, f)) {
        bnd = true;
        bnd_id = prob->BndID(c.Face(f));
      } else cf = M.find_neighbour_cell(c, f);
      SpaceTimeViscoElasticFaceElement felem(*disc, M, c, f, prob->nL());
      int cf_deg = RHS.GetDoF().get_space_deg(*cf);
      int f1 = M.find_neighbour_face_id(c.Face(f), cf);
      SpaceTimeViscoElasticFaceElement felem_1(*disc, M, cf, f1, prob->nL());
      DGRowEntries M_cf(M, *c, *cf, false);
      row r_cf = M.find_row(cf());

      const double c_s_f = prob->C_s(cf());
      const double c_p_f = prob->C_p(cf());
      const double rho_f = prob->Rho(cf());

      const double alpha1 = 1.0 / (rho * c_p + rho_f * c_p_f);
      const double alpha2 = rho * c_p * alpha1;
      const double alpha3 = rho_f * c_p_f * alpha1;
      const double alpha4 = (rho * c_p * rho_f * c_p_f) * alpha1;

      const double alpha5 = 1.0 / (rho * c_s + rho_f * c_s_f);
      const double alpha6 = rho * c_s * alpha5;
      const double alpha7 = rho_f * c_s_f * alpha5;
      const double alpha8 = (rho * c_s * rho_f * c_s_f) * alpha5;

      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        double w = felem.QWeight(q);
        VectorField Nq = felem.QNormal(q);
        VectorField Tq = felem.QTangent(q);
        for (int i = 0; i < felem.i_dimension(); ++i) {
          if (felem.is_zero_test(q, i)) continue;
          VectorField V_i = felem.VelocityTestSpace(q, i);
          Tensor S_i = felem.StressTestSpace(q, i);
          const double SFlux_NN = SFlux_c(S_i, Nq, Nq);
          const double SFlux_TN = SFlux_c(S_i, Tq, Nq);
          const double VFlux_N = VFlux_c(V_i, Nq);
          const double VFlux_T = VFlux_c(V_i, Tq);
          for (int j = 0; j < felem.j_dimension() - r_c.n(); ++j) {
            if (felem.is_zero(q, j)) continue;
            Tensor S_j = felem.Stress(q, j);
            VectorField V_j = felem.Velocity(q, j);
            M_c(i, j)
                += w * (alpha1 * SFlux_c(S_j, Nq, Nq) * SFlux_NN
                        + alpha2 * SFlux_c(S_j, Nq, Nq) * VFlux_N
                        + alpha3 * VFlux_c(V_j, Nq) * SFlux_NN
                        + alpha4 * VFlux_c(V_j, Nq) * VFlux_N
                        + alpha5 * SFlux_c(S_j, Tq, Nq) * SFlux_TN
                        + alpha6 * SFlux_c(S_j, Tq, Nq) * VFlux_T
                        + alpha7 * VFlux_c(V_j, Tq) * SFlux_TN
                        + alpha8 * VFlux_c(V_j, Tq) * VFlux_T);
          }
          for (int j = 0; j < felem_1.j_dimension() - r_cf.n(); ++j) {
            if (felem_1.is_zero(q1, j)) continue;
            Tensor S1_j = felem_1.Stress(q1, j);
            VectorField V1_j = felem_1.Velocity(q1, j);
            M_cf(i, j)
                -=
                w * (alpha1 * SFlux_cf(S1_j, Nq, Nq, bnd, bnd_id) * SFlux_NN
                     + alpha2 * SFlux_cf(S1_j, Nq, Nq, bnd, bnd_id) * VFlux_N
                     + alpha3 * VFlux_cf(V1_j, Nq, bnd, bnd_id) * SFlux_NN
                     + alpha4 * VFlux_cf(V1_j, Nq, bnd, bnd_id) * VFlux_N
                     + alpha5 * SFlux_cf(S1_j, Tq, Nq, bnd, bnd_id) * SFlux_TN
                     + alpha6 * SFlux_cf(S1_j, Tq, Nq, bnd, bnd_id) * VFlux_T
                     + alpha7 * VFlux_cf(V1_j, Tq, bnd, bnd_id) * SFlux_TN
                     + alpha8 * VFlux_cf(V1_j, Tq, bnd, bnd_id) * VFlux_T);
          }
        }
      }
    }
    t_face += Date() - Start_face;
  }
  u_0.Accumulate();
  Vector Bu_0 = M * u_0;
  RHS -= Bu_0;

  t_cell.Max();
  t_face.Max();
  vout(3) << "      cell init assemble time " << t_cell << endl;
  vout(3) << "      face init assemble time " << t_face << endl;

  M = 0;
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    Date Start_cell;
    SpaceTimeViscoElasticElement elem(*disc, M, c, prob->nL());
    //mout << OUT(c()) << OUT(elem.j_dimension()) << OUT(elem.nQ()) << t_cell << endl;
    cell c_prev = M.find_previous_cell(c);
    row r_c = M.find_row(c());
    int c_deg = RHS.GetDoF().get_space_deg(*c);
    int c_prev_deg = RHS.GetDoF().get_space_deg(*c_prev);
    vector<vector<Scalar>> T_Matrix(RHS.GetDoF().get_m(c()));
    vector<Scalar> vc_diag(RHS.GetDoF().get_m(c()));
    vector<Scalar> vc(RHS.GetDoF().get_m(c_prev()));
    if (c_deg != c_prev_deg) {
      for (int i = 0; i < vc_diag.size(); ++i)
        T_Matrix[i].resize(vc.size());
      for (int i = 0; i < vc.size(); ++i) {
        for (int k = 0; k < vc.size(); ++k)
          vc[k] = 0.0;
        vc[i] = 1.;
        PX_P0(vc, vc_diag, c_prev_deg, c_deg, RHS);
        for (int k = 0; k < vc_diag.size(); ++k)
          T_Matrix[k][i] = vc_diag[k];
      }
    }
    DGRowEntries M_c(M, *c, *c, false);
    DGRowEntries M_c_prev(M, *c, *c_prev, true);
    const double rho = prob->Rho(c());
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);

      VectorField V_rhs(prob->F(QP, COMPONENT::V_X), prob->F(QP, COMPONENT::V_Y));
      Tensor S_rhs(prob->F(QP, COMPONENT::S0_x), prob->F(QP, COMPONENT::S0_xy),
                   prob->F(QP, COMPONENT::S0_xy), prob->F(QP, COMPONENT::S0_y));
      if (c.dim() != 2) Exit("TODO");

      for (int i = 0; i < elem.i_dimension(); ++i) {
        int var_i = elem.variable(i);
        VectorField V_i = elem.VelocityTestSpace(q, i);
        Tensor S_i = elem.StressTestSpace(q, i);

        if (prob->HasRHS())
          RHS(c(), i) += w * (V_rhs * V_i + Frobenius(S_rhs, S_i));

        for (int j = 0; j < elem.j_dimension() - r_c.n(); ++j) {
          if (c_prev == c) break;
          int var_j = elem.variable(j);
          double mu = prob->Mu(c(), var_j);
          //double kappa = problem->Kappa( c(),var_j );
          double tau = prob->Tau(c(), var_j);
          double s = 1 / (2.0 * mu);
          //double ss = 3/kappa - 1/(6.0 * mu);
          double lambda = prob->Lambda(c(), var_j);
          double ss = -lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));

          VectorField DtV_j = elem.DtVelocity(q, j);
          Tensor DtS_j = elem.DtStress(q, j);
          Tensor S_j = elem.Stress(q, j);
          VectorField DivS_j = elem.DivStress(q, j);
          Tensor E_j = elem.Strain(q, j);

          if (var_i != var_j) DtS_j = Zero;

          if (c_deg == c_prev_deg) {
            M_c_prev(i, j) +=
                w * ((rho * (DtV_j * V_i) + s * Frobenius(DtS_j, S_i)
                      + ss * trace(DtS_j) * trace(S_i))
                     - ((DivS_j * V_i) + Frobenius(E_j, S_i)));
            if (var_j > 1 && var_i == var_j)
              M_c_prev(i, j) += w / tau * (
                  s * Frobenius(S_j, S_i) + ss * trace(S_j) * trace(S_i));
          } else {
            for (int k = 0; k < vc.size(); ++k) {
              M_c_prev(i, k) +=
                  w * ((rho * (DtV_j * V_i) + s * Frobenius(DtS_j, S_i)
                        + ss * trace(DtS_j) * trace(S_i))
                       - ((DivS_j * V_i) + Frobenius(E_j, S_i)))
                  * T_Matrix[j][k];
              if (var_j > 1 && var_i == var_j)
                M_c_prev(i, k) += w / tau * (
                    s * Frobenius(S_j, S_i)
                    + ss * trace(S_j) * trace(S_i)) * T_Matrix[j][k];
            }
          }
        }
        for (int j = elem.j_dimension() - r_c.n(); j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          double mu = prob->Mu(c(), var_j);
          //double kappa = problem->Kappa( c(),var_j );
          double tau = prob->Tau(c(), var_j);
          double s = 1 / (2.0 * mu);
          //double ss = 3/kappa - 1/(6.0 * mu);
          double lambda = prob->Lambda(c(), var_j);
          double ss = -lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));

          VectorField DtV_j = elem.DtVelocity(q, j);
          Tensor DtS_j = elem.DtStress(q, j);
          Tensor S_j = elem.Stress(q, j);
          VectorField DivS_j = elem.DivStress(q, j);
          Tensor E_j = elem.Strain(q, j);

          if (var_i != var_j) DtS_j = Zero;

          M_c(i, j - elem.j_dimension() + r_c.n())
              += w * ((rho * (DtV_j * V_i) + s * Frobenius(DtS_j, S_i)
                       + ss * trace(DtS_j) * trace(S_i))
                      - ((DivS_j * V_i) + Frobenius(E_j, S_i)));
          if (var_j > 1 && var_i == var_j)
            M_c(i, j - elem.j_dimension() + r_c.n()) += w / tau * (
                s * Frobenius(S_j, S_i) + ss * trace(S_j) * trace(S_i));
        }
      }
    }

    t_cell += Date() - Start_cell;
    Date Start_face;

    const double c_s = prob->C_s(c());
    const double c_p = prob->C_p(c());
    for (int f = 0; f < c.Faces() - 2; ++f) {
      cell cf = c;
      bool bnd = false;
      //int bnd_id = M.BoundaryFaces::Part(c.Face(f));
      int bnd_id = -1;
      if (M.OnBoundary(*c, f)) {
        bnd = true;
        bnd_id = prob->BndID(c.Face(f));
      } else cf = M.find_neighbour_cell(c, f);
      SpaceTimeViscoElasticFaceElement *felem
          = new SpaceTimeViscoElasticFaceElement(*disc, M, c, f, prob->nL());
      int cf_deg = RHS.GetDoF().get_space_deg(*cf);
      cell cf_prev = M.find_previous_cell(cf);
      int f1 = M.find_neighbour_face_id(c.Face(f), cf);
      SpaceTimeViscoElasticFaceElement *felem_1
          = new SpaceTimeViscoElasticFaceElement(*disc, M, cf, f1, prob->nL());
      int cf_prev_deg = RHS.GetDoF().get_space_deg(*cf_prev);

      vector<vector<Scalar>> T_Matrix_f_prev(RHS.GetDoF().get_m(cf()));
      vector<Scalar> vcf(RHS.GetDoF().get_m(cf()));
      vector<Scalar> vcf_prev(RHS.GetDoF().get_m(cf_prev()));

      if (cf_deg != cf_prev_deg) {
        for (int i = 0; i < vcf.size(); ++i)
          T_Matrix_f_prev[i].resize(vcf_prev.size());
        for (int i = 0; i < vcf_prev.size(); ++i) {
          for (int k = 0; k < vcf_prev.size(); ++k)
            vcf_prev[k] = 0.0;
          vcf_prev[i] = 1.;
          PX_P0(vcf_prev, vcf, cf_prev_deg, cf_deg, RHS);
          for (int k = 0; k < vcf.size(); ++k)
            T_Matrix_f_prev[k][i] = vcf[k];
        }
      }

      DGRowEntries M_cf(M, *c, *cf, false);
      DGRowEntries M_cf_prev(M, *c, *cf_prev, true);
      row r_cf = M.find_row(cf());

      const double c_s_f = prob->C_s(cf());
      const double c_p_f = prob->C_p(cf());
      const double rho_f = prob->Rho(cf());

      const double alpha1 = 1.0 / (rho * c_p + rho_f * c_p_f);
      const double alpha2 = rho * c_p * alpha1;
      const double alpha3 = rho_f * c_p_f * alpha1;
      const double alpha4 = (rho * c_p * rho_f * c_p_f) * alpha1;

      const double alpha5 = 1.0 / (rho * c_s + rho_f * c_s_f);
      const double alpha6 = rho * c_s * alpha5;
      const double alpha7 = rho_f * c_s_f * alpha5;
      const double alpha8 = (rho * c_s * rho_f * c_s_f) * alpha5;

      for (int q = 0; q < felem->nQ(); ++q) {
        int q1 = find_q_id(felem->QPoint(q), *felem_1);
        double w = felem->QWeight(q);
        VectorField Nq = felem->QNormal(q);
        VectorField Tq = felem->QTangent(q);
        Point QP = felem->QPoint(q);

        VectorField V_bnd = prob->ut_Velocity(QP);
        Tensor S_bnd = prob->ut_Stress(QP);
        if (c.dim() != 2) Exit("TODO");

        for (int i = 0; i < felem->i_dimension(); ++i) {
          if (felem->is_zero_test(q, i)) continue;
          Tensor S_i = felem->StressTestSpace(q, i);
          VectorField V_i = felem->VelocityTestSpace(q, i);
          const double SFlux_NN = SFlux_c(S_i, Nq, Nq);
          const double SFlux_TN = SFlux_c(S_i, Tq, Nq);
          const double VFlux_N = VFlux_c(V_i, Nq);
          const double VFlux_T = VFlux_c(V_i, Tq);

          if (bnd_id == 1) {
            RHS(c(), i) -= 2 * w * (alpha3 * VFlux_c(V_bnd, Nq) * SFlux_NN
                                    + alpha4 * VFlux_c(V_bnd, Nq) * VFlux_N
                                    + alpha7 * VFlux_c(V_bnd, Tq) * SFlux_TN
                                    + alpha8 * VFlux_c(V_bnd, Tq) * VFlux_T);
          } else if (bnd_id == 2) {
            RHS(c(), i) -= 2 * w * (alpha1 * SFlux_c(S_bnd, Nq, Nq) * SFlux_NN
                                    + alpha2 * SFlux_c(S_bnd, Nq, Nq) * VFlux_N
                                    + alpha5 * SFlux_c(S_bnd, Tq, Nq) * SFlux_TN
                                    + alpha6 * SFlux_c(S_bnd, Tq, Nq) * VFlux_T);
          }
          for (int j = 0; j < felem->j_dimension() - r_c.n(); ++j) {
            if (cf_prev == cf) break;
            if (felem->is_zero(q, j)) continue;
            Tensor S_j = felem->Stress(q, j);
            VectorField V_j = felem->Velocity(q, j);
            if (c_deg == c_prev_deg) {
              M_c_prev(i, j) += w * (alpha1 * SFlux_c(S_j, Nq, Nq) * SFlux_NN
                          + alpha2 * SFlux_c(S_j, Nq, Nq) * VFlux_N
                          + alpha3 * VFlux_c(V_j, Nq) * SFlux_NN
                          + alpha4 * VFlux_c(V_j, Nq) * VFlux_N
                          + alpha5 * SFlux_c(S_j, Tq, Nq) * SFlux_TN
                          + alpha6 * SFlux_c(S_j, Tq, Nq) * VFlux_T
                          + alpha7 * VFlux_c(V_j, Tq) * SFlux_TN
                          + alpha8 * VFlux_c(V_j, Tq) * VFlux_T);
            } else {
              for (int k = 0; k < vc.size(); ++k)
                M_c_prev(i, k)+= w * (alpha1 * SFlux_c(S_j, Nq, Nq) * SFlux_NN
                            + alpha2 * SFlux_c(S_j, Nq, Nq) * VFlux_N
                            + alpha3 * VFlux_c(V_j, Nq) * SFlux_NN
                            + alpha4 * VFlux_c(V_j, Nq) * VFlux_N
                            + alpha5 * SFlux_c(S_j, Tq, Nq) * SFlux_TN
                            + alpha6 * SFlux_c(S_j, Tq, Nq) * VFlux_T
                            + alpha7 * VFlux_c(V_j, Tq) * SFlux_TN
                            + alpha8 * VFlux_c(V_j, Tq) * VFlux_T)
                       * T_Matrix[j][k];
            }
          }
          for (int j = felem->j_dimension() - r_c.n(); j < felem->j_dimension(); ++j) {
            if (felem->is_zero(q, j)) continue;
            Tensor S_j = felem->Stress(q, j);
            VectorField V_j = felem->Velocity(q, j);
            M_c(i, j - (felem->j_dimension() - r_c.n()))
                += w * (alpha1 * SFlux_c(S_j, Nq, Nq) * SFlux_NN
                        + alpha2 * SFlux_c(S_j, Nq, Nq) * VFlux_N
                        + alpha3 * VFlux_c(V_j, Nq) * SFlux_NN
                        + alpha4 * VFlux_c(V_j, Nq) * VFlux_N
                        + alpha5 * SFlux_c(S_j, Tq, Nq) * SFlux_TN
                        + alpha6 * SFlux_c(S_j, Tq, Nq) * VFlux_T
                        + alpha7 * VFlux_c(V_j, Tq) * SFlux_TN
                        + alpha8 * VFlux_c(V_j, Tq) * VFlux_T);
          }
          for (int j = 0; j < felem_1->j_dimension() - r_cf.n(); ++j) {
            if (cf_prev == cf) break;
            if (felem_1->is_zero(q1, j)) continue;
            Tensor S1_j = felem_1->Stress(q1, j);
            VectorField V1_j = felem_1->Velocity(q1, j);
            if (cf_deg == cf_prev_deg) {
              M_cf_prev(i, j) -= w
                     * (alpha1 * SFlux_cf(S1_j, Nq, Nq, bnd, bnd_id) * SFlux_NN
                        + alpha2 * SFlux_cf(S1_j, Nq, Nq, bnd, bnd_id)
                          * VFlux_N
                        + alpha3 * VFlux_cf(V1_j, Nq, bnd, bnd_id) * SFlux_NN
                        + alpha4 * VFlux_cf(V1_j, Nq, bnd, bnd_id) * VFlux_N
                        + alpha5 * SFlux_cf(S1_j, Tq, Nq, bnd, bnd_id)
                          * SFlux_TN
                        + alpha6 * SFlux_cf(S1_j, Tq, Nq, bnd, bnd_id)
                          * VFlux_T
                        + alpha7 * VFlux_cf(V1_j, Tq, bnd, bnd_id) * SFlux_TN
                        + alpha8 * VFlux_cf(V1_j, Tq, bnd, bnd_id) * VFlux_T);
            } else {
              for (int k = 0; k < vcf_prev.size(); ++k)
                M_cf_prev(i, k) -= w * (alpha1 * SFlux_cf(S1_j, Nq, Nq, bnd, bnd_id)
                            * SFlux_NN
                            + alpha2 * SFlux_cf(S1_j, Nq, Nq, bnd, bnd_id)
                              * VFlux_N
                            + alpha3 * VFlux_cf(V1_j, Nq, bnd, bnd_id) * SFlux_NN
                            + alpha4 * VFlux_cf(V1_j, Nq, bnd, bnd_id) * VFlux_N
                            + alpha5 * SFlux_cf(S1_j, Tq, Nq, bnd, bnd_id)
                              * SFlux_TN
                            + alpha6 * SFlux_cf(S1_j, Tq, Nq, bnd, bnd_id)
                              * VFlux_T
                            + alpha7 * VFlux_cf(V1_j, Tq, bnd, bnd_id) * SFlux_TN
                            + alpha8 * VFlux_cf(V1_j, Tq, bnd, bnd_id) * VFlux_T)
                       * T_Matrix_f_prev[j][k];
            }
          }
          for (int j = felem_1->j_dimension() - r_cf.n(); j < felem_1->j_dimension(); ++j) {
            if (felem_1->is_zero(q1, j)) continue;
            Tensor S1_j = felem_1->Stress(q1, j);
            VectorField V1_j = felem_1->Velocity(q1, j);
            M_cf(i, j - (felem_1->j_dimension() - r_cf.n())) -=
                w * (alpha1 * SFlux_cf(S1_j, Nq, Nq, bnd, bnd_id) * SFlux_NN
                     + alpha2 * SFlux_cf(S1_j, Nq, Nq, bnd, bnd_id) * VFlux_N
                     + alpha3 * VFlux_cf(V1_j, Nq, bnd, bnd_id) * SFlux_NN
                     + alpha4 * VFlux_cf(V1_j, Nq, bnd, bnd_id) * VFlux_N
                     + alpha5 * SFlux_cf(S1_j, Tq, Nq, bnd, bnd_id) * SFlux_TN
                     + alpha6 * SFlux_cf(S1_j, Tq, Nq, bnd, bnd_id) * VFlux_T
                     + alpha7 * VFlux_cf(V1_j, Tq, bnd, bnd_id) * SFlux_TN
                     + alpha8 * VFlux_cf(V1_j, Tq, bnd, bnd_id) * VFlux_T);
          }
        }
      }
      delete felem;
      delete felem_1;
    }
    t_face += Date() - Start_face;

  }
  t_cell.Max();
  t_face.Max();
  mout << "      cell assemble time " << t_cell << endl;
  mout << "      face assemble time " << t_face << endl;
  mout << "   assemble time " << Date() - Start << " finish assemble " << Date();
}

void STPGViscoElasticAssemble::SystemAddDoubleD(Matrix &M) const {
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    SpaceTimeViscoElasticElement elem(*disc, M, c, prob->nL());
    cell c_prev = M.find_previous_cell(c);
    row r_c = M.find_row(c());
    int c_deg = M.GetDoF().get_space_deg(*c);
    int c_prev_deg = M.GetDoF().get_space_deg(*c_prev);
    vector<vector<Scalar>> T_Matrix(M.GetDoF().get_m(c()));
    vector<Scalar> vc_diag(M.GetDoF().get_m(c()));
    vector<Scalar> vc(M.GetDoF().get_m((c_prev())));
    if (c_deg != c_prev_deg) {
      for (int i = 0; i < vc_diag.size(); ++i)
        T_Matrix[i].resize(vc.size());
      for (int i = 0; i < vc.size(); ++i) {
        for (int k = 0; k < vc.size(); ++k)
          vc[k] = 0.0;
        vc[i] = 1.;
        PX_P0(vc, vc_diag, c_prev_deg, c_deg, M);
        for (int k = 0; k < vc_diag.size(); ++k)
          T_Matrix[k][i] = vc_diag[k];
      }
    }
    DGRowEntries M_c(M, *c, *c, false);
    DGRowEntries M_c_prev(M, *c, *c_prev, true);
    const double rho = prob->Rho(c());
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);

      for (int i = 0; i < elem.i_dimension(); ++i) {
        int var_i = elem.variable(i);
        if (var_i < 2) continue;
        Tensor S_i = elem.StressTestSpace(q, i);
        for (int j = 0; j < elem.j_dimension() - r_c.n(); ++j) {
          if (c_prev == c) break;
          int var_j = elem.variable(j);
          if (var_i != var_j) continue;

          double mu = prob->Mu(c(), var_j);
          double tau = prob->Tau(c(), var_j);
          double s = 1 / (2.0 * mu);
          double lambda = prob->Lambda(c(), var_j);
          double ss = -lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));
          Tensor S_j = elem.Stress(q, j);

          double value = w / tau * (s * Frobenius(S_j, S_i) + ss * trace(S_j) * trace(S_i));

          if (c_deg == c_prev_deg) {
            M_c_prev(i, j) += value;
          } else {
            for (int k = 0; k < vc.size(); ++k) {
              M_c_prev(i, k) += value * T_Matrix[j][k];
            }
          }
        }
        for (int j = elem.j_dimension() - r_c.n(); j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          if (var_i != var_j) continue;

          double mu = prob->Mu(c(), var_j);
          double tau = prob->Tau(c(), var_j);
          double s = 1 / (2.0 * mu);
          double lambda = prob->Lambda(c(), var_j);
          double ss = -lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));

          Tensor S_j = elem.Stress(q, j);
          M_c(i, j - elem.j_dimension() + r_c.n()) += w / tau * (
              s * Frobenius(S_j, S_i) + ss * trace(S_j) * trace(S_i));
        }
      }
    }
  }
}

pair<double, double> STPGViscoElasticAssemble::DiscNorm(const Matrix &L, const Vector &u) const {
  Exit("TODO");
}

double STPGViscoElasticAssemble::L2Error(const Vector &u) const {
  Exit("TODO");
}

double STPGViscoElasticAssemble::GNError(const Vector &) const {
  mout << "WARNING ... GNError not implemented!" << endl;
  return 0.0;
}

double STPGViscoElasticAssemble::L2Norm(const Vector &) const {Exit("not implemented"); }

void STPGViscoElasticAssemble::DualRHS_linear
    (Vector &RHS, const Vector &U, Point a, Point b) const {
  RHS = 0.0;

  if (a.t() != b.t()) { // volume evaluation
    Exit("TODO");
  } else { // evaluation at some time point
    for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {
      if (c.max() != a.t()) continue;

      SpaceTimeViscoElasticElement elem(*disc, RHS, c, prob->nL());

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
        Tensor j_S;
        for (int i = 0; i < prob->Dim() - c.dim(); ++i) {
          Exit("check components")
          if (i < c.dim()) {
            COMPONENT comp = prob->GetComponents()[c.dim() + i];
            j_S[i][i] = prob->ut_dual(QP, comp, a, b);
          } else {
            COMPONENT comp = prob->GetComponents()[2 * c.dim()];
            j_S[0][1] = prob->ut_dual(QP, comp, a, b);
            j_S[1][0] = prob->ut_dual(QP, comp, a, b);
            if (c.dim() == 3) Exit("not implemented")
          }
        }

        int time_deg = U.GetDoF().get_time_deg(*c);
        double scale = 0.0;
        for (int j = time_deg * (elem.j_dimension() / (time_deg + 1));
             j < elem.j_dimension(); ++j)
          scale += elem.Velocity(q, j)[0];

        for (int j = time_deg * (elem.j_dimension() / (time_deg + 1));
             j < elem.j_dimension(); ++j) {
          Tensor S_j = elem.Stress(q, j);
          VectorField V_j = elem.Velocity(q, j);
          RHS(r_RHS)[j - elem.j_dimension() + r_RHS.n()] -=
              w * (Frobenius(j_S, S_j) + j_V * V_j) / scale;
        }
      }
    }
  }
}

void STPGViscoElasticAssemble::DualRHS_quadratic
    (Vector &RHS, const Vector &U, Point a, Point b) const {
  Exit("TODO");
}

double STPGViscoElasticAssemble::DualErrorEstimateCell
    (const cell &c, const Vector &u, const Vector &u_star) const {
  Exit("TODO");
}

double STPGViscoElasticAssemble::DualErrorEstimateJumpCell
    (const cell &c, const Vector &u, const Vector &u_star) const {
  double eta = 0.0;
  //Date Start_cell;

  SpaceTimeViscoElasticElement elem(*disc, u, c, prob->nL());
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

  /*double h = dist(c[0],c[1]);
   *   if (dist(c[0],c[2])>h)
   *       h = dist(c[0],c[2]);
   *   if (dist(c[1],c[2])>h)
   *       h = dist(c[1],c[2]);*/
  double h = max(double(dist(c[0], c[1])), double(dist(c[0], c[2])));
  h = max(h, dist(c[1], c[2]));
  //     h = sqrt(h * h + (c.max()-c.min()) * (c.max()-c.min()));
  double I = c.max() - c.min();

  double rho_1_V = 0.0;
  double rho_1_S = 0.0;
  double omega_1_V = 0.0;
  double omega_1_S = 0.0;

  double rho_2_V = 0.0;
  double rho_2_S = 0.0;
  double omega_2_V = 0.0;
  double omega_2_S = 0.0;

  int deg = u.GetDoF().get_space_deg(*c);
  int time_deg = u.GetDoF().get_time_deg(*c);

  if (c.min() != 0) {
      c_prev = u.find_cell_or_overlap_cell(f.Left());
  } else {
    vector<Point> c_nodal = u.GetDoF().GetNodalPoints(*c);
    c_nodal_value.resize(c_nodal.size());
    int cnt = 0;
    for (int ti = 0; ti < time_deg; ++ti) {
      for (int pi = 0; pi < prob->Dim(); ++pi) {
        COMPONENT comp = prob->GetComponents()[pi];
        for (int i = 0; i < u.GetDoF().NodalPointsLocal(deg); ++i) {
          c_nodal_value[cnt] = prob->ut(c_nodal[cnt].CopyWithT(0.0), comp);
          cnt++;
        }
      }
    }
  }

  int prev_deg = u.GetDoF().get_space_deg(*c_prev);
  int prev_time_deg = u.GetDoF().get_time_deg(*c_prev);
  int shift = (prev_time_deg - 1) * (r_prev.n() / prev_time_deg);

  vector<Scalar> vc(u.GetDoF().get_m(c()));
  vector<Scalar> vc_prev(u.GetDoF().get_m(c_prev()));
  vector<vector<Scalar> > T_Matrix(vc.size());

  if (deg != prev_deg) {
    for (int i = 0; i < vc.size(); ++i)
      T_Matrix[i].resize(vc_prev.size());
    for (int i = 0; i < vc_prev.size(); ++i) {
      for (int k = 0; k < vc_prev.size(); ++k)
        vc_prev[k] = 0.0;
      vc_prev[i] = 1.;
      PX_P0(vc_prev, vc, prev_deg, deg, u);
      for (int k = 0; k < vc.size(); ++k)
        T_Matrix[k][i] = vc[k];
    }
  }

  const double rho = prob->Rho(c());
  for (int q = 0; q < elem.nQ(); ++q) {
    VectorField V = 0.0;
    Tensor S = 0.0;

    double w = elem.QWeight(q);
    Point QP = elem.QPoint(q);

    for (int j = 0; j < elem.j_dimension() - r.n(); ++j) {
      int var_j = elem.variable(j);
      double mu = prob->Mu(c(), var_j);
      //double kappa = problem->Kappa( c(),var_j );
      double tau = prob->Tau(c(), var_j);
      double s = 1 / (2.0 * mu);
      //double ss = 3/kappa - 1/(6.0 * mu);
      double lambda = prob->Lambda(c(), var_j);
      double ss = -lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));

      Tensor DtS_j = elem.DtStress(q, j);
      Tensor TraceDtS_j = 0;
      TraceDtS_j[0][0] = trace(DtS_j);
      TraceDtS_j[1][1] = trace(DtS_j);
      if (c.dim() == 3)
        TraceDtS_j[2][2] = trace(DtS_j);
      VectorField DivS_j = elem.DivStress(q, j);
      VectorField DtV_j = elem.DtVelocity(q, j);
      Tensor E_j = elem.Strain(q, j);
      if (prev_deg != deg) {
        for (int k = 0; k < vc_prev.size(); ++k) {
          V += u_r_prev[k + shift] * (rho * DtV_j - DivS_j) * T_Matrix[j][k];
          S += T_Matrix[j][k] * u_r_prev[k + shift]
               * (s * DtS_j - ss * TraceDtS_j - E_j);
        }
      } else {
        if (c.min() != 0.0) {
          V += u_r_prev[j + shift] * (rho * DtV_j - DivS_j);
          S += u_r_prev[j + shift] * (s * DtS_j - ss * TraceDtS_j - E_j);
        } else {
          V += c_nodal_value[j] * (rho * DtV_j - DivS_j);
          S += c_nodal_value[j] * (s * DtS_j - ss * TraceDtS_j - E_j);
        }
      }
    }
    for (int j = elem.j_dimension() - r.n(); j < elem.j_dimension(); ++j) {
      int var_j = elem.variable(j);
      double mu = prob->Mu(c(), var_j);
      //double kappa = problem->Kappa( c(),var_j );
      double tau = prob->Tau(c(), var_j);
      double s = 1 / (2.0 * mu);
      //double ss = 3/kappa - 1/(6.0 * mu);
      double lambda = prob->Lambda(c(), var_j);
      double ss = -lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));

      Tensor DtS_j = elem.DtStress(q, j);
      Tensor TraceDtS_j = 0;
      TraceDtS_j[0][0] = trace(DtS_j);
      TraceDtS_j[1][1] = trace(DtS_j);
      if (c.dim() == 3)
        TraceDtS_j[2][2] = trace(DtS_j);
      VectorField DivS_j = elem.DivStress(q, j);
      VectorField DtV_j = elem.DtVelocity(q, j);
      Tensor E_j = elem.Strain(q, j);
      S +=
          u_r[j - elem.j_dimension() + r.n()] * (s * DtS_j - ss * TraceDtS_j - E_j);
      V += u_r[j - elem.j_dimension() + r.n()] * (rho * DtV_j - DivS_j);
    }
    rho_1_S += w * Frobenius(S, S);
    rho_1_V += w * V * V;
  }

  vector<VectorField> dual_V_projection(time_deg);
  vector<Tensor> dual_S_projection(time_deg);
  for (int l = 0; l < time_deg; ++l) {
    dual_V_projection[l] = 0.0;
    dual_S_projection[l] = 0.0;
  }

  double K_vol = 0;
  for (int q = 0; q < elem.nQ(); ++q) {
    K_vol += elem.QWeight(q);
    for (int i = 0; i < elem.i_dimension(); ++i) {
      for (int l = 0; l < time_deg; ++l) {
        dual_V_projection[l] +=
            elem.QWeight(q) * u_star_r[i % (elem.i_dimension() / time_deg)
                                       + l * elem.i_dimension() / time_deg]
            * elem.VelocityTestSpace(q, i);
        dual_S_projection[l] +=
            elem.QWeight(q) * u_star_r[i % (elem.i_dimension() / time_deg)
                                       + l * elem.i_dimension() / time_deg]
            * elem.StressTestSpace(q, i);
      }
    }
  }
  for (int l = 0; l < time_deg; ++l) {
    dual_V_projection[l] *= 1. / K_vol;
    dual_S_projection[l] *= 1. / K_vol;
  }

  const Shape &time_shape(disc->GetTimeShape(u.GetDoF().get_time_deg(*c) - 1));

  vector<double> V_Fluxes(c.Faces() - 2);
  vector<double> S_Fluxes(c.Faces() - 2);
  vector<double> V_dual_Fluxes(c.Faces() - 2);
  vector<double> S_dual_Fluxes(c.Faces() - 2);

  const double c_s = prob->C_s(c());
  const double c_p = prob->C_p(c());
  const double mu = prob->Mu(c());
  const double lambda = prob->Lambda(c());

  for (int f = 0; f < c.Faces() - 2; ++f) {
    V_Fluxes[f] = 0;
    V_dual_Fluxes[f] = 0;
    S_Fluxes[f] = 0;
    S_dual_Fluxes[f] = 0;

    cell cf = c;
    bool bnd = false;
    if (u.OnBoundary(*c, f)) bnd = false;
    else cf = u.find_neighbour_cell(c, f);

    if (bnd) continue;

    double hf = dist(c[f], c[(f + 1) % 3]);
    hf = sqrt(hf * hf + (c.max() - c.min()) * (c.max() - c.min()));

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
    const Scalar *u_star_rf = u_star(rf);

    int cf_prev_time_deg = u.GetDoF().get_time_deg(*cf_prev);

    int cf_deg = u.GetDoF().get_space_deg(*cf);
    int cf_prev_deg = u.GetDoF().get_space_deg(*cf_prev);

    int f1 = 0;
    Point f_c = u.find_face(c.Face(f))();
    Point f_cf = Origin;
    for (; f1 < cf.Faces() - 2; ++f1) {
      f_cf = u.find_face(cf.Face(f1))();
      Point pfcf = f_cf.CopyWithT(0.0);
      Point pfc = f_c.CopyWithT(0.0);
      if (pfc == pfcf) break;
    }

    SpaceTimeViscoElasticFaceElement felem(*disc, u, c, f, prob->nL());
    SpaceTimeViscoElasticFaceElement felem_1(*disc, u, cf, f1, prob->nL());
    SpaceTimeViscoElasticElement elem_1(*disc, u, cf, prob->nL());

    if (c.min() != 0) {
        cf_prev = u.find_cell_or_overlap_cell(ff.Left());
    } else {
      vector<Point> cf_nodal = u.GetDoF().GetNodalPoints(*cf);
      cf_nodal_value.resize(u.GetDoF().get_m(cf()));

      int cnt = 0;
      for (int pi = 0; pi < prob->Dim(); ++pi) {
        COMPONENT comp = prob->GetComponents()[pi];
        for (int i = 0; i < u.GetDoF().NodalPointsLocal(cf_deg); ++i) {
          cf_nodal_value[cnt] = prob->ut(cf_nodal[cnt].CopyWithT(0.0), comp);
          cnt++;
        }
      }
    }

    vector<Scalar> vcf(u.GetDoF().get_m(cf()));
    vector<Scalar> vcf_prev(u.GetDoF().get_m(cf_prev()));
    vector<vector<Scalar> > T_Matrix_f_prev(vcf.size());

    if (cf_deg != cf_prev_deg) {
      for (int i = 0; i < vcf.size(); ++i)
        T_Matrix_f_prev[i].resize(vcf_prev.size());
      for (int i = 0; i < vcf_prev.size(); ++i) {
        for (int k = 0; k < vcf_prev.size(); ++k)
          vcf_prev[k] = 0.0;
        vcf_prev[i] = 1.;
        PX_P0(vcf_prev, vcf, cf_prev_deg, cf_deg, u);
        for (int k = 0; k < vcf.size(); ++k)
          T_Matrix_f_prev[k][i] = vcf[k];
      }
    }

    double eta_v = 0.;
    int shift_1 = (cf_prev_time_deg - 1) * (rf_prev.n() / cf_prev_time_deg);

    const double c_s_f = prob->C_s(cf());
    const double c_p_f = prob->C_p(cf());
    const double mu_f = prob->Mu(cf());
    const double lambda_f = prob->Lambda(cf());
    const double alpha1 = ((2.0 * mu_f + lambda_f) * c_p)
                 / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
    const double alpha2 = ((2.0 * mu_f + lambda_f) * (2.0 * mu + lambda))
                          / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
    const double alpha3 = (mu_f * c_s) / (c_s * mu_f + c_s_f * mu);
    const double alpha4 = (mu_f * mu) / (c_s * mu_f + c_s_f * mu);
    const double alpha5 =
        (c_p_f * c_p) / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
    const double
        alpha6 = (c_p_f * (2.0 * mu + lambda))
                 / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
    const double alpha7 = (c_s_f * c_s) / (c_s * mu_f + c_s_f * mu);
    const double alpha8 = (c_s_f * mu) / (c_s * mu_f + c_s_f * mu);

    double fxI = 0;
    for (int q = 0; q < felem.nQ(); ++q) {
      int q1 = find_q_id(felem.QPoint(q), felem_1);
      VectorField V_Flux = 0.0;
      Tensor S_Flux = 0.0;

      double w = felem.QWeight(q);
      VectorField Nq = felem.QNormal(q);
      VectorField Tq = felem.QTangent(q);

      fxI += w;

      for (int j = 0; j < felem.j_dimension(); ++j) {
        if (felem.is_zero(q, j)) continue;
        Tensor S_j = 0.0;
        VectorField V_j = 0.0;
        if (j < felem.j_dimension() - r.n()) {
          if (deg != prev_deg) {
            for (int k = 0; k < vc_prev.size(); ++k) {
              S_j +=
                  T_Matrix[j][k] * u_r_prev[k + shift] * felem.Stress(q, j);
              V_j += u_r_prev[k + shift] * felem.Velocity(q, j)
                     * T_Matrix[j][k];
            }
          } else {
            S_j += u_r_prev[j + shift] * felem.Stress(q, j);
            V_j += u_r_prev[j + shift] * felem.Velocity(q, j);
          }
        } else {
          S_j += u_r[j - felem.j_dimension() + r.n()] * felem.Stress(q, j);
          V_j += u_r[j - felem.j_dimension() + r.n()] * felem.Velocity(q, j);
        }
        S_Flux +=
            double(alpha1 * VFlux_c(V_j, Nq) + alpha5 * SFlux_c(S_j, Nq, Nq))
            * Tensor(Product(Nq, Nq))
            + double(
                alpha3 * VFlux_c(V_j, Tq) + alpha7 * SFlux_c(S_j, Tq, Nq))
              * Tensor(Product(Tq, Nq));
        V_Flux -=
            double(alpha2 * VFlux_c(V_j, Nq) + alpha6 * SFlux_c(S_j, Nq, Nq)) * Nq
            + double(
                alpha4 * VFlux_c(V_j, Tq) + alpha8 * SFlux_c(S_j, Tq, Nq))
              * Tq;
      }
      for (int j = 0; j < felem_1.j_dimension(); ++j) {
        if (felem_1.is_zero(q1, j)) continue;
        Tensor S1_j = 0.0;
        VectorField V1_j = 0.0;
        if (j < felem_1.j_dimension() - rf.n()) {
          if (cf_deg != cf_prev_deg) {
            for (int k = 0; k < vcf_prev.size(); ++k) {
              S1_j += T_Matrix_f_prev[j % vcf.size()][k]
                      * u_rf_prev[k + shift_1] * felem_1.Stress(q1, j);
              V1_j +=
                  u_rf_prev[k + shift_1] * felem_1.Velocity(q1, j)
                  * T_Matrix_f_prev[j % vcf.size()][k];
            }
          } else {
            S1_j += u_rf_prev[j + shift_1] * felem_1.Stress(q1, j);
            V1_j += u_rf_prev[j + shift_1] * felem_1.Velocity(q1, j);
          }
        } else {
          S1_j +=
              u_rf[j - felem_1.j_dimension() + rf.n()] * felem_1.Stress(q1, j);
          V1_j += u_rf[j - felem_1.j_dimension() + rf.n()]
                  * felem_1.Velocity(q1, j);
        }

        S_Flux -= double(alpha1 * VFlux_cf(V1_j, Nq, bnd)
                         + alpha5 * SFlux_cf(S1_j, Nq, Nq, bnd))
                  * Tensor(Product(Nq, Nq))
                  + double(alpha3 * VFlux_cf(V1_j, Tq, bnd)
                           + alpha7 * SFlux_cf(S1_j, Tq, Nq, bnd))
                    * Tensor(Product(Tq, Nq));
        V_Flux -= double(alpha2 * VFlux_cf(V1_j, Nq, bnd)
                         + alpha6 * SFlux_cf(S1_j, Nq, Nq, bnd)) * Nq
                  + double(alpha4 * VFlux_cf(V1_j, Tq, bnd)
                           + alpha8 * SFlux_cf(S1_j, Tq, Nq, bnd)) * Tq;
      }
      V_Fluxes[f] += w * V_Flux * V_Flux;
      S_Fluxes[f] += w * Frobenius(S_Flux, S_Flux);
    }

    int time_deg_cf = u.GetDoF().get_time_deg(*cf);

    vector<VectorField> dual_V_projection_cf(time_deg_cf);
    vector<Tensor> dual_S_projection_cf(time_deg_cf);
    for (int l = 0; l < time_deg_cf; ++l) {
      dual_V_projection_cf[l] = 0.0;
      dual_S_projection_cf[l] = 0.0;
    }

    double Kf_vol = 0;

    for (int q = 0; q < elem_1.nQ(); ++q) {
      Kf_vol += elem_1.QWeight(q);
      for (int i = 0; i < elem_1.i_dimension(); ++i) {
        for (int l = 0; l < time_deg_cf; ++l) {
          dual_V_projection_cf[l] += elem_1.QWeight(q)
                                     * u_star_rf[i % (elem_1.i_dimension() / time_deg_cf)
                                                 + l * elem_1.i_dimension() / time_deg_cf]
                                     * elem_1.VelocityTestSpace(q, i);
          dual_S_projection_cf[l] += elem_1.QWeight(q)
                                     * u_star_rf[i % (elem_1.i_dimension() / time_deg_cf)
                                                 + l * elem_1.i_dimension() / time_deg_cf]
                                     * elem_1.StressTestSpace(q, i);
        }
      }
    }

    for (int l = 0; l < time_deg_cf; ++l) {
      dual_V_projection_cf[l] *= 1. / Kf_vol;
      dual_S_projection_cf[l] *= 1. / Kf_vol;
    }

    const Shape &time_shape_cf(disc->GetTimeShape(u.GetDoF().get_time_deg(*cf) - 1));

    double f_vol = fxI / I;

    auto &timeQuad = disc->GetTimeQuad(u.GetDoF().get_time_deg(*cf), true);

    for (int q = 0; q < timeQuad.size(); ++q) {
      VectorField values_dual_V = 0.0;
      VectorField values_dual_S = 0.0;
      for (int l = 0; l < time_shape.size(); ++l) {
        double timeValue = time_shape(timeQuad.QPoint(q), l);
        values_dual_V += dual_V_projection[l] * timeValue;
        values_dual_S += dual_S_projection[l] * felem.QNormal(0) * timeValue;
      }

      VectorField values_dual_V_cf = 0.0;
      VectorField values_dual_S_cf = 0.0;
      for (int l = 0; l < time_shape_cf.size(); ++l) {
        double timeValue = time_shape_cf(timeQuad.QPoint(q), l);
        values_dual_V_cf += dual_V_projection_cf[l] * timeValue;
        values_dual_S_cf += dual_S_projection_cf[l] * felem.QNormal(0) * timeValue;
      }

      V_dual_Fluxes[f] += f_vol * (values_dual_V_cf - values_dual_V)
                          * (values_dual_V_cf - values_dual_V)
                          * timeQuad.Weight(q) * I;
      S_dual_Fluxes[f] += f_vol * (values_dual_S_cf - values_dual_S)
                          * (values_dual_S_cf - values_dual_S)
                          * timeQuad.Weight(q) * I;
    }
  }

  double eta_f = 0;
  for (int f = 0; f < c.Faces() - 2; ++f) {
    omega_1_V += V_dual_Fluxes[f];
    omega_1_S += S_dual_Fluxes[f];
    eta_f += sqrt(V_Fluxes[f]) * sqrt(V_dual_Fluxes[f])
             + sqrt(S_Fluxes[f]) * sqrt(S_dual_Fluxes[f]);
  }

  //     mout<<sqrt(rho_1_P)<<" : "<<sqrt(h)<<" : "<<omega_1_P<<" : "<<sqrt(rho_1_V)<<" : "<<sqrt(h)<<" : "<<omega_1_V <<" : "<< 0.5 * eta_f<<endl;
  eta = sqrt(rho_1_S) * sqrt(h) * sqrt(omega_1_S)
        + sqrt(rho_1_V) * sqrt(h) * sqrt(omega_1_V) + 0.5 * eta_f;
  return eta;
}

double STPGViscoElasticAssemble::ErrorEstimateCell
    (const cell &c, const Vector &u) const {
  Exit("Not implemented");
}

double STPGViscoElasticAssemble::Goal_Functional
    (const Vector &u, Point a, Point b) const {
  Scalar e = 0.0;

  if (a.t() != b.t()) {
    Exit("ToDo");
  } else { // evaluation at some time point
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != a.t()) continue;
      SpaceTimeViscoElasticElement elem(*disc, u, c, prob->nL());
      row r = u.find_row(c());
      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        Point QP = elem.QPoint(q).CopyWithT(0.0);
        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;
        double w = elem.SpaceQWeight(q);
        Exit("add ut_dual_Velocity and ut_dual_Stress to Problem-Class")
        VectorField j_V = (prob->ut_dual(QP, COMPONENT::V_X, a, b),
                        prob->ut_dual(QP, COMPONENT::V_Y, a, b));
        Tensor j_S(prob->ut_dual(QP, COMPONENT::S0_x, a, b), prob->ut_dual(QP, COMPONENT::S0_xy, a, b),
                   prob->ut_dual(QP, COMPONENT::S0_xy, a, b), prob->ut_dual(QP, COMPONENT::S0_y, a, b));
        if (c.dim() == 3) Exit("not implemented")

        int time_deg = u.GetDoF().get_time_deg(*c);
        double scale = 0.0;
        const int startj = time_deg * (elem.j_dimension() / (time_deg + 1));
        for (int j = startj; j < elem.j_dimension(); ++j) {
          scale += elem.Velocity(q, j)[0];
        }
        Tensor S = 0.0;
        VectorField V = 0.0;
        for (int j = startj; j < elem.j_dimension(); ++j) {
          Tensor S_j = elem.Stress(q, j);
          VectorField V_j = elem.Velocity(q, j);
          S += u(r, j - elem.j_dimension() + r.n()) * S_j;
          V += u(r, j - elem.j_dimension() + r.n()) * V_j;
        }
        e += w * (Frobenius(j_S, S) + j_V * V) / scale;
      }
    }
  }
  measure(u);
  return abs(PPM->SumOnCommSplit(e, u.CommSplit()));
}

double STPGViscoElasticAssemble::Energy(const Vector &u, Point a, Point b) const {
  Exit("TODO");
}

void STPGViscoElasticAssemble::measure(const Vector &u) const {
  int N = -1;
  Config::Get("receiver_count", N, -1);

  if (N == -1) return;

  if (N <= 0) {
    mout << "******************************************************************* *"
         << endl;
    mout << "* ERROR: invalid number of recievers - WILL NOT STORE MEASUREMENTS! *"
         << endl;
    mout << "*********************************************************************"
         << endl;
    return;
  }

  Date Start;
  mout << "Start Measure at " << Start << endl;

  Point r_first = Origin;
  Config::Get("first_receiver", r_first);
  Point r_last = Origin;
  Config::Get("last_receiver", r_last);
  Point time_config = Origin;
  Config::Get("receive_times", time_config);
  list<Point> receiver_positions;
  for (int i = 0; i < N; ++i)
    receiver_positions.push_back(r_first + i * (r_last - r_first) / (N - 1));
  FWIObservationSpecification
      spec(time_config[0], time_config[1], time_config[2], receiver_positions);

  //const char* TMP_SPEC_FILE = "test_spec.txt";
  //spec.writeToFile(TMP_SPEC_FILE);
  //FWIObservationSpecification tmp_spec(TMP_SPEC_FILE);
  //FWISeismogram s(tmp_spec, STATE_COMPS);
  //auto& spec = s.getSpecification();



  double radius = 0.3;
  Config::Get("measurement_radius", radius);
  int STATE_COMPS = prob->Dim();
  FWISeismogram s(spec, STATE_COMPS);

  for (auto it = spec.begin(); it != spec.end(); ++it) {
    Point rec(*it);
    for (int i = 0; i < spec.nSteps(); ++i) {
      double t = spec.step(i);
      rec = rec.CopyWithT(t);
      FWIMeasurementKernel_L1 phi(u.dim(), rec, radius);

      RVector X(0.0, STATE_COMPS);

      for (cell c = u.cells(); c != u.cells_end(); ++c) {

        if (!cellAndSphereOverlap(c, rec, radius)) continue;

        SpaceTimeViscoElasticElement elem(*disc, u, c, prob->nL());
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

        int prev_time_deg = u.GetDoF().get_time_deg(*c_prev);

        vector<Scalar> vc(u.GetDoF().get_m(c()));
        vector<Scalar> vc_prev(u.GetDoF().get_m(c_prev()));

        vector<vector<Scalar> > T_Matrix(vc.size());
        for (int i = 0; i < vc.size(); ++i)
          T_Matrix[i].resize(vc_prev.size());
        for (int i = 0; i < vc_prev.size(); ++i) {
          for (int k = 0; k < vc_prev.size(); ++k)
            vc_prev[k] = 0.0;
          vc_prev[i] = 1.;
          PX_P0(vc_prev, vc, u.GetDoF().get_space_deg(*c_prev), u.GetDoF().get_space_deg(*c), u);
          for (int k = 0; k < vc.size(); ++k)
            T_Matrix[k][i] = vc[k];
        }

        for (int q = 0; q < elem.nQ(); ++q) {
          const Point QP = elem.QPoint(q);
          const double w = elem.QWeight(q);
          const double phiX = phi(QP);

          if (phiX == 0) continue;

          for (int j = 0; j < elem.j_dimension(); ++j) {

            VectorField V_j = elem.Velocity(q, j);
            Tensor S_j = elem.Stress(q, j);

            if (j < elem.j_dimension() - r.n()) {
              if (c.min() == 0) {
                for (int i = 0; i < c.dim(); ++i) {
                  Exit("check components!")
                  COMPONENT comp = prob->GetComponents()[i];
                  V_j[i] *= prob->ut(c_nodal[j], comp) * w * phiX;
                  COMPONENT comp2 = prob->GetComponents()[c.dim() + i];
                  S_j[i][i] *= prob->ut(c_nodal[j], comp2) * w * phiX;
                  for (int ii = 0; ii < i; ++ii) {
                    COMPONENT comp = prob->GetComponents()[2 * c.dim() + ii + i - 1];
                    S_j[i][ii] *= prob->ut(c_nodal[j], comp) * w * phiX;
                    S_j[ii][i] *= prob->ut(c_nodal[j], comp) * w * phiX;
                  }
                }
              } else
                for (int k = 0; k < vc_prev.size(); ++k) {
                  const int idx = j + (prev_time_deg - 1) * (r_prev.n() / prev_time_deg);
                  const double tmp = T_Matrix[j][k] * u(r_prev, idx) * w * phiX;
                  V_j *= tmp;
                  S_j *= tmp;
                }
            } else {
              const double tmp = u(r, j - elem.j_dimension() + r.n()) * w * phiX;
              V_j *= tmp;
              S_j *= tmp;
            }

            for (int i = 0; i < c.dim(); ++i)
              X[i] += V_j[i];
            if (elem.variable(j))
              for (int i = 0; i < c.dim(); ++i) {
                X[c.dim() + (elem.variable(j) - 1) * elem.td() + i] += S_j[i][i];
                for (int ii = 0; ii < i; ++ii) {
                  int XXX = c.dim() + (elem.variable(j) - 1) * elem.td()
                            + c.dim() + ii + i - 1;
                  X[XXX] += S_j[i][ii];
                }
              }
          }
        }
      }

      rec = rec.CopyWithT(0.0);
      s.addMeasurement(rec, t, X);
    }
  }

  mout << "     collect measure" << endl;
  s.collect();     // collects seismogram on master
  //s.broadcast(); // broadcasts to all processes

  int spaceDeg(0);
  Config::Get("degree", spaceDeg);
  int timeDeg(1);
  Config::Get("time_degree", timeDeg);
  int level = 0;
  Config::Get("level", level);
  const int InvDtCell = 1.0 / (u.cells().max() - u.cells().min());

  const string sname = "seismogram_lvl_" + to_string(level)
                       + "_sDeg_" + to_string(spaceDeg)
                       + "_tDeg_" + to_string(timeDeg)
                       + "_dtCell_1:" + to_string(InvDtCell);

  mout << "     write measure to " << sname << endl;
  s.writeToFile(sname);
  vout(1) << "TIME:        measuring time: " << Date() - Start << endl;
}

