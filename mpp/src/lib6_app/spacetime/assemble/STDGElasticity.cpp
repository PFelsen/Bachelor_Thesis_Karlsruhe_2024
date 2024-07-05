#include "STDGElasticity.hpp"
#include "GMRES.hpp"
#include "TimeDate.hpp"

using namespace std;


void TDGElasticityAssemble::System(Matrix &M, Vector &RHS) const {
  Date Start;
  Time t_cell;
  Time t_face;
  M = 0;
  Vector u_0(0.0, RHS);
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    cell c_prev = M.find_previous_cell(c);
    if (c_prev != c) continue;
    Date Start_cell;
    SpaceTimeElasticityElement elem(*disc, M, c);
    row r_c = M.find_row(c());
    int c_deg = RHS.GetDoF().get_space_deg(*c);
    DGRowEntries M_c(M, *c, *c, false);
    double rho = prob->Rho(c());
    double mu = prob->Mu(c());
    double lambda = prob->Lambda(c());
    double s = 1 / (2.0 * mu);
    double ss = lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      for (int i = 0; i < elem.i_dimension(); ++i) {
        VectorField V_i = elem.VelocityTestSpace(q, i);
        Tensor S_i = elem.StressTestSpace(q, i);
        for (int j = 0; j < elem.j_dimension() - r_c.n(); ++j) {
          VectorField DtV_j = elem.DtVelocity(q, j);
          Tensor DtS_j = elem.DtStress(q, j);
          VectorField DivS_j = elem.DivStress(q, j);
          Tensor E_j = elem.Strain(q, j);
          M_c(i, j) += w * ((rho * (DtV_j * V_i)
                             + (s * Frobenius(S_i, DtS_j)
                                - ss * trace(S_i) * trace(DtS_j)))
                            - ((V_i * DivS_j) + Frobenius(S_i, E_j)));
        }
      }
    }

    int time_deg = RHS.GetDoF().get_time_deg(*c);
    vector<Point> c_nodal = RHS.GetDoF().GetNodalPoints(*c);
    int cnt = 0;
    for (int ti = 0; ti < time_deg; ++ti) {
      for (int pi = 0; pi < prob->Dim(); ++pi) {
        COMPONENT comp = prob->GetComponents()[pi];
        for (int i = 0; i < RHS.GetDoF().NodalPointsLocal(c_deg); ++i) {
          u_0(r_c, cnt) = prob->ut(c_nodal[cnt].CopyWithT(0.0), comp);
          cnt++;
        }
      }
    }

    t_cell += Date() - Start_cell;
    Date Start_face;
    for (int f = 0; f < c.Faces() - 2; ++f) {
      cell cf = c;
      bool bnd = M.OnBoundary(*c, f);
      bnd_face bface = M.find_bnd_face(c.Face(f));
      int bnd_id = (bface == M.bnd_faces_end()) ? -1 : bface.Part();

      if (!bnd) cf = M.find_neighbour_cell(c, f);
      SpaceTimeElasticityFaceElement felem(*disc, M, c, f);
      int cf_deg = RHS.GetDoF().get_space_deg(*cf);
      int f1 = M.find_neighbour_face_id(c.Face(f), cf);

//   different code for cf = c !

      SpaceTimeElasticityFaceElement felem_1(*disc, M, cf, f1);
      DGRowEntries M_cf(M, *c, *cf, false);
      row r_cf = M.find_row(cf());
      double c_s = prob->C_s(c());
      double c_p = prob->C_p(c());
      double c_s_f = prob->C_s(cf());
      double c_p_f = prob->C_p(cf());
      double mu_f = prob->Mu(cf());
      double lambda_f = prob->Lambda(cf());

      const double tmpVar =
          1.0 / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));

      double alpha1 = ((2.0 * mu_f + lambda_f) * c_p) * tmpVar;
      double alpha2 = ((2.0 * mu_f + lambda_f) * (2.0 * mu + lambda)) * tmpVar;
      double alpha3 = (mu_f * c_s) / (c_s * mu_f + c_s_f * mu);
      double alpha4 = (mu_f * mu) / (c_s * mu_f + c_s_f * mu);
      double alpha5 = (c_p_f * c_p) * tmpVar;
      double alpha6 = (c_p_f * (2.0 * mu + lambda)) * tmpVar;
      double alpha7 = (c_s_f * c_s) / (c_s * mu_f + c_s_f * mu);
      double alpha8 = (c_s_f * mu) / (c_s * mu_f + c_s_f * mu);
      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        double w = felem.QWeight(q);
        VectorField Nq = felem.QNormal(q);
        VectorField Tq = felem.QTangent(q);
        for (int i = 0; i < felem.i_dimension(); ++i) {
          Tensor S_i = felem.StressTestSpace(q, i);
          VectorField V_i = felem.VelocityTestSpace(q, i);
          double SFlux_NN = SFlux_c(S_i, Nq, Nq);
          double SFlux_TN = SFlux_c(S_i, Tq, Nq);
          for (int j = 0; j < felem.j_dimension() - r_c.n(); ++j) {
            Tensor S_j = felem.Stress(q, j);
            VectorField V_j = felem.Velocity(q, j);
            M_c(i, j) += w * (alpha1 * SFlux_NN * VFlux_c(V_j, Nq)
                              + alpha3 * SFlux_TN * VFlux_c(V_j, Tq)
                              + alpha2 * VFlux_c(V_i, Nq) * VFlux_c(V_j, Nq)
                              + alpha4 * VFlux_c(V_i, Tq) * VFlux_c(V_j, Tq)
                              + alpha5 * SFlux_NN * SFlux_c(S_j, Nq, Nq)
                              + alpha7 * SFlux_TN * SFlux_c(S_j, Tq, Nq)
                              + alpha6 * VFlux_c(V_i, Nq) * SFlux_c(S_j, Nq, Nq)
                              + alpha8 * VFlux_c(V_i, Tq) * SFlux_c(S_j, Tq, Nq));
          }
          for (int j = 0; j < felem_1.j_dimension() - r_cf.n(); ++j) {
            Tensor S1_j = felem_1.Stress(q1, j);
            VectorField V1_j = felem_1.Velocity(q1, j);
            M_cf(i, j) -= w
                          * (alpha1 * SFlux_NN * VFlux_cf(V1_j, Nq, bnd)
                             + alpha3 * SFlux_TN * VFlux_cf(V1_j, Tq, bnd)
                             + alpha2 * VFlux_c(V_i, Nq) * VFlux_cf(V1_j, Nq, bnd)
                             + alpha4 * VFlux_c(V_i, Tq) * VFlux_cf(V1_j, Tq, bnd)
                             + alpha5 * SFlux_NN * SFlux_cf(S1_j, Nq, Nq, bnd)
                             + alpha7 * SFlux_TN * SFlux_cf(S1_j, Tq, Nq, bnd)
                             + alpha6 * VFlux_c(V_i, Nq) * SFlux_cf(S1_j, Nq, Nq, bnd)
                             + alpha8 * VFlux_c(V_i, Tq)
                               * SFlux_cf(S1_j, Tq, Nq, bnd));
          }
        }
      }
    }
    t_face += Date() - Start_face;
  }
  u_0.Accumulate();
  Vector Bu_0 = M * u_0;
  RHS -= Bu_0;

  M = 0;
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    Date Start_cell;
    SpaceTimeElasticityElement elem(*disc, M, c);
    cell c_prev = M.find_previous_cell(c);
    row r_c = M.find_row(c());
    int c_deg = RHS.GetDoF().get_space_deg(*c);
    int c_prev_deg = RHS.GetDoF().get_space_deg(*c_prev);
    vector<vector<Scalar> > T_Matrix(RHS.GetDoF().get_m(c()));
    vector<Scalar> vc_diag(RHS.GetDoF().get_m(c()));
    vector<Scalar> vc(RHS.GetDoF().get_m(c_prev()));
    if (c_deg != c_prev_deg) {
      for (int i = 0; i < vc_diag.size(); ++i)
        T_Matrix[i].resize(vc.size());

      for (int i = 0; i < vc.size(); ++i) {
        vc = ft::create_unit_vec(vc.size(), i);
        PX_P0(vc, vc_diag, c_prev_deg, c_deg, RHS);

        for (int k = 0; k < vc_diag.size(); ++k)
          T_Matrix[k][i] = vc_diag[k];
      }
    }
    DGRowEntries M_c(M, *c, *c, false);
    DGRowEntries M_c_prev(M, *c, *c_prev, true);

    double rho = prob->Rho(c());
    double mu = prob->Mu(c());
    double lambda = prob->Lambda(c());
    double s = 1 / (2.0 * mu);
    double ss = lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);

      Tensor S_rhs = prob->F_Stress(QP);
      VectorField V_rhs = prob->F_Velocity(QP);
      if (c.dim() != 2) Exit("TODO");

      for (int i = 0; i < elem.i_dimension(); ++i) {
        VectorField V_i = elem.VelocityTestSpace(q, i);
        Tensor S_i = elem.StressTestSpace(q, i);

        RHS(c(), i) += w * (Frobenius(S_i, S_rhs) + V_i * V_rhs);

        for (int j = 0; j < elem.j_dimension() - r_c.n(); ++j) {
          if (c_prev == c) break;
          VectorField DtV_j = elem.DtVelocity(q, j);
          Tensor DtS_j = elem.DtStress(q, j);
          VectorField DivS_j = elem.DivStress(q, j);
          Tensor E_j = elem.Strain(q, j);
          if (c_deg == c_prev_deg) {
            M_c_prev(i, j) += w * ((rho * (DtV_j * V_i)
                                    + (s * Frobenius(S_i, DtS_j)
                                       - ss * trace(S_i) * trace(DtS_j)))
                                   - ((V_i * DivS_j) + Frobenius(S_i, E_j)));
          } else {
            for (int k = 0; k < vc.size(); ++k)
              M_c_prev(i, k) +=
                  w * ((rho * (DtV_j * V_i)
                        + (s * Frobenius(S_i, DtS_j)
                           - ss * trace(S_i) * trace(DtS_j)))
                       - ((V_i * DivS_j) + Frobenius(S_i, E_j)))
                  * T_Matrix[j][k];
          }
        }
        for (int j = elem.j_dimension() - r_c.n(); j < elem.j_dimension(); ++j) {
          VectorField DtV_j = elem.DtVelocity(q, j);
          Tensor DtS_j = elem.DtStress(q, j);
          VectorField DivS_j = elem.DivStress(q, j);
          Tensor E_j = elem.Strain(q, j);
          M_c(i, j - elem.j_dimension() + r_c.n()) +=
              w * ((rho * (DtV_j * V_i)
                    + (s * Frobenius(S_i, DtS_j)
                       - ss * trace(S_i) * trace(DtS_j)))
                   - ((V_i * DivS_j) + Frobenius(S_i, E_j)));
        }
      }
    }
    t_cell += Date() - Start_cell;
    Date Start_face;
    for (int f = 0; f < c.Faces() - 2; ++f) {
      cell cf = c;
      bool bnd = M.OnBoundary(*c, f);
      int bnd_id = -1;
      if (bnd) {
        bnd_face bface = M.find_bnd_face(c.Face(f));
        bnd_id = (bface == M.bnd_faces_end()) ? -1 : bface.Part();
      } else cf = M.find_neighbour_cell(c, f);
      SpaceTimeElasticityFaceElement *felem = new SpaceTimeElasticityFaceElement(*disc, M, c, f);
      int cf_deg = RHS.GetDoF().get_space_deg(*cf);
      cell cf_prev = M.find_previous_cell(cf);
      int f1 = M.find_neighbour_face_id(c.Face(f), cf);

      SpaceTimeElasticityFaceElement *felem_1 = new SpaceTimeElasticityFaceElement(*disc, M, cf,
                                                                                   f1);
      int cf_prev_deg = RHS.GetDoF().get_space_deg(*cf_prev);

      vector<vector<Scalar> > T_Matrix_f_prev(RHS.GetDoF().get_m(cf()));
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
      double c_s = prob->C_s(c());
      double c_p = prob->C_p(c());
      double c_s_f = prob->C_s(cf());
      double c_p_f = prob->C_p(cf());
      double mu_f = prob->Mu(cf());
      double lambda_f = prob->Lambda(cf());

      const double tmpVar = 1.0 / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));

      double alpha1 = ((2.0 * mu_f + lambda_f) * c_p) * tmpVar;
      double alpha2 = ((2.0 * mu_f + lambda_f) * (2.0 * mu + lambda)) * tmpVar;
      double alpha3 = (mu_f * c_s) / (c_s * mu_f + c_s_f * mu);
      double alpha4 = (mu_f * mu) / (c_s * mu_f + c_s_f * mu);
      double alpha5 = (c_p_f * c_p) * tmpVar;
      double alpha6 = (c_p_f * (2.0 * mu + lambda)) * tmpVar;
      double alpha7 = (c_s_f * c_s) / (c_s * mu_f + c_s_f * mu);
      double alpha8 = (c_s_f * mu) / (c_s * mu_f + c_s_f * mu);

      for (int q = 0; q < felem->nQ(); ++q) {
        int q1 = find_q_id(felem->QPoint(q), *felem_1);
        double w = felem->QWeight(q);
        VectorField Nq = felem->QNormal(q);
        VectorField Tq = felem->QTangent(q);
        Point QP = felem->QPoint(q);

        Tensor S_bnd = prob->ut_Stress(QP);
        VectorField V_bnd = prob->ut_Velocity(QP);

        for (int i = 0; i < felem->i_dimension(); ++i) {
          Tensor S_i = felem->StressTestSpace(q, i);
          VectorField V_i = felem->VelocityTestSpace(q, i);
          double SFlux_NN = SFlux_c(S_i, Nq, Nq);
          double SFlux_TN = SFlux_c(S_i, Tq, Nq);

          if (bnd_id == 1)
            RHS(c(), i) -= 2 * w * (alpha1 * SFlux_NN * VFlux_c(V_bnd, Nq)
                                    + alpha2 * VFlux_c(V_i, Nq) * VFlux_c(V_bnd, Nq)
                                    + alpha3 * SFlux_TN * VFlux_c(V_bnd, Tq)
                                    + alpha4 * VFlux_c(V_i, Tq) * VFlux_c(V_bnd, Tq));
          else if (bnd_id == 2)
            RHS(c(), i) -= 2 * w * (alpha5 * SFlux_NN * SFlux_c(S_bnd, Nq, Nq)
                                    + alpha6 * VFlux_c(V_i, Nq) * SFlux_c(S_bnd, Nq, Nq)
                                    + alpha7 * SFlux_TN * SFlux_c(S_bnd, Tq, Nq)
                                    + alpha8 * VFlux_c(V_i, Tq) * SFlux_c(S_bnd, Tq, Nq));

          for (int j = 0; j < felem->j_dimension() - r_c.n(); ++j) {
            if (cf_prev == cf) break;
            Tensor S_j = felem->Stress(q, j);
            VectorField V_j = felem->Velocity(q, j);
            if (c_deg == c_prev_deg) {
              M_c_prev(i, j) +=
                  w * (alpha1 * SFlux_NN * VFlux_c(V_j, Nq)
                       + alpha3 * SFlux_TN * VFlux_c(V_j, Tq)
                       + alpha2 * VFlux_c(V_i, Nq) * VFlux_c(V_j, Nq)
                       + alpha4 * VFlux_c(V_i, Tq) * VFlux_c(V_j, Tq)
                       + alpha5 * SFlux_NN * SFlux_c(S_j, Nq, Nq)
                       + alpha7 * SFlux_TN * SFlux_c(S_j, Tq, Nq)
                       + alpha6 * VFlux_c(V_i, Nq) * SFlux_c(S_j, Nq, Nq)
                       + alpha8 * VFlux_c(V_i, Tq) * SFlux_c(S_j, Tq, Nq));
            } else {
              for (int k = 0; k < vc.size(); ++k)
                M_c_prev(i, k) +=
                    w * (alpha1 * SFlux_NN * VFlux_c(V_j, Nq)
                         + alpha3 * SFlux_TN * VFlux_c(V_j, Tq)
                         + alpha2 * VFlux_c(V_i, Nq) * VFlux_c(V_j, Nq)
                         + alpha4 * VFlux_c(V_i, Tq) * VFlux_c(V_j, Tq)
                         + alpha5 * SFlux_NN * SFlux_c(S_j, Nq, Nq)
                         + alpha7 * SFlux_TN * SFlux_c(S_j, Tq, Nq)
                         + alpha6 * VFlux_c(V_i, Nq) * SFlux_c(S_j, Nq, Nq)
                         + alpha8 * VFlux_c(V_i, Tq)
                           * SFlux_c(S_j, Tq, Nq)) * T_Matrix[j][k];
            }
          }
          for (int j = felem->j_dimension() - r_c.n(); j < felem->j_dimension();
               ++j) {
            Tensor S_j = felem->Stress(q, j);
            VectorField V_j = felem->Velocity(q, j);
            M_c(i, j - (felem->j_dimension() - r_c.n())) +=
                w * (alpha1 * SFlux_NN * VFlux_c(V_j, Nq)
                     + alpha3 * SFlux_TN * VFlux_c(V_j, Tq)
                     + alpha2 * VFlux_c(V_i, Nq) * VFlux_c(V_j, Nq)
                     + alpha4 * VFlux_c(V_i, Tq) * VFlux_c(V_j, Tq)
                     + alpha5 * SFlux_NN * SFlux_c(S_j, Nq, Nq)
                     + alpha7 * SFlux_TN * SFlux_c(S_j, Tq, Nq)
                     + alpha6 * VFlux_c(V_i, Nq) * SFlux_c(S_j, Nq, Nq)
                     + alpha8 * VFlux_c(V_i, Tq) * SFlux_c(S_j, Tq, Nq));
          }
          for (int j = 0; j < felem_1->j_dimension() - r_cf.n(); ++j) {
            if (cf_prev == cf) break;
            Tensor S1_j = felem_1->Stress(q1, j);
            VectorField V1_j = felem_1->Velocity(q1, j);
            if (cf_deg == cf_prev_deg) {
              M_cf_prev(i, j) -=
                  w * (alpha1 * SFlux_NN * VFlux_cf(V1_j, Nq, bnd)
                       + alpha3 * SFlux_TN * VFlux_cf(V1_j, Tq, bnd)
                       + alpha2 * VFlux_c(V_i, Nq) * VFlux_cf(V1_j, Nq, bnd)
                       + alpha4 * VFlux_c(V_i, Tq) * VFlux_cf(V1_j, Tq, bnd)
                       + alpha5 * SFlux_NN * SFlux_cf(S1_j, Nq, Nq, bnd)
                       + alpha7 * SFlux_TN * SFlux_cf(S1_j, Tq, Nq, bnd)
                       + alpha6 * VFlux_c(V_i, Nq)
                         * SFlux_cf(S1_j, Nq, Nq, bnd)
                       + alpha8 * VFlux_c(V_i, Tq)
                         * SFlux_cf(S1_j, Tq, Nq, bnd));
            } else {
              for (int k = 0; k < vcf_prev.size(); ++k)
                M_cf_prev(i, k) -=
                    w * (alpha1 * SFlux_NN * VFlux_cf(V1_j, Nq, bnd)
                         + alpha3 * SFlux_TN * VFlux_cf(V1_j, Tq, bnd)
                         + alpha2 * VFlux_c(V_i, Nq)
                           * VFlux_cf(V1_j, Nq, bnd)
                         + alpha4 * VFlux_c(V_i, Tq)
                           * VFlux_cf(V1_j, Tq, bnd)
                         + alpha5 * SFlux_NN * SFlux_cf(S1_j, Nq, Nq, bnd)
                         + alpha7 * SFlux_TN * SFlux_cf(S1_j, Tq, Nq, bnd)
                         + alpha6 * VFlux_c(V_i, Nq)
                           * SFlux_cf(S1_j, Nq, Nq, bnd)
                         + alpha8 * VFlux_c(V_i, Tq)
                           * SFlux_cf(S1_j, Tq, Nq, bnd))
                    * T_Matrix_f_prev[j][k];
            }
          }
          for (int j = felem_1->j_dimension() - r_cf.n();
               j < felem_1->j_dimension(); ++j) {
            Tensor S1_j = felem_1->Stress(q1, j);
            VectorField V1_j = felem_1->Velocity(q1, j);
            M_cf(i, j - (felem_1->j_dimension() - r_cf.n())) -= w
                                                                * (alpha1 * SFlux_NN *
                                                                   VFlux_cf(V1_j, Nq, bnd)
                                                                   + alpha3 * SFlux_TN *
                                                                     VFlux_cf(V1_j, Tq, bnd)
                                                                   + alpha2 * VFlux_c(V_i, Nq) *
                                                                     VFlux_cf(V1_j, Nq, bnd)
                                                                   + alpha4 * VFlux_c(V_i, Tq) *
                                                                     VFlux_cf(V1_j, Tq, bnd)
                                                                   + alpha5 * SFlux_NN *
                                                                     SFlux_cf(S1_j, Nq, Nq, bnd)
                                                                   + alpha7 * SFlux_TN *
                                                                     SFlux_cf(S1_j, Tq, Nq, bnd)
                                                                   + alpha6 * VFlux_c(V_i, Nq) *
                                                                     SFlux_cf(S1_j, Nq, Nq, bnd)
                                                                   + alpha8 * VFlux_c(V_i, Tq)
                                                                     * SFlux_cf(S1_j, Tq, Nq, bnd));
          }
        }
      }
      delete felem;
      delete felem_1;
    }
    t_face += Date() - Start_face;
  }
  int dof_sum = 0;
  for (row r = RHS.rows(); r != RHS.rows_end(); ++r)
    if (RHS.find_cell(r()) != RHS.cells_end())
      dof_sum += r.n();
  vout(2) << "   Problem Size: " << PPM->SumOnCommSplit(dof_sum, RHS.CommSplit()) << endl;;
  vout(2) << "      cell assemble time " << t_cell << endl;
  vout(2) << "      face assemble time " << t_face << endl;
  vout(1) << "   assemble time " << Date() - Start << endl;
}

std::pair<double, double> TDGElasticityAssemble::DiscNorm(const Matrix &L, const Vector &u) const {
  Matrix M(u);
  M = 0;
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    SpaceTimeElasticityElement elem(*disc, M, c);
    cell c_prev = M.find_previous_cell(c);
    row r_c = M.find_row(c());
    int c_deg = u.GetDoF().get_space_deg(*c);
    int c_prev_deg = u.GetDoF().get_space_deg(*c_prev);
    vector<vector<Scalar> > T_Matrix(u.GetDoF().get_m(c()));
    vector<Scalar> vc_diag(u.GetDoF().get_m(c()));
    vector<Scalar> vc(u.GetDoF().get_m(c_prev()));
    if (c_deg != c_prev_deg) {
      for (int i = 0; i < vc_diag.size(); ++i)
        T_Matrix[i].resize(vc.size());

      for (int i = 0; i < vc.size(); ++i) {
        for (int k = 0; k < vc.size(); ++k)
          vc[k] = 0.0;
        vc[i] = 1.;
        PX_P0(vc, vc_diag, c_prev_deg, c_deg, u);

        for (int k = 0; k < vc_diag.size(); ++k)
          T_Matrix[k][i] = vc_diag[k];
      }
    }
    DGRowEntries M_c(M, *c, *c, false);
    DGRowEntries M_c_prev(M, *c, *c_prev, true);
    double rho = prob->Rho(c());
    double mu = prob->Mu(c());
    double lambda = prob->Lambda(c());
    double s = 1 / (2.0 * mu);
    double ss = lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      for (int i = 0; i < elem.i_dimension(); ++i) {
        VectorField V_i = elem.VelocityTestSpace(q, i);
        Tensor S_i = elem.StressTestSpace(q, i);
        for (int j = 0; j < elem.j_dimension() - r_c.n(); ++j) {
          if (c_prev == c) break;
          VectorField DtV_j = elem.DtVelocity(q, j);
          Tensor DtS_j = elem.DtStress(q, j);
          VectorField DivS_j = elem.DivStress(q, j);
          Tensor E_j = elem.Strain(q, j);
          if (c_deg == c_prev_deg) {
            M_c_prev(i, j) += w
                              * ((rho * (DtV_j * V_i) + (s * Frobenius(S_i, DtS_j)
                                                         - ss * trace(S_i) * trace(DtS_j)))
                                 - ((V_i * DivS_j) + Frobenius(S_i, E_j)));
          } else {
            for (int k = 0; k < vc.size(); ++k)
              M_c_prev(i, k) +=
                  w * ((rho * (DtV_j * V_i)
                        + (s * Frobenius(S_i, DtS_j)
                           - ss * trace(S_i) * trace(DtS_j)))
                       - ((V_i * DivS_j) + Frobenius(S_i, E_j)))
                  * T_Matrix[j][k];
          }
        }
        for (int j = elem.j_dimension() - r_c.n(); j < elem.j_dimension(); ++j) {
          VectorField DtV_j = elem.DtVelocity(q, j);
          Tensor DtS_j = elem.DtStress(q, j);
          VectorField DivS_j = elem.DivStress(q, j);
          Tensor E_j = elem.Strain(q, j);
          M_c(i, j - elem.j_dimension() + r_c.n()) +=
              w * ((rho * (DtV_j * V_i) + (s * Frobenius(S_i, DtS_j)
                                           - ss * trace(S_i) * trace(DtS_j)))
                   - ((V_i * DivS_j) + Frobenius(S_i, E_j)));
        }
      }
    }
  }

  Vector Mu(M * u);
  Vector Lu(L * u);

  Preconditioner *PC = GetPC("PointBlockGaussSeidel");
  LinearSolver S(PC);
  S(M);

  Vector mLu(u);
  mLu = S * Lu;

  double discNormW = u * Mu;
  double discNormV = Lu * mLu;

  mout << "discrete norm ||u-U||_W  = " << scientific //<< setprecision(10)
       << sqrt(discNormW) << endl;
  mout << "discrete norm ||u-U||_V  = " << scientific //<< setprecision(10)
       << sqrt(discNormW + discNormV) << endl;

  return pair<double, double>(std::sqrt(discNormW), std::sqrt(discNormW + discNormV));
}

double TDGElasticityAssemble::L2Error(Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeElasticityElement elem(*disc, u, c);
    row r = u.find_row(c());

    cell c_prev = c;
    vector<Scalar> c_nodal_value;

    int deg = u.GetDoF().get_space_deg(*c);
    int time_deg = u.GetDoF().get_time_deg(*c);

    if (c.min() != 0) {
      face ff = u.find_face(c.Face(c.Faces() - 2));
      c_prev = u.find_cell(ff.Left());
      if (c_prev == u.cells_end())
        c_prev = u.find_overlap_cell(ff.Left());
    } else {
      vector<Point> c_nodal = u.GetDoF().GetNodalPoints(*c);
      c_nodal_value.resize(u.GetDoF().get_m(c()));

      int cnt = 0;
      for (int pi = 0; pi < prob->Dim(); ++pi) {
        COMPONENT comp = prob->GetComponents()[pi];
        for (int i = 0; i < u.GetDoF().NodalPointsLocal(deg); ++i) {
          c_nodal_value[cnt] = prob->ut(c_nodal[cnt].CopyWithT(0.0), comp);
          cnt++;
        }
      }
    }

    row r_prev = u.find_row(c_prev());

    int prev_deg = u.GetDoF().get_space_deg(*c_prev);
    int prev_time_deg = u.GetDoF().get_time_deg(*c_prev);
    int shift = (prev_time_deg - 1) * (r_prev.n() / prev_time_deg);

    const Scalar *u_r_prev = u(r_prev);
    const Scalar *u_r = u(r);

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

    for (int q = 0; q < elem.nQ(); ++q) {
      Tensor S = Zero;
      VectorField V = zero;

      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);

      for (int j = 0; j < elem.j_dimension() - r.n(); ++j) {
        Tensor phi_S_j = elem.Stress(q, j);
        VectorField phi_V_j = elem.Velocity(q, j);

        if (prev_deg != deg) {
          for (int k = 0; k < vc_prev.size(); ++k) {
            S += T_Matrix[j][k] * u_r_prev[k + shift] * phi_S_j;
            V += T_Matrix[j][k] * u_r_prev[k + shift] * phi_V_j;
          }
        } else {
          if (c.min() != 0.0) {
            S += u_r_prev[j + shift] * phi_S_j;
            V += u_r_prev[j + shift] * phi_V_j;
          } else {
            S += c_nodal_value[j] * phi_S_j;
            V += c_nodal_value[j] * phi_V_j;
          }
        }
      }
      for (int j = elem.j_dimension() - r.n(); j < elem.j_dimension(); ++j) {
        Tensor phi_S_j = elem.Stress(q, j);
        VectorField phi_V_j = elem.Velocity(q, j);
        S += u_r[j - elem.j_dimension() + r.n()] * phi_S_j;
        V += u_r[j - elem.j_dimension() + r.n()] * phi_V_j;
      }

      if (c.dim() == 2) {
        S -= prob->ut_Stress(QP);
        V -= prob->ut_Velocity(QP);
      } else Exit("ToDo")

      nrm += w * Frobenius(S, S);
      nrm += w * V * V;
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double TDGElasticityAssemble::GNError(Vector &) const {
  return -1.0;
}

double TDGElasticityAssemble::L2Norm(Vector &) const {Exit("not implemented"); }

void TDGElasticityAssemble::DualRHS_linear(Vector &RHS,
                                           const Vector &U,
                                           Point a,
                                           Point b) const {
  RHS = 0.0;

  if (a.t() != b.t()) { // volume evaluation
    for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {
      face f_diag = RHS.find_face(c.Face(c.Faces() - 2));
      cell c_prev = RHS.find_cell(f_diag.Left());
      if (c_prev == RHS.cells_end())
        c_prev = RHS.find_overlap_cell(f_diag.Left());
      if (c_prev == RHS.overlap_end())
        c_prev = c;

      SpaceTimeElasticityElement elem(*disc, RHS, c);

      row r_RHS = RHS.find_row(c());
      row r_RHS_prev = RHS.find_row(c_prev());

      for (int q = 0; q < elem.nQ(); ++q) {
        double w = elem.QWeight(q);
        Point QP = elem.QPoint(q);

        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;
        if (a.t() > QP.t() || QP.t() > b.t()) continue;

        Point QP_x = QP.CopyWithT(0.0);

        Tensor j_S = prob->ut_dual_Velocity(QP, a, b);
        VectorField j_V = prob->ut_dual_Velocity(QP, a, b);

        for (int j = 0; j < elem.j_dimension(); ++j) {
          Tensor S_j = elem.Stress(q, j);
          VectorField V_j = elem.Velocity(q, j);
          if (j < elem.j_dimension() - r_RHS.n()) {
            if (c_prev != c) {
              const int idx = 2 * r_RHS.n() - elem.j_dimension() + j;
              RHS(r_RHS_prev)[idx] -= w * (Frobenius(j_S, S_j) + j_V * V_j);
            }
          } else {
            const int idx = j - elem.j_dimension() + r_RHS.n();
            RHS(r_RHS)[idx] -= w * (Frobenius(j_S, S_j) + j_V * V_j);
          }
        }
      }
    }
  } else { // evaluation at some time point
    for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {
      if (c.max() != a.t()) continue;

      SpaceTimeElasticityElement elem(*disc, RHS, c);

      row r_RHS = RHS.find_row(c());
      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        double w = elem.SpaceQWeight(q);
        Point QP = elem.QPoint(q).CopyWithT(0.0);

        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;

        Tensor j_S = prob->ut_dual_Stress(QP, a, b);
        VectorField j_V = prob->ut_dual_Velocity(QP, a, b);

        int time_deg = U.GetDoF().get_time_deg(*c);
        double scale = 0.0;
        const int startj = time_deg * (elem.j_dimension() / (time_deg + 1));
        for (int j = startj; j < elem.j_dimension(); ++j) {
          scale += elem.Velocity(q, j)[0];
        }

        for (int j = startj; j < elem.j_dimension(); ++j) {
          Tensor S_j = elem.Stress(q, j);
          VectorField V_j = elem.Velocity(q, j);
          const int idx = j - elem.j_dimension() + r_RHS.n();
          RHS(r_RHS)[idx] -= w * (Frobenius(j_S, S_j) + j_V * V_j) / scale;
        }
      }
    }
  }
}

void TDGElasticityAssemble::DualRHS_quadratic(Vector &RHS,
                                              const Vector &U,
                                              Point a,
                                              Point b) const { //TODO
  RHS = 0;

  if (a.t() != b.t()) { // volume evaluation
    for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {

      face f_diag = RHS.find_face(c.Face(c.Faces() - 2));
      cell c_prev = RHS.find_cell(f_diag.Left());
      if (c_prev == RHS.cells_end())
        c_prev = RHS.find_overlap_cell(f_diag.Left());
      if (c_prev == RHS.overlap_end())
        c_prev = c;

      SpaceTimeElasticityElement elem(*disc, RHS, c);

      row r_RHS = RHS.find_row(c());
      row r_RHS_prev = RHS.find_row(c_prev());

      for (int q = 0; q < elem.nQ(); ++q) {
        double w = elem.QWeight(q);
        Point QP = elem.QPoint(q);

        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;
        if (a.t() > QP.t() || QP.t() > b.t()) continue;

        for (int i = 0; i < elem.j_dimension(); ++i) {
          Tensor S_i = elem.Stress(q, i);
          VectorField V_i = elem.Velocity(q, i);
          for (int j = 0; j < elem.j_dimension(); ++j) {
            Tensor S_j = elem.Stress(q, j);
            VectorField V_j = elem.Velocity(q, j);

            if (j < elem.j_dimension() - r_RHS.n()) {
              if (c_prev != c) {
                RHS(r_RHS_prev)[2 * r_RHS.n() - elem.j_dimension() + j] -=
                    w * U(r_RHS_prev)[i + 2 * r_RHS.n()
                                      - elem.j_dimension()]
                    * (Frobenius(S_i, S_j) + V_i * V_j);
              }
            } else {
              RHS(r_RHS)[j - elem.j_dimension() + r_RHS.n()] -=
                  w * U(r_RHS)[i - elem.j_dimension() + r_RHS.n()]
                  * (Frobenius(S_i, S_j) + V_i * V_j);
            }
          }
        }
      }
    }
  } else { // evaluation at some time point
    for (cell c = RHS.cells(); c != RHS.cells_end(); ++c) {
      if (c.max() != a.t()) continue;

      SpaceTimeElasticityElement elem(*disc, RHS, c);

      row r_RHS = RHS.find_row(c());

      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        double w = elem.SpaceQWeight(q);
        Point QP = elem.QPoint(q).CopyWithT(0.0);

        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;

        int time_deg = U.GetDoF().get_time_deg(*c);

        Scalar scale = 0.0;
        for (int j = time_deg * (elem.j_dimension() / (time_deg + 1));
             j < elem.j_dimension(); ++j)
          scale += elem.Velocity(q, j)[0];

        for (int i = time_deg * (elem.j_dimension() / (time_deg + 1));
             i < elem.j_dimension(); ++i) {
          Tensor S_i = (1. / scale) * elem.Stress(q, i);
          VectorField V_i = (1. / scale) * elem.Velocity(q, i);

          for (int j = time_deg * (elem.j_dimension() / (time_deg + 1));
               j < elem.j_dimension(); ++j) {
            Tensor S_j = (1. / scale) * elem.Stress(q, j);
            VectorField V_j = (1. / scale) * elem.Velocity(q, j);

            RHS(r_RHS)[j - elem.j_dimension() + r_RHS.n()] -=
                w * U(r_RHS)[i - elem.j_dimension() + r_RHS.n()]
                * (Frobenius(S_i, S_j) + V_i * V_j);
          }
        }
      }
    }
  }
}

double TDGElasticityAssemble::DualErrorEstimateCell(const cell &c,
                                                    const Vector &u,
                                                    const Vector &u_star) const {
  double eta = 0.0;
  //Date Start_cell;

  SpaceTimeElasticityElement elem(*disc, u, c);
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
  if (dist(c[0],c[2])>h)
      h = dist(c[0],c[2]);
  if (dist(c[1],c[2])>h)
      h = dist(c[1],c[2]);*/
  double h = max(double(dist(c[0], c[1])), double(dist(c[0], c[2])));
  h = max(h, dist(c[1], c[2]));
  h = sqrt(h * h + (c.max() - c.min()) * (c.max() - c.min()));

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

  vector<Scalar> vc(get_m(deg));
  vector<Scalar> vc_prev(get_m(prev_deg));
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

  double rho = prob->Rho(c());
  double mu = prob->Mu(c());
  double lambda = prob->Lambda(c());
  double s = 1 / (2.0 * mu);
  double ss = lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));

  for (int q = 0; q < elem.nQ(); ++q) {
    Tensor S = 0.0;
    VectorField V = 0.0;

    double w = elem.QWeight(q);
    Point QP = elem.QPoint(q);

    for (int j = 0; j < elem.j_dimension() - r.n(); ++j) {
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
          S += T_Matrix[j][k] * u_r_prev[k + shift] * (s * DtS_j - ss * TraceDtS_j - E_j);
          V += u_r_prev[k + shift] * (rho * DtV_j - DivS_j) * T_Matrix[j][k];
        }
      } else {
        if (c.min() != 0.0) {
          S += u_r_prev[j + shift] * (s * DtS_j - ss * TraceDtS_j - E_j);
          V += u_r_prev[j + shift] * (rho * DtV_j - DivS_j);
        } else {
          S += c_nodal_value[j] * (s * DtS_j - ss * TraceDtS_j - E_j);
          V += c_nodal_value[j] * (rho * DtV_j - DivS_j);
        }
      }
    }
    for (int j = elem.j_dimension() - r.n(); j < elem.j_dimension(); ++j) {
      Tensor DtS_j = elem.DtStress(q, j);
      Tensor TraceDtS_j = 0;
      TraceDtS_j[0][0] = trace(DtS_j);
      TraceDtS_j[1][1] = trace(DtS_j);
      if (c.dim() == 3)
        TraceDtS_j[2][2] = trace(DtS_j);
      VectorField DivS_j = elem.DivStress(q, j);
      VectorField DtV_j = elem.DtVelocity(q, j);
      Tensor E_j = elem.Strain(q, j);
      S += u_r[j - elem.j_dimension() + r.n()] * (s * DtS_j - ss * TraceDtS_j - E_j);
      V += u_r[j - elem.j_dimension() + r.n()] * (rho * DtV_j - DivS_j);
    }
    rho_1_S += w * Frobenius(S, S);
    rho_1_V += w * V * V;

    Tensor dual_S = 0.0;
    VectorField dual_V = 0.0;
    for (int i = 0; i < elem.i_dimension(); ++i) {
      dual_S += u_star_r[i] * elem.StressTestSpace(q, i);
      dual_V += u_star_r[i] * elem.VelocityTestSpace(q, i);
    }
    omega_1_S += w * Frobenius(dual_S, dual_S);
    omega_1_V += w * dual_V * dual_V;
  }

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
    cf_prev = u.find_cell_or_overlap_cell(ff.Left());
    row rf_prev = u.find_row(cf_prev());

    const Scalar *u_rf = u(rf);
    const Scalar *u_rf_prev = u(rf_prev);

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

    SpaceTimeElasticityFaceElement felem(*disc, u, c, f);
    SpaceTimeElasticityFaceElement felem_1(*disc, u, cf, f1);

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

    vector<Scalar> vcf(get_m(cf_deg));
    vector<Scalar> vcf_prev(get_m(cf_prev_deg));
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

    for (int q = 0; q < felem.nQ(); ++q) {
      int q1 = find_q_id(felem.QPoint(q), felem_1);
      Tensor S_Flux = 0.0;
      VectorField V_Flux = 0.0;

      double w = felem.QWeight(q);
      VectorField Nq = felem.QNormal(q);
      VectorField Tq = felem.QTangent(q);

      double c_s = prob->C_s(c());
      double c_p = prob->C_p(c());
      double c_s_f = prob->C_s(cf());
      double c_p_f = prob->C_p(cf());
      double mu_f = prob->Mu(cf());
      double lambda_f = prob->Lambda(cf());
      double alpha1 = ((2.0 * mu_f + lambda_f) * c_p)
                      / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
      double alpha2 = ((2.0 * mu_f + lambda_f) * (2.0 * mu + lambda))
                      / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
      double alpha3 = (mu_f * c_s) / (c_s * mu_f + c_s_f * mu);
      double alpha4 = (mu_f * mu) / (c_s * mu_f + c_s_f * mu);
      double alpha5 = (c_p_f * c_p)
                      / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
      double alpha6 = (c_p_f * (2.0 * mu + lambda))
                      / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
      double alpha7 = (c_s_f * c_s) / (c_s * mu_f + c_s_f * mu);
      double alpha8 = (c_s_f * mu) / (c_s * mu_f + c_s_f * mu);

      for (int j = 0; j < felem.j_dimension(); ++j) {
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
            + double(alpha3 * VFlux_c(V_j, Tq) + alpha7 * SFlux_c(S_j, Tq, Nq))
              * Tensor(Product(Tq, Nq));
        V_Flux -=
            double(alpha2 * VFlux_c(V_j, Nq) + alpha6 * SFlux_c(S_j, Nq, Nq)) * Nq
            + double(alpha4 * VFlux_c(V_j, Tq) + alpha8 * SFlux_c(S_j, Tq, Nq)) * Tq;
      }
      for (int j = 0; j < felem_1.j_dimension(); ++j) {
        Tensor S1_j = 0.0;
        VectorField V1_j = 0.0;
        if (j < felem_1.j_dimension() - rf.n()) {
          if (cf_deg != cf_prev_deg) {
            for (int k = 0; k < vcf_prev.size(); ++k) {
              S1_j += T_Matrix_f_prev[j % vcf.size()][k]
                      * u_rf_prev[k + shift_1] * felem_1.Stress(q1, j);
              V1_j += u_rf_prev[k + shift_1] * felem_1.Velocity(q1, j)
                      * T_Matrix_f_prev[j % vcf.size()][k];
            }
          } else {
            S1_j += u_rf_prev[j + shift_1] * felem_1.Stress(q1, j);
            V1_j += u_rf_prev[j + shift_1] * felem_1.Velocity(q1, j);
          }
        } else {
          S1_j += u_rf[j - felem_1.j_dimension() + rf.n()] * felem_1.Stress(q1, j);
          V1_j += u_rf[j - felem_1.j_dimension() + rf.n()] * felem_1.Velocity(q1, j);
        }

        S_Flux -= double(alpha1 * VFlux_cf(V1_j, Nq, bnd) + alpha5 * SFlux_cf(S1_j, Nq, Nq, bnd))
                  * Tensor(Product(Nq, Nq))
                  + double(alpha3 * VFlux_cf(V1_j, Tq, bnd)
                           + alpha7 * SFlux_cf(S1_j, Tq, Nq, bnd))
                    * Tensor(Product(Tq, Nq));
        V_Flux -= double(alpha2 * VFlux_cf(V1_j, Nq, bnd) + alpha6 * SFlux_cf(S1_j, Nq, Nq, bnd)) * Nq
                  + double(alpha4 * VFlux_cf(V1_j, Tq, bnd)
                           + alpha8 * SFlux_cf(S1_j, Tq, Nq, bnd)) * Tq;
      }
      rho_2_S += w * 1. / h * Frobenius(S_Flux, S_Flux);
      rho_2_V += w * 1. / h * V_Flux * V_Flux;

      Tensor dual_S_Flux = 0.0;
      VectorField dual_V_Flux = 0.0;
      for (int i = 0; i < felem.i_dimension(); ++i) {
        dual_S_Flux += u_star_r[i] * felem.StressTestSpace(q, i);
        dual_V_Flux += u_star_r[i] * felem.VelocityTestSpace(q, i);
      }
      omega_2_S += w * h * Frobenius(dual_S_Flux, dual_S_Flux);
      omega_2_V += w * h * dual_V_Flux * dual_V_Flux;
    }
  }
  //est_face_time += Date() - Start_face;

  eta = sqrt(rho_1_V + rho_1_S) * sqrt(omega_1_V + omega_1_S)
        + sqrt(rho_2_V + rho_2_S) * sqrt(omega_2_V + omega_2_S);
  return eta;
}

double TDGElasticityAssemble::DualErrorEstimateJumpCell
    (const cell &c, const Vector &u, const Vector &u_star) const {
  double eta = 0.0;
  //Date Start_cell;

  SpaceTimeElasticityElement elem(*disc, u, c);
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
    face f = u.find_face(c.Face(c.Faces() - 2));
    c_prev = u.find_cell(f.Left());
    if (c_prev == u.cells_end())
      c_prev = u.find_overlap_cell(f.Left());
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

  vector<Scalar> vc(get_m(deg));
  vector<Scalar> vc_prev(get_m(prev_deg));
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

  double rho = prob->Rho(c());
  double mu = prob->Mu(c());
  double lambda = prob->Lambda(c());
  double s = 1 / (2.0 * mu);
  double ss = lambda / (2.0 * mu * (2.0 * lambda + 2.0 * mu));

  for (int q = 0; q < elem.nQ(); ++q) {
    Tensor S = 0.0;
    VectorField V = 0.0;

    double w = elem.QWeight(q);
    Point QP = elem.QPoint(q);

    for (int j = 0; j < elem.j_dimension() - r.n(); ++j) {
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
          S += T_Matrix[j][k] * u_r_prev[k + shift] * (s * DtS_j - ss * TraceDtS_j - E_j);
          V += u_r_prev[k + shift] * (rho * DtV_j - DivS_j) * T_Matrix[j][k];
        }
      } else {
        if (c.min() != 0.0) {
          S += u_r_prev[j + shift] * (s * DtS_j - ss * TraceDtS_j - E_j);
          V += u_r_prev[j + shift] * (rho * DtV_j - DivS_j);
        } else {
          S += c_nodal_value[j] * (s * DtS_j - ss * TraceDtS_j - E_j);
          V += c_nodal_value[j] * (rho * DtV_j - DivS_j);
        }
      }
    }
    for (int j = elem.j_dimension() - r.n(); j < elem.j_dimension(); ++j) {
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

    SpaceTimeElasticityFaceElement felem(*disc, u, c, f);
    SpaceTimeElasticityFaceElement felem_1(*disc, u, cf, f1);
    SpaceTimeElasticityElement elem_1(*disc, u, cf);

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

    vector<Scalar> vcf(get_m(cf_deg));
    vector<Scalar> vcf_prev(get_m(cf_prev_deg));
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

    double fxI = 0;

    for (int q = 0; q < felem.nQ(); ++q) {
      int q1 = find_q_id(felem.QPoint(q), felem_1);
      VectorField V_Flux = 0.0;
      Tensor S_Flux = 0.0;

      double w = felem.QWeight(q);
      VectorField Nq = felem.QNormal(q);
      VectorField Tq = felem.QTangent(q);

      fxI += w;

      double c_s = prob->C_s(c());
      double c_p = prob->C_p(c());
      double c_s_f = prob->C_s(cf());
      double c_p_f = prob->C_p(cf());
      double mu_f = prob->Mu(cf());
      double lambda_f = prob->Lambda(cf());
      double alpha1 = ((2.0 * mu_f + lambda_f) * c_p)
          / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
      double alpha2 = ((2.0 * mu_f + lambda_f) * (2.0 * mu + lambda))
                      / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
      double alpha3 = (mu_f * c_s) / (c_s * mu_f + c_s_f * mu);
      double alpha4 = (mu_f * mu) / (c_s * mu_f + c_s_f * mu);
      double alpha5 = (c_p_f * c_p)
                      / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
      double alpha6 = (c_p_f * (2.0 * mu + lambda))
                   / (c_p * (2.0 * mu_f + lambda_f) + c_p_f * (2.0 * mu + lambda));
      double alpha7 = (c_s_f * c_s) / (c_s * mu_f + c_s_f * mu);
      double alpha8 = (c_s_f * mu) / (c_s * mu_f + c_s_f * mu);

      for (int j = 0; j < felem.j_dimension(); ++j) {
        Tensor S_j = 0.0;
        VectorField V_j = 0.0;
        if (j < felem.j_dimension() - r.n()) {
          if (deg != prev_deg) {
            for (int k = 0; k < vc_prev.size(); ++k) {
              S_j += T_Matrix[j][k] * u_r_prev[k + shift] * felem.Stress(q, j);
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
        S_Flux +=double(alpha1 * VFlux_c(V_j, Nq) + alpha5 * SFlux_c(S_j, Nq, Nq))
            * Tensor(Product(Nq, Nq))
            + double(
                alpha3 * VFlux_c(V_j, Tq) + alpha7 * SFlux_c(S_j, Tq, Nq))
              * Tensor(Product(Tq, Nq));
        V_Flux -= double(alpha2 * VFlux_c(V_j, Nq) + alpha6 * SFlux_c(S_j, Nq, Nq)) * Nq
            + double(
                alpha4 * VFlux_c(V_j, Tq) + alpha8 * SFlux_c(S_j, Tq, Nq))
              * Tq;
      }
      for (int j = 0; j < felem_1.j_dimension(); ++j) {
        Tensor S1_j = 0.0;
        VectorField V1_j = 0.0;
        if (j < felem_1.j_dimension() - rf.n()) {
          if (cf_deg != cf_prev_deg) {
            for (int k = 0; k < vcf_prev.size(); ++k) {
              S1_j += T_Matrix_f_prev[j % vcf.size()][k]
                      * u_rf_prev[k + shift_1] * felem_1.Stress(q1, j);
              V1_j += u_rf_prev[k + shift_1] * felem_1.Velocity(q1, j)
                  * T_Matrix_f_prev[j % vcf.size()][k];
            }
          } else {
            S1_j += u_rf_prev[j + shift_1] * felem_1.Stress(q1, j);
            V1_j += u_rf_prev[j + shift_1] * felem_1.Velocity(q1, j);
          }
        } else {
          S1_j += u_rf[j - felem_1.j_dimension() + rf.n()] * felem_1.Stress(q1, j);
          V1_j += u_rf[j - felem_1.j_dimension() + rf.n()] * felem_1.Velocity(q1, j);
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

    vector<VectorField> dual_V_projection_cf(time_deg_cf, 0.0);
    vector<Tensor> dual_S_projection_cf(time_deg_cf, 0.0);

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

    const Shape &time_shape_cf(disc->GetTimeShape(time_deg_cf - 1));

    double f_vol = fxI / I;


    auto &timeQuad = disc->GetTimeQuad(time_deg_cf - 1);

    for (int q = 0; q < timeQuad.size(); ++q) {
      VectorField values_dual_V = 0.0;
      VectorField values_dual_S = 0.0;
      for (int l = 0; l < time_shape.size(); ++l) {
        values_dual_V += dual_V_projection[l] * time_shape.values()[q][l];
        values_dual_S +=
            dual_S_projection[l] * felem.QNormal(0) * time_shape.values()[q][l];
      }

      VectorField values_dual_V_cf = 0.0;
      VectorField values_dual_S_cf = 0.0;
      for (int l = 0; l < time_shape_cf.size(); ++l) {
        values_dual_V_cf += dual_V_projection_cf[l] * time_shape_cf.values()[q][l];
        values_dual_S_cf += dual_S_projection_cf[l] * felem.QNormal(0)
                            * time_shape_cf.values()[q][l];
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

double TDGElasticityAssemble::ErrorEstimateCell(const cell &c,
                                                const Vector &u) const { // TODO
  Exit("Not implemented");
}

double TDGElasticityAssemble::Goal_Functional(const Vector &u, Point a, Point b) const {
  Scalar e = 0.0;

  if (a.t() != b.t()) { // volume evaluation
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      SpaceTimeElasticityElement elem(*disc, u, c);
      row r = u.find_row(c());

      cell c_prev = c;
      vector<Point> c_nodal;

      if (c.min() != 0) {
        face f = u.find_face(c.Face(c.Faces() - 2));
        c_prev = u.find_cell_or_overlap_cell(f.Left());
      } else {
        c_nodal = u.GetDoF().GetNodalPoints(*c);
        for (int i = 0; i < c_nodal.size(); ++i) c_nodal[i] = c_nodal[i].CopyWithT(0.0);
      }

      row r_prev = u.find_row(c_prev());

      int prev_time_deg = u.GetDoF().get_time_deg(*c_prev);

      vector<Scalar> vc(u.GetDoF().get_m(c()));
      vector<Scalar> vc_prev(u.GetDoF().get_m(c_prev()));

      vector<vector<Scalar>> T_Matrix(vc.size());
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
        Point QP = elem.QPoint(q);
        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;
        if (a.t() > QP.t() || QP.t() > b.t()) continue;
        double w = elem.QWeight(q);

        Tensor j_S = prob->ut_dual_Stress(QP, a, b);
        VectorField j_V = prob->ut_dual_Velocity(QP, a, b);

        Scalar U = 0.0;
        for (int j = 0; j < elem.j_dimension(); ++j) {
          Tensor S_j = elem.Stress(q, j);
          VectorField V_j = elem.Velocity(q, j);

          if (j < elem.j_dimension() - r.n()) {
            if (c.min() == 0) {
              //for (int i = 0; i < prob->SpaceDim(); ++i)
                //U += prob->ut(c_nodal[j], COMPONENT::) * sqrt(Frobenius(S_j, S_j) + V_j * V_j);
              Exit("check");
            } else
              for (int k = 0; k < vc_prev.size(); ++k) {
                const int idx = j + (prev_time_deg - 1) * (r_prev.n() / prev_time_deg);
                U += u(r_prev, idx) * (Frobenius(j_S, S_j) + j_V * V_j) * T_Matrix[j][k];
              }
          } else {
            U += u(r, j - elem.j_dimension() + r.n()) * (Frobenius(j_S, S_j) + j_V * V_j);
          }
        }
        e += w * U;
      }
    }
  } else { // evaluation at some time point
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != a.t()) continue;

      SpaceTimeElasticityElement elem(*disc, u, c);

      row r = u.find_row(c());

      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        Point QP = elem.QPoint(q).CopyWithT(0.0);
        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;
        double w = elem.SpaceQWeight(q);

        Tensor j_S = prob->ut_dual_Stress(QP, a, b);
        VectorField j_V = prob->ut_dual_Velocity(QP, a, b);

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
  return abs(PPM->SumOnCommSplit(e, u.CommSplit()));
}

double TDGElasticityAssemble::Energy(const Vector &u,
                                     Point a,
                                     Point b) const {//Exit("CHECK");
  Scalar e = 0;

  if (a.t() != b.t()) { // volume evaluation
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      SpaceTimeElasticityElement elem(*disc, u, c);
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
        for (int i = 0; i < c_nodal.size(); ++i) c_nodal[i] = c_nodal[i].CopyWithT(0);
      }

      row r_prev = u.find_row(c_prev());

      int time_deg = u.GetDoF().get_time_deg(*c_prev);

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
        double w = elem.QWeight(q);
        Point QP = elem.QPoint(q);
        VectorField V = 0.0;
        Tensor S = 0.0;
        if ((a[0] <= QP[0]) && (QP[0] <= b[0]))
          if ((a[1] <= QP[1]) && (QP[1] <= b[1]))
            if (a.t() <= QP.t() && QP.t() <= b.t())
              for (int j = 0; j < elem.j_dimension(); ++j) {
                VectorField V_j = elem.Velocity(q, j);
                Tensor S_j = elem.Stress(q, j);

                if (j < elem.j_dimension() - r.n()) {
                  if (c.min() == 0) {
                    Exit("check!")
                    /*S[0] += prob->ut(c_nodal[j], 0) * S_j[0];
                    S[1] += prob->ut(c_nodal[j], 1) * S_j[1];
                    S[2] += prob->ut(c_nodal[j], 2) * S_j[2];
                    V[0] += prob->ut(c_nodal[j], 3) * V_j[0];
                    V[1] += prob->ut(c_nodal[j], 4) * V_j[1];
                     */
                  } else
                    for (int k = 0; k < vc_prev.size(); ++k) {
                      const int idx = j + (time_deg - 1) * (r_prev.n() / time_deg);
                      V += (u(r_prev,idx) * T_Matrix[j][k]) * V_j;
                      S += (u(r_prev, idx) * T_Matrix[j][k]) * S_j;
                    }
                } else {
                  V += u(r, j - elem.j_dimension() + r.n()) * V_j;
                  S += u(r, j - elem.j_dimension() + r.n()) * S_j;
                }
              }
        e += 0.5 * w * (V * V + Frobenius(S, S));
      }
    }
  } else { // evaluation at some time point
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != a.t()) continue;

      SpaceTimeElasticityElement elem(*disc, u, c);

      row r = u.find_row(c());

      for (int q = 0; q < elem.nSpaceQ(); ++q) {
        double w = elem.SpaceQWeight(q);
        Point QP = elem.QPoint(q).CopyWithT(0.0);

        if ((a[0] > QP[0]) || (QP[0] > b[0])) continue;
        if ((a[1] > QP[1]) || (QP[1] > b[1])) continue;

        int time_deg = u.GetDoF().get_time_deg(*c);

        Scalar scale = 0.0;
        for (int j = time_deg * (elem.j_dimension() / (time_deg + 1));
             j < elem.j_dimension(); ++j)
          scale += elem.Velocity(q, j)[0];

        VectorField V = 0.0;
        Tensor S = 0.0;
        const int startj = time_deg * (elem.j_dimension() / (time_deg + 1));
        for (int j = startj; j < elem.j_dimension(); ++j) {
          VectorField V_j = elem.Velocity(q, j);
          Tensor S_j = elem.Stress(q, j);
          V += (u(r, j - elem.j_dimension() + r.n()) / scale) * V_j;
          S += (u(r, j - elem.j_dimension() + r.n()) / scale) * S_j;
        }
        e += 0.5 * w * (V * V + Frobenius(S, S));
      }
    }
  }
  return abs(PPM->SumOnCommSplit(e, u.CommSplit()));
}