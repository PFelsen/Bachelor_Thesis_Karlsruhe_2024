#include "STGLViscoAcousticAssemble.hpp"
#include "STDGDGViscoAcousticElementGLGL.hpp"
#include "STDGDGViscoAcousticFaceElementGLGL.hpp"

std::unique_ptr<Operator> STGLViscoAcousticAssemble::GetMatrixFreeOperator() const {
  return std::make_unique<STGLViscoAcousticSystemOperator>(*this);
}

void STGLViscoAcousticSystemOperator::multiply_plus(Vector &b, const Vector &x) const {
  mout.StartBlock("MFapply");
  int numL = problem.nL();
  for (cell c = b.cells(); c != b.cells_end(); ++c) {
    STDGDGViscoAcousticElementGLGL elem(b, c, numL);
    row r_c = b.find_row(c());
    const double rho = problem.Rho(*c,  c());

    vector<double> kappaInv(2 + numL);
    vector<double> kappaTauInv(2 + numL);
    kappaInv[1] = 1.0 / problem.Kappa_i(*c, c(), 0);
    for (int j = 2; j < 2 + numL; j++) {
      kappaInv[j] = 1.0 / problem.Kappa_i(*c, c(), j - 1);
      kappaTauInv[j] = kappaInv[j] / problem.Tau_i(*c, c(), j - 1);
    }


    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      for (int i = 0; i < elem.i_dimension(); ++i) {
        int var_i = elem.variable(i);
        VectorField V_i = elem.Velocity(q, i);
        double P_i = elem.Pressure(q, i);
        for (int j = 0; j < elem.j_dimension(); ++j) {
          int var_j = elem.variable(j);
          VectorField DtV_j = elem.DtVelocity(q, j);
          Scalar DivV_j = elem.DivVelocity(q, j);
          Scalar P_j = elem.Pressure(q, j);
          Scalar DtP_j = elem.DtPressure(q, j);
          VectorField GradP_j = elem.GradPressure(q, j);
          if (var_i != var_j) DtP_j = 0.0;
          b(r_c, i) += x(r_c, j) * w * ((rho * (DtV_j * V_i)
                                         + kappaInv[var_j] * (DtP_j * P_i))
                                        - ((GradP_j * V_i) + (DivV_j * P_i)));
          if (var_j > 1 && var_i == var_j) {
            b(r_c, i) += x(r_c, j) * w * kappaTauInv[var_j] * (P_j * P_i);
          }
        }
      }
    }
  }

  for (cell c = b.cells(); c != b.cells_end(); ++c) {
    const double rho = problem.Rho(*c, c());
    vector<double> kappaInv(2 + numL);
    vector<double> kappaTauInv(2 + numL);
    kappaInv[1] = 1.0 / problem.Kappa(*c, c());
    for (int j = 2; j < 2 + numL; j++) {
      kappaInv[j] = 1.0 / problem.Kappa_i(*c, c(), j - 1);
      kappaTauInv[j] = kappaInv[j] / problem.Tau_i(*c, c(), j - 1);
    }
    row r_c = b.find_row(c());
    const double z_K = sqrt(problem.Rho(*c, c()) * problem.Kappa(*c, c()));
    for (int f = 0; f < c.Faces() - 2; ++f) {
      bool bnd = b.OnBoundary(*c, f);
      int bnd_id = bnd ? problem.BndID(c.Face(f)) : -1;
      cell cf = bnd ? c : b.find_neighbour_cell(c, f);
      int f1 = b.find_neighbour_face_id(c.Face(f), cf);
      STDGDGViscoAcousticFaceElementGLGL felem(assemble.GetDisc(), b, c, f, numL);
      STDGDGViscoAcousticFaceElementGLGL felem_1(assemble.GetDisc(), b, cf, f1, numL);
      row r_cf = b.find_row(cf());

      const double z_Kf = sqrt(problem.Rho(*c, cf()) * problem.Kappa(*c, cf()));
      const double alpha1 = 1.0 / (z_K + z_Kf);
      const double alpha2 = z_Kf * z_K * alpha1;
      const double alpha3 = z_Kf * alpha1;
      const double alpha4 = z_K * alpha1;

      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        double w = felem.QWeight(q);
        VectorField Nq = felem.QNormal(q);

        for (int i = 0; i < felem.i_dimension(); ++i) {
          if (felem.is_zero_test(q, i)) continue;
          VectorField V_i = felem.Velocity(q, i);
          double P_i = felem.Pressure(q, i);

          for (int j = 0; j < felem.j_dimension(); ++j) {
            if (felem.is_zero(q, j)) continue;
            VectorField V_j = felem.Velocity(q, j);
            double P_j = felem.Pressure(q, j);

            double val = (alpha1 * P_j * P_i
                          + alpha2 * VFlux_c(V_j, Nq) * VFlux_c(V_i, Nq)
                          + alpha3 * VFlux_c(V_j, Nq) * P_i
                          + alpha4 * P_j * VFlux_c(V_i, Nq));

            b(r_c, i) += x(r_c, j) * w * val;
          }
          for (int j = 0; j < felem_1.j_dimension(); ++j) {
            if (felem_1.is_zero(q1, j)) continue;
            VectorField V1_j = felem_1.Velocity(q1, j);
            double P1_j = felem_1.Pressure(q1, j);

            b(r_c, i) -= x(r_cf, j) * w * (alpha1 * PFlux_cf(P1_j, bnd, bnd_id) * P_i
                                           + alpha2 * VFlux_cf(V1_j, Nq, bnd, bnd_id) *
                                             VFlux_c(V_i, Nq)
                                           + alpha3 * VFlux_cf(V1_j, Nq, bnd, bnd_id) * P_i
                                           + alpha4 * PFlux_cf(P1_j, bnd, bnd_id) *
                                             VFlux_c(V_i, Nq));
          }
        }
      }
    }
    cell c_prev = b.find_previous_cell(c);
    row r_cprev = b.find_row(c_prev());
    STDGDGViscoAcousticFaceElementGLGL felem(assemble.GetDisc(), b, c, c.Faces() - 2, numL, "dg");
    STDGDGViscoAcousticFaceElementGLGL felem_1(assemble.GetDisc(), b, c_prev, c_prev.Faces() - 1, numL,
                                                     "dg", c);
    for (int q = 0; q < felem.nQ(); ++q) {
      double w = felem.QWeight(q);
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
          b(r_c, i) += x(r_c, i) * w * (rho * (V_j * V_i) + kappaInv[var_i] * P_j * P_i);
        }

        if (c != c_prev && c.min() > 0) {
          for (int j = 0; j < felem_1.j_dimension(); ++j) {
            int var_j = felem_1.variable(j);
            if (felem_1.is_zero(q, j) || var_i != var_j) continue;
            VectorField V1_j = felem_1.Velocity(q, j);
            double P1_j = felem_1.Pressure(q, j);
            b(r_c, i) -= x(r_cprev, j) * w * (rho * (V1_j * V_i) + kappaInv[var_j] * P1_j * P_i);
          }
        }
      }
    }
  }
  b.Collect();
  mout.EndBlock();
}


