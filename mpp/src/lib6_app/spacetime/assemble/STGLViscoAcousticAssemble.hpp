#ifndef STGLGLVISCOACOUSTIC_HPP
#define STGLGLVISCOACOUSTIC_HPP

#include "STDGViscoAcousticAssemble.hpp"
#include "SpaceTimeDiscretization_DGDG_GLGL.hpp"
#include "STDGDGViscoAcousticElementGLGL.hpp"
#include "STDGDGViscoAcousticFaceElementGLGL.hpp"

class STGLViscoAcousticSystemOperator;

class STGLViscoAcousticAssemble : public STDGViscoAcousticAssemble {
protected:
public:
  STGLViscoAcousticAssemble(const Meshes &meshes, DegreePair degree,
                            const string &problemName) :
      STDGViscoAcousticAssemble(meshes, degree, problemName) {
    disc = std::make_shared<STDiscretization_DGDG_GLGL>(meshes, degree, prob->Dim(), false);
    Config::Get<double>("tau_zero_inv", tau_zero_inv);
  }

  STGLViscoAcousticAssemble(const Meshes &meshes, DegreePair degree,
                            const std::shared_ptr<AcousticProblem> problem) :
      STDGViscoAcousticAssemble(meshes, degree, problem) {
    disc = std::make_shared<STDiscretization_DGDG_GLGL>(meshes, degree, prob->Dim(), false);
    Config::Get<double>("tau_zero_inv", tau_zero_inv);
  }

  const char* Name() const override { return "STGLViscoAcousticAssemble"; }

  bool HasMatrixFreeOperator() const override { return true; }

  std::unique_ptr<Operator> GetMatrixFreeOperator() const override;

  void System(Matrix &M, Vector &RHS) const override {

    Date Start;
    mout.PrintInfo("STGLViscoAcousticAssemble", verbose,
                   PrintInfoEntry<int>("Problem Size", M.pSize(), 2),
                   PrintInfoEntry<Date>("start assemble", Start, 2)
    );
    Time t_cell;
    Time t_face;
    M = 0;
    for (cell c = M.cells(); c != M.cells_end(); ++c) {
      Date Start_cell;
      STDGDGViscoAcousticElementGLGL elem(M, c, prob->nL());
      row r_c = M.find_row(c());
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
        Point QP = elem.QPoint(q);
        double t = QP.t();
        double w = elem.QWeight(q);
        VectorField V_rhs = prob->F_Velocity(t, *c, QP);
        Scalar P_rhs = prob->F_Pressure(t, *c, QP);

        for (int i = 0; i < elem.i_dimension(); ++i) {
          int var_i = elem.variable(i);
          Velocity V_i = elem.Velocity(q, i);
          double P_i = elem.Pressure(q, i);
          if (prob->HasRHS()) {
            RHS(c(), i) += w * (V_rhs * V_i + P_rhs * P_i);
          }
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

      t_cell += Date() - Start_cell;
      Date Start_face;

      const double z_K = sqrt(prob->Rho(*c, c()) * prob->Kappa(*c, c()));
      for (int f = 0; f < c.Faces() - 2; ++f) {
        bool bnd = M.OnBoundary(*c, f);
        int bnd_id = bnd ? prob->BndID(c.Face(f)) : -1;
        cell cf = bnd ? c : M.find_neighbour_cell(c, f);
        STDGDGViscoAcousticFaceElementGLGL felem(*disc, M, c, f, prob->nL());
        int f1 = M.find_neighbour_face_id(c.Face(f), cf);
        STDGDGViscoAcousticFaceElementGLGL felem_1(*disc, M, cf, f1, prob->nL());


        DGRowEntries M_cf(M, *c, *cf, false);

        const double z_Kf = sqrt(prob->Rho(*cf, cf()) * prob->Kappa(*cf, cf()));
        const double alpha1 = 1.0 / (z_K + z_Kf);
        const double alpha2 = z_Kf * z_K * alpha1;
        const double alpha3 = z_Kf * alpha1;
        const double alpha4 = z_K * alpha1;
        for (int q = 0; q < felem.nQ(); ++q) {
          int q1 = find_q_id(felem.QPoint(q), felem_1);
          Point QP = felem.QPoint(q);
          double w = felem.QWeight(q);
          VectorField Nq = felem.QNormal(q);

          Scalar VN_bnd = Nq * prob->ut_Velocity(QP.t(), QP, *c);
          Scalar P_bnd = prob->ut_Pressure(QP.t(), QP, *c);

          for (int i = 0; i < felem.i_dimension(); ++i) {
            if (felem.is_zero_test(q, i)) continue;
            VectorField V_i = felem.Velocity(q, i);
            double P_i = felem.Pressure(q, i);
            double VFlux_c_Vi = VFlux_c(V_i, Nq);

            if (bnd_id == 1) {
              RHS(c(), i) += 2 * w * (alpha2 * P_i + alpha3 * VFlux_c_Vi) * P_bnd;
            } else if (bnd_id == 2 || bnd_id == 0) {
              RHS(c(), i) += 2 * w * (alpha1 * VFlux_c_Vi + alpha4 * P_i) * VN_bnd;
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
      STDGDGViscoAcousticFaceElementGLGL felem(*disc, M, c, c.Faces() - 2, prob->nL(), "dg");
      STDGDGViscoAcousticFaceElementGLGL felem_1(*disc, M, c_prev, c_prev.Faces() - 1,
                                                 prob->nL(), "dg", c);

      DGRowEntries M_c_prev(M, *c, *c_prev, false);
      for (int q = 0; q < felem.nQ(); ++q) {
        Point QP = felem.QPoint(q);
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
            M_c(i, j) += w * (rho * (V_j * V_i) + kappaInv[var_j] * P_j * P_i);
          }
          if (c.min() == 0.0) {
            VectorField V1_j = prob->ut_Velocity(QP.t(), QP, *c);
            double P1_j = prob->ut_Pressure(QP.t(), QP, *c);
            RHS(c(), i) += w * (rho * (V1_j * V_i) + kappaInv[var_i] * P1_j * P_i);
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

    t_cell.Max();
    t_face.Max();

    mout.PrintInfo("STGLViscoAcousticAssemble", verbose,
                   PrintInfoEntry<Time>("cell assemble time", t_cell),
                   PrintInfoEntry<Time>("face assemble time", t_face),
                   PrintInfoEntry<Time>("assemble time", Date() - Start),
                   PrintInfoEntry<Date>("finish assemble", Date()));
  }

};

class STGLViscoAcousticSystemOperator : public Operator {
private:
  const STGLViscoAcousticAssemble &assemble;
  const STDiscretization &disc;
  const AcousticProblem &problem;
public:
  STGLViscoAcousticSystemOperator(const STGLViscoAcousticAssemble &assemble,
                                  bool assembleDiagonal = false)
      : assemble(assemble), disc(assemble.GetDisc()), problem(assemble.GetProblem()) {}

  void multiply_plus(Vector &b, const Vector &x) const override;
};

#endif //STGLGLVISCOACOUSTIC_HPP
