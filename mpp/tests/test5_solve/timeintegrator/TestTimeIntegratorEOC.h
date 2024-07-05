//
// Created by cr on 14.11.21.
//

#ifndef FWI_TESTTIMEINTEGRATOREOC_H
#define FWI_TESTTIMEINTEGRATOREOC_H

#include <Assemble.hpp>
#include <Plotting.hpp>

class TestAssembleEOC : public ILinearTimeAssemble {
public:
  TestAssembleEOC(int i) : ILinearTimeAssemble(0, 0.25, 0.25 * std::pow(2, -(1 + i))) {}

  vector<double> t_vec;
  vector<Vector> u_vec;

  void SetInitialValue(Vector &u) override {
    t_vec.push_back(Time());
    u.Clear();
    Solution(Time(), u);
  }

  const char *Name() const override { return "TestAssemble"; };

  void Initialize(Vector &u) const override {}

  // Identity Matrix
  void MassMatrix(Matrix &massMatrix) const override {
    massMatrix = 0.0;
    const Mesh &Me = massMatrix.GetMesh();
    //        MOUT( Me.MeshWidthStr());
    auto h = Me.MeshWidth();
    for (cell c = massMatrix.cells(); c != massMatrix.cells_end(); ++c) {
      const Cell &C = *c;
      DGRowEntries A_c(massMatrix, C, C);
      A_c(0, 0) = h.first;
      A_c(1, 1) = h.first;
    }
  }

  void SystemMatrix(Matrix &systemMatrix) const override {
    systemMatrix = 0.0;
    for (cell c = systemMatrix.cells(); c != systemMatrix.cells_end(); ++c) {
      const Cell &C = *c;

      DGRowEntries A_c(systemMatrix, C, C);
      for (int f = 0; f < c.Faces(); ++f) {
        if (systemMatrix.OnBoundary(C, f)) {
          if (f == 0) {
            A_c(0, 1) = 0.5;
            A_c(1, 0) = -0.5;
            A_c(1, 1) = -1;

          } else if (f == 1) {
            A_c(0, 1) = -0.5;
            A_c(1, 0) = 0.5;
            A_c(1, 1) = -1.5;
          } else {
            Exit("error!!");
          }
        } else {
          cell cf = systemMatrix.find_neighbour_cell(c, f);
          const Cell &CF = *cf;
          DGRowEntries A_cf(systemMatrix, C, CF);
          double direc = (C()[0] < CF()[0]) ? 1.0 : -1.0;
          A_c(0, 0) -= 0.5;
          A_cf(0, 0) += 0.5;
          A_cf(0, 1) += 0.5 * direc;
          A_c(1, 1) -= 0.5;
          A_cf(1, 1) += 0.5;
          A_cf(1, 0) += 0.5 * direc;
        }
      }
    }
    systemMatrix.CreateSparse();
  }

  void RHS(double t, Vector &rhs) const override {
    rhs.Clear();
    auto h = rhs.GetMesh().MeshWidth();
    for (cell c = rhs.cells(); c != rhs.cells_end(); ++c) {
      const Cell &C = *c;
      double x = C()[0];
      rhs(C(), 0) = h.first * -pDerivX(x, t) + h.first * vDerivT(x, t);
      rhs(C(), 1) = h.first * -vDerivX(x, t) + h.first * pDerivT(x, t);
      for (int f = 0; f < C.Faces(); ++f) {
        if (rhs.OnBoundary(C, f)) {
          if (f == 0) {
            double y = pSolution(0, t);
            rhs(C(), 0) -= y;
            rhs(C(), 1) += y;

          } else if (f == 1) {
            double y = pSolution(1, t);
            rhs(C(), 0) += y;
            rhs(C(), 1) += y;
          }
        }
      }
    }
  }

  void FinishTimeStep(const Vector &u) override {
    //        MOUT(u);
    //        mpp::plot("test") << u << mpp::endp;
  }

  double pDerivX(double x, double t) const { return cos(t); }

  //    double pDerivX (double x,double t) const {
  //        return -sin(x)*cos(t);
  //    }
  double pDerivT(double x, double t) const { return -sin(t) * x; }

  //    double pDerivT (double x,double t) const {
  //        return -sin(t)*cos(x);
  //    }
  double vDerivX(double x, double t) const { return sin(t); }

  //    double vDerivX (double x,double t) const {
  //        return cos(x)*sin(t);
  //    }
  double vDerivT(double x, double t) const { return cos(t) * x; }

  //    double vDerivT (double x,double t) const {
  //        return cos(t)*sin(x);
  //    }
  double pSolution(double x, double t) const {
    //        return cos(x)*cos(t);
    return x * cos(t);
  }

  double vSolution(double x, double t) const {
    //        return sin(x)*sin(t);
    return x * sin(t);
  }

  virtual void Solution(double t, Vector &uEx) {
    for (cell c = uEx.cells(); c != uEx.cells_end(); ++c) {
      const Cell &C = *c;
      double x = C()[0];
      uEx(C(), 0) = vSolution(x, t);
      uEx(C(), 1) = pSolution(x, t);
    }
  };
};

class TestTimeIntegratorEOC {};


#endif // FWI_TESTTIMEINTEGRATOREOC_H
