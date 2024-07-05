#ifndef _DGACUSTIC_H_
#define _DGACUSTIC_H_

#include "Assemble.hpp"
#include "SpaceTimeTools.hpp"
#include "CImg.hpp"
#include "LinearSolver.hpp"
#include "Itertools.hpp"
#include "ProblemBase.hpp"
#include "ObservationSpecification.hpp"
#include "SeismogramData.hpp"

#include <unordered_map>
#include <iostream>

class STAssemble : public ILinearAssemble {
protected:
  virtual STDiscretization &getDisc() = 0;

  int ErrorEstimatorVerbose = 0;
public:
  STAssemble() : ILinearAssemble() {
    Config::Get("ErrorEstimatorVerbose", ErrorEstimatorVerbose);
  }

  ~STAssemble() override = default;

  virtual void SetProblem(std::shared_ptr<IProblem> problem) { THROW("Implement SetProblem in STAssemble"); }

  virtual const IProblem& GetAbstractProblem() const = 0;

  virtual const STDiscretization &GetDisc() const = 0;

  virtual std::shared_ptr<const STDiscretization> GetSharedDisc() const = 0;

  void PrintInfo() const override {
    mout.PrintInfo("Assemble", verbose,
                   PrintInfoEntry("Name", Name(), 1),
                   PrintInfoEntry("Discretization Name", GetDisc().DiscName(), 1));
  }

  virtual SeismogramData measure(const ObservationSpecification &specification, const Vector &u) const {
    throw NotImplementedException();
  }

  virtual void PlotSingleSolution(const Vector &u, const std::string &filename) const = 0;

  virtual void PlotParameters(const Vector &u, std::string suffix) const{};

  const char *Name() const { return "STAssemble"; }

  virtual void Initialize(Vector &u) const {Exit("Not Implemented!")}

  virtual void RHS(double t, Vector &rhs) const override {};

  virtual void MassMatrix(Matrix &massMatrix) const {};

  virtual void SystemMatrix(Matrix &systemMatrix) const {};

  virtual void System(Matrix &, Vector &) const {};

  virtual void RHS(Vector &) const {}

  virtual void get_projected_exact_solution(Vector &) const {};

  virtual void get_exact_solution(Vector &) const {};

  virtual void MakeAdjointMatrix(Matrix &B) const {
    SystemAddDoubleD(B);
    B.Transpose();
    B *= -1.;
  }

  virtual void SystemAddDoubleD(Matrix &) const {};

  virtual std::pair<double, double> DiscNorm(const Matrix &, const Vector &) const { return {}; };

  virtual double L1Error(const Vector &) const { return -1; };

  virtual double L1Norm(const Vector &) const { return -1; }

  virtual double DGError(const Vector &) const { return 0.0; };

  virtual double L2Error(const Vector &) const { return 0.0; };

  virtual double LInfError(const Vector &) const { return -1; };

  virtual double GNError(const Vector &) const { return 0.0; };

  virtual double VNorm(const Matrix &L, const Vector &) const { return -1.0; };

  virtual double MhalfNorm(const Vector &U) const { return -1; }

  virtual double MhalfInvLNorm(const Vector &U) const { return -1; }

  virtual double L2Norm(const Vector &) const { return 0.0; };

  virtual double L2ScalarProduct(const Vector &A, const Vector &B) const {
    Warning("Implement L2ScalarProduct for " + std::string(Name()));
    return 0.0;
  };

  virtual double EnergyNorm(const Vector &) const {Exit("Not Implemented!"); }

  virtual double DGSemiNorm(const Vector &u) const { return -1; }

  virtual double DGNormError(const Vector &u) { return -1; }

  virtual double DGNorm(const Vector &u, const Matrix &M) const { return -1; }

  virtual double L2SpaceNormAtEndTime(const Vector &u) const {
    double endTime = u.GetMesh().slice(u.GetMesh().slices() - 1);
    return L2SpaceNormAtTime(u, endTime);
  }

  virtual double L2SpaceNormAtTime(const Vector &u, double time) const { return -1; }

  virtual double L2SpaceNormAtTimeError(const Vector &u, double time) const { return -1; }

  virtual void DualRHS_linear(Vector &, const Vector &, Point, Point) const {}

  virtual void DualRHS_quadratic(Vector &, const Vector &, Point, Point) const {}

  virtual double Energy(const Vector &u, Point x0, Point x1) const { return 0.0; }

  virtual void Energy(const Vector &u, vector<double> &Energy) const { };

  virtual double Goal_Functional(const Vector &, Point, Point) const { return 0.0; }

  virtual void DualMatrix(Matrix &A, const Matrix &B) const {
    A = B;
    A.Transpose();
    A *= -1.;
  }


  virtual double MminushalfLuMinusF(const cell &c, const Vector &U) const { return 0.0; }

  virtual double MhalfLuMinusF(const Vector &U, const std::string &filename) const;

  virtual double MhalfLu_exactMinusF(LevelPair pair) const { return -1.0; }

  virtual double DualErrorEstimateJumpCell(const cell &c,
                                           const Vector &u,
                                           const Vector &u_star) const { return 0.0; };

  virtual double ResidualErrorEstimateJumpCell(const cell &c,
                                               const Vector &u) const { return 0.0; };

  virtual double L2ReliableResidualErrorEstimateJumpCell(const cell &c,
                                                         const Vector &U,
                                                         const Vector &conf) const { return 0.0; };

  virtual double DGErrorEstimateJumpCell(const cell &c,
                                         const Vector &u,
                                         const Vector &conf,
                                         double eta_res_c,
                                         double eta_conf_c) const { return 0.0; };

  virtual double DualErrorEstimateJump(const Vector &u,
                                       const Vector &u_star,
                                       Vector &Eta) const {
    Eta = 0;
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      Eta(c(), 0) = DualErrorEstimateJumpCell(c, u, u_star);
    }
    Eta.MakeAdditive();
    Eta.Accumulate();

    double max = Eta.Max();
    mout << " || Eta ||_infty = " << max << endl
         << " || Eta ||_2 = " << norm(Eta) << endl;
    return max;
  }

  virtual double ResidualErrorEstimateJump(const Vector &u, Vector &Eta) const {
    const int verbose = ErrorEstimatorVerbose;
    Eta = 0;
    for (cell c = u.cells(); c != u.cells_end(); ++c)
      Eta(c(), 0) = ResidualErrorEstimateJumpCell(c, u);
    double max = Eta.Max();
    vout(2) << " || Eta ||_infty = " << max << endl
            << " || Eta ||_2 = " << norm(Eta) << endl;
    return max;
  }

  virtual double L2ReliableResidualErrorEstimateJump(const Vector &U,
                                                     const Vector &conf,
                                                     Vector &Eta) const {
    const int verbose = ErrorEstimatorVerbose;
    Eta = 0;
    for (cell c = U.cells(); c != U.cells_end(); ++c)
      Eta(c(), 0) = L2ReliableResidualErrorEstimateJumpCell(c, U, conf);
    double max = Eta.Max();
    vout(2) << " || Eta ||_infty = " << max << endl
            << " || Eta ||_2 = " << norm(Eta) << endl;
    return max;
  }

  virtual double DGErrorEstimateJump(const Vector &u,
                                     Vector &Eta,
                                     const Vector &conf,
                                     const Vector &residual_ee,
                                     const Vector &conf_ee) const {
    const int verbose = ErrorEstimatorVerbose;
    Eta = 0;
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      Eta(c(), 0) = DGErrorEstimateJumpCell(c,
                                            u,
                                            conf,
                                            residual_ee(c(), 0),
                                            conf_ee(c(), 0));
    }
    double max = Eta.Max();
    vout(2) << " || Eta ||_infty = " << max << endl
            << " || Eta ||_2 = " << norm(Eta) << endl;
    return max;
  }

  virtual double Energy(const Vector &U) const {
    vector<double> energy_per_slice(U.GetMesh().steps());
    double energy = 0.0;
    for (int i = 0; i < U.GetMesh().steps(); ++i) {
      energy += energy_per_slice[i];
    }
    return energy;
  }

  virtual bool HasConformingLagrangeInterpolation() {
    return false;
  }

  virtual void
  ConformingLagrangeInterpolation(const Vector &u, Vector &conf, int degree = 1) const {
    Warning("Implement ConformingLagrangeInterpolation.");
    conf = 0.0;
  }


  virtual void ConformingInterpolation(const Vector &u, Vector &conf) const {
    Exit("Implement ConformingInterpolation.");
  };

  virtual void ConformingProjectionWithPenalty(const Vector &u,
                                               Vector &u_conf,
                                               Vector &r_conf) const {
    Exit("Implement ConformingProjectionWithPenalty.");
  };

  void AdaptQuadrature(LevelPair levels, int degree, bool adaptCellQuad) {
    bool AdaptQuadrature = false;
    Config::Get("AdaptQuadrature", AdaptQuadrature);
    if (AdaptQuadrature || !adaptCellQuad) {
      getDisc().adaptQuadrature(levels, degree, adaptCellQuad);
    }
  }

  void
  AdaptPolynomialDegrees(const Vector &Eta, double theta, double theta_min,
                         const std::string &refine_by) {
    //increase_poly_degs(getDisc(), Eta, theta, theta_min, refine_by);
    std::unordered_map<Point, DegreePair> polyDist = createPolynomialDistribution(
        getDisc(), Eta, theta, theta_min, refine_by);
    mout << "AdaptPolynomialDegrees:" << Eta.Level() << endl;
    getDisc().CreateAdaptivityLevel(polyDist, Eta.Level().NextInAdaptivity());
    //getDisc().communicate(Eta.Level().NextInAdaptivity());
    //RedistributeCells(Eta.Level());
  }

  virtual void RedistributeCells(LevelPair levels) = 0;

  void AdaptPolynomialDegreesUniformly(LevelPair levels) {
    increase_poly_degs_uniform(getDisc());
    getDisc().communicate(levels);
  }

  virtual bool HasMatrixFreeOperator() const { return false; }

  virtual std::unique_ptr<Operator> GetMatrixFreeOperator() const {
    THROW("Not implemented for: " + std::string(Name()))
  }

  virtual void PlotSolution(const Vector &u) {}

  virtual void PlotDualSolution(const Vector &dual_u, const Vector &eta, int run) {}


  void PlotErrorEstimator(const Vector &vector, std::string filename);
};

std::unique_ptr<STAssemble> CreateSTAssemble(const string &modelName,
                                             const Meshes &meshes,
                                             DegreePair degree,
                                             const string &probName);

template<typename FaceElement>
inline int find_q_id(const Point &Qf_c, const FaceElement &felem_1) {
  for (int q1 = 0; q1 < felem_1.nQ(); ++q1) {
    if (Qf_c == felem_1.QPoint(q1)) {
      return q1;
    }
  }
  std::stringstream s;
  s << "Looking for Q Point: " << Qf_c << endl;
  for (int q1 = 0; q1 < felem_1.nQ(); ++q1) {
    s << q1 << " " << felem_1.QPoint(q1) << endl;
  }
  std::cout << s.str();
  THROW("Error, should not get here.");
  Exit("Error, should not get here.");
}




#endif
