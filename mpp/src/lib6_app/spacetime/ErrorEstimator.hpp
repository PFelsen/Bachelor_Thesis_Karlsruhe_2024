#ifndef SPACETIME_ERRORESTIMATOR_HPP
#define SPACETIME_ERRORESTIMATOR_HPP

#include "LinearSolver.hpp"
#include "STAssemble.hpp"

#include <memory>

class ErrorEstimator {
protected:
  STAssemble &assemble;
  LinearSolver &solver;
  int verbose = 0;
public:
  ErrorEstimator(STAssemble &assemble, LinearSolver &solver) : assemble(assemble),
                                                               solver(solver) {
    Config::Get("ErrorEstimatorVerbose", verbose);
  };

  virtual Vector EstimateError(const Matrix &B, const Vector &U, int run) const = 0;

  virtual ~ErrorEstimator() = default;
};

std::unique_ptr<ErrorEstimator>
CreateErrorEstimator(const std::string &name, STAssemble &assemble, LinearSolver &solver);



class ResidualErrorEstimator : public ErrorEstimator {
public:
  Vector EstimateError(const Matrix &B, const Vector &U, int run) const override {
    mout.StartBlock("ResidualErrorEstimate");
    Vector Eta(0.0, U);
    assemble.ResidualErrorEstimateJump(U, Eta);
    mout.EndBlock(verbose <= 0);
    assemble.PlotErrorEstimator(Eta, "Residual_Est_Error"+ std::to_string(run));
    return Eta;
  }

  ResidualErrorEstimator(STAssemble &assemble, LinearSolver &solver)
      : ErrorEstimator(assemble, solver) {}
};

class L2ReliableResidualErrorEstimator : public ErrorEstimator {
public:
  Vector EstimateError(const Matrix &B, const Vector &U, int run) const override {
    mout.StartBlock("ConfErrorEstimate");
    Vector Eta(0.0, U);
    Vector conf(0.0, U);
    assemble.ConformingLagrangeInterpolation(U, conf);
    assemble.L2ReliableResidualErrorEstimateJump(U, conf, Eta);
    mout.EndBlock(verbose <= 0);;
    assemble.PlotErrorEstimator(Eta, "L2_Est_Error"+ std::to_string(run));
    return Eta;
  }


  L2ReliableResidualErrorEstimator(STAssemble &assemble, LinearSolver &solver)
      : ErrorEstimator(assemble, solver) {}
};


class DGErrorEstimator : public ErrorEstimator {
public:
  Vector EstimateError(const Matrix &B, const Vector &U, int run) const override {
    mout.StartBlock("DGErrorEstimate1");

    L2ReliableResidualErrorEstimator ee_reliable(assemble, solver);
    ResidualErrorEstimator ee_residual(assemble, solver);
    Vector eta_reliable = ee_reliable.EstimateError(B, U, run);
    eta_reliable.MakeAdditive();
    Vector eta_residual = ee_residual.EstimateError(B, U, run);
    eta_residual.MakeAdditive();
    Vector Eta = EstimateError(B, U, run, eta_residual, eta_reliable);
    mout.EndBlock(verbose <= 0);

    return Eta;
  }

  virtual Vector EstimateError(const Matrix &B, const Vector &U, int run,
                               const Vector &residual_ee, const Vector &conf_ee) const{
    mout.StartBlock("DGErrorEstimate2");
    Vector conf(0.0, U);
    assemble.ConformingLagrangeInterpolation(U, conf);
    Vector Eta(0.0, U);
    assemble.DGErrorEstimateJump(U, Eta, conf, residual_ee, conf_ee);
    mout.EndBlock(verbose <= 0);
    assemble.PlotErrorEstimator(Eta, "DG_Est_Error"+ std::to_string(run));
    return Eta;
  }


  DGErrorEstimator(STAssemble &assemble, LinearSolver &solver)
      : ErrorEstimator(assemble, solver) {}
};


class DualErrorEstimator : public ErrorEstimator {
protected:
  void fill_Dual_RHS(Vector &Dual_RHS, const Vector &U) const {
    Point roi_min;
    Point roi_max;
    std::string dual_functional = "none";
    Config::Get("roi_min", roi_min);
    Config::Get("roi_max", roi_max);
    Config::Get("dual_functional", dual_functional);
    if (dual_functional == "linear")
      assemble.DualRHS_linear(Dual_RHS, U, roi_min, roi_max);
    if (dual_functional == "quadratic")
      assemble.DualRHS_quadratic(Dual_RHS, U, roi_min, roi_max);
    Dual_RHS.Collect();
  }

public:
  DualErrorEstimator(STAssemble &assemble, LinearSolver &solver) : ErrorEstimator(assemble,
                                                                                  solver) {}

  Vector EstimateError(const Matrix &M, const Vector &U, int run) const override  {
    Vector Dual_U(0.0, U);
    Vector Dual_RHS(0.0, U);
    fill_Dual_RHS(Dual_RHS, U);
    Matrix B(M);
    B = M;

    assemble.MakeAdjointMatrix(B);
    solver.GetPreconditioner().transpose();
    solver(B, true);
    Dual_U = solver * Dual_RHS;
    solver.GetPreconditioner().transpose();

    Vector Eta(0.0, Dual_U);
    assemble.DualErrorEstimateJump(U, Dual_U, Eta);

    assemble.PlotDualSolution(Dual_U, Eta, run);
    return Eta;
  }
};




#endif //SPACETIME_ERRORESTIMATOR_HPP
