#ifndef SPACETIME_STCALLBACKSTRATEGY_HPP
#define SPACETIME_STCALLBACKSTRATEGY_HPP

#include "LinearSolver.hpp"
#include "ErrorEstimator.hpp"
#include "STAssemble.hpp"
#include "SpaceTimePlotting.hpp"

class ErrorEstimatorCallbackStrategy : public CallbackStrategy {

  const ErrorEstimator &ee;
  const STAssemble &assemble;

  bool plot = false;
  bool calcDGError = false;
  bool printChange = false;
  double epsilon1 = 0.0;
  double epsilon2 = 0.0;

  int iteration = 0;
  double e_error = 1;
  double e_change = 0.0;

  int firstActivationIteration = -1;

public:
  ErrorEstimatorCallbackStrategy(const ErrorEstimator &ee, const STAssemble &assemble)
  : ee(ee), assemble(assemble) {
    Config::Get("EECallbackPlot", plot);
    Config::Get("EECallbackEpsilon1", epsilon1);
    Config::Get("EECallbackEpsilon2", epsilon2);
    Config::Get("EECallbackCalculateDGError", calcDGError);
    Config::Get("EECallbackPrintChange", printChange);
  }

  bool isActive(int restartIter, int totalIter, double residual) override {
    iteration = totalIter;
    if (residual < epsilon1){
      if (firstActivationIteration < 0) {
        firstActivationIteration = totalIter;
      }
      return true;
    }
    return false;
  }
  double checkCurrentSolution(const Matrix &M, const Vector &u) override {
    Vector error = ee.EstimateError(M, u, 0);
    if (plot) {
      printVTK_Eta(u, assemble, "Eta_Iter" + std::to_string(iteration));
    }
    if (calcDGError) {
      double dg_error = assemble.DGError(u);
      mout << "DG(" << iteration << ")= " << dg_error << endl;
    }
    double current_e_error = norm(error);

    e_change = e_error / current_e_error;
    if(printChange) {
      mout << "E/E(" << iteration << ")= " << e_change << endl;
    }
    e_error = current_e_error;

    return e_error;
  }

  bool isFinished(double residual) const override{
    return firstActivationIteration != iteration && residual < epsilon2;
  }

};

class STPlotCallbackStrategy : public CallbackStrategy {

  const STAssemble &assemble;
  int iteration = 0;
  int plotFrequency = 1;

public:
  STPlotCallbackStrategy(const STAssemble &assemble)
      : assemble(assemble) {
    Config::Get("PlotCallbackFrequency", plotFrequency);
    globalPlotBackendCreator = std::make_unique<DefaultSpaceTimeBackend>;
  }

  bool isActive(int restartIter, int totalIter, double residual) override {
    iteration = totalIter;
    return iteration % plotFrequency == 0;
  }
  double checkCurrentSolution(const Matrix &M, const Vector &u) override {
    assemble.PlotSingleSolution(u, "Iter" + std::to_string(iteration));
    return -1;
  }

  bool isFinished(double residual) const override{
    return false;
  }

};

#endif //SPACETIME_STCALLBACKSTRATEGY_HPP
