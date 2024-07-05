#ifndef STTRANSPORTMAIN_HPP
#define STTRANSPORTMAIN_HPP

#include "STDGTransportAssemble.hpp"
#include "TransportProblems.hpp"
#include "MeshesCreator.hpp"
#include "PDESolver.hpp"
#include "LinearSolver.hpp"
#include "ErrorEstimator.hpp"
#include "AdaptiveConfig.hpp"

class STTransportMain : public PDESolver<ITransportProblem> {
private:

  std::unique_ptr<STDGTransportAssemble> assemble = nullptr;

  std::unique_ptr<Matrix> M = nullptr;
  std::unique_ptr<LinearSolver> solver = nullptr;


  std::vector<double> CreateTimestepsFromConfig(const PDESolverConfig &mc) {
    double T = mc.T;
    double dt = mc.dt * pow(2, mc.level);

    if (ceil(T / mc.dt) != pow(2, mc.level) * ceil(T / dt)) {
      THROW("Timesteps do not fit for coarse level in spacetime");
    }

    vector<double> timesteps;
    for (size_t i = 0; i <= ceil(T / dt); i++) {
      timesteps.push_back(i * dt);
    }
    return timesteps;
  }

protected:
  void run(Solution &solution) const override;

  void createAssemble(std::shared_ptr<ITransportProblem> problem) override;

  void computeValues(Solution &solution) const override;

public:
  explicit STTransportMain(const PDESolverConfig &conf) : GenericPDESolver(conf) { }

  STTransportMain() : STTransportMain(PDESolverConfig{}) {}

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return assemble->GetSharedDisc();
  }

  std::string Name() const override {
    return "STTransport";
  }

  void PrintValues(const Solution &solution) override;

  void PlotVtu(const Vector &u) const override;

  void plotVtu(Solution &solution) const override {PlotVtu(solution.vector);};

  bool AdaptiveStep(Vector &u, int step);

  void EstimateErrorAndApplyAdaptivity(const Vector &U, AdaptiveConfig conf) {
    mout.StartBlock("Adaptivity");
    auto errorEstimator = CreateErrorEstimator(conf.data.errorEstimator,
                                               *assemble,
                                               *solver);
    Vector Eta = errorEstimator->EstimateError(*M, U, U.Level().adaptivityLevel);
    double theta = conf.data.theta;
    double theta_min = conf.data.theta_min;
    double theta_factor = conf.data.theta_factor;
    std::string refine_by = std::ref(conf.data.refine_by);
    theta *= pow(theta_factor, Eta.Level().adaptivityLevel);
    assemble->AdaptPolynomialDegrees(Eta, theta, theta_min, refine_by);
    assemble->AdaptQuadrature(Eta.Level().NextInAdaptivity(), -1, false);
    mout.EndBlock();
  }

};

#endif