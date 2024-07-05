#ifndef MLUQ_MAINSGD_HPP
#define MLUQ_MAINSGD_HPP


#include "ScalarElement.hpp"
#include "MultiLevelEstimator.hpp"
#include "SingleLevelEstimator.hpp"
#include "Plotting.hpp"

class StochasticGradientDescent {
private:
  //problem parameter
  int verbose = 1;
  std::string Quantity = "L2";
  std::string modelName = "LagrangeElliptic";
  std::string problemName = "OCStochasticLaplace2D";
  std::string protocolName = "OCStochasticLaplace2D";
  int level = 4;
  int plotting = 0;
  int iterationCounter = 0;

  //OCP parameter
  double u_a, u_b, targetCost;

  // SGD parameter
  double alpha;
  double stepsize;

  //breakpoint arguments
  int maxSteps = 100;
  std::string breakpointCriteria;


  //ADAM parameter
  std::string descentType;
  double beta1 = 0.9;
  double beta2 = 0.99;
  std::unique_ptr<Vector> m_k, v_k;
  double epsilonADAM = 0.00000001;
  double gammaADAM;

  // Evaluation parameter
  std::string secondaryQOI;
  bool overkill;//Todo: Rename

  // step size parameter
  std::string stepsizeRule;
  double theta, nu, gamma;

  // others
  std::string name = "StochasticGradientDescent";
  std::unique_ptr<double> sum_of_stepsizes;

public:
  StochasticGradientDescent() {
    Config::Get("SGDVerbose", verbose);
    Config::Get("Quantity", Quantity);
    Config::Get("Model", modelName);
    Config::Get("Problem", problemName);
    Config::Get("Protocol", protocolName);
    Config::Get("level", level);
    Config::Get("VtuPlot", plotting);

    Config::Get("u_a", u_a);
    Config::Get("u_b", u_b);
    Config::Get("targetCost", targetCost);

    Config::Get("maxSteps", maxSteps);
    Config::Get("alpha", alpha);

    Config::Get("descentType", descentType);
    Config::Get("gammaADAM", gammaADAM);

    Config::Get("secondaryQOI", secondaryQOI);
    Config::Get("overkill", overkill);

    Config::Get("stepsizeRule", stepsizeRule);
    Config::Get("theta", theta);
    Config::Get("nu", nu);
    Config::Get("gamma", gamma);

  }

  double Compute();

  void Method() {
    Compute();
  }



  void AveragingSolution(const SampleSolution &u_k, SampleSolution &u_T) const;

  double Stepsize() const {return stepsize; };

  void SetStepsize();

  SampleSolution Update(SampleSolution &control, const SampleSolution &descent) const;

  SampleSolution DescentEstimator(SampleSolution &gradient) const;

  SampleSolution AdmissibleMapping(SampleSolution &control) const;

  void PlotControl(const SampleSolution &control) const;

  int Iteration() const {return iterationCounter; };

  void UpdateIteration() {iterationCounter += 1; };

  bool Breakpoint() const;

  void UpdateFunctionQOI(vector<double> &j_k_l2_norm, const SampleID &id, std::shared_ptr<MultiSampleFEM> &msFEM, const SampleSolution &control) const;
};

#endif //MLUQ_MAINSGD_HPP
