#ifndef MLUQ_Metropolis_H
#define MLUQ_Metropolis_H

#include "PDEModel.hpp"

class Measurement {
public:
  RVector mean;

  RMatrix precision;

  Measurement(RVector measurement, RMatrix precision);

  double LogLikelihood(RVector solution);
};

class MetropolisKernel {
private:
  std::shared_ptr<EllipticPDEModel> pdeSolver;

  std::shared_ptr<Measurement> measurement;

  std::unique_ptr<ProposalGenerator> proposal_kernel;

  RVector observation;

  RVector new_observation;

  bool adaptive;
public:
  double potential;
  MetropolisKernel(std::shared_ptr<EllipticPDEModel> model, std::shared_ptr<Measurement> measurement,
                   ProposalConfig conf);
  void Step(double temperature = 1);
  ProposalGenerator &GetProposalKernel();
  RVector GetObservation();
};

#endif // MLUQ_Metropolis_H
