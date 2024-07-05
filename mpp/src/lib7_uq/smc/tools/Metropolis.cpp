#include "Metropolis.hpp"
#include "Random.hpp"

MetropolisKernel::MetropolisKernel(std::shared_ptr<EllipticPDEModel> model,
                                   std::shared_ptr<Measurement> measurement, ProposalConfig conf) :
    measurement(measurement), pdeSolver(model), observation(measurement->mean.size()),
    new_observation(measurement->mean.size()) {

  proposal_kernel = CreateProposalGeneratorUnique(conf, pdeSolver->GetProblem()
                                                            .GetInput()
                                                            .GetRandomField()
                                                            .GetParameterVector());
  adaptive = conf.adaptive;

  pdeSolver->Run(observation);
  potential = measurement->LogLikelihood(observation);
  new_observation = observation;
}

void MetropolisKernel::Step(double temperature) {
  if (adaptive) { proposal_kernel->AdaptStepSize(); }
  proposal_kernel->Propose();
  pdeSolver->Run(new_observation);
  double new_potential = measurement->LogLikelihood(new_observation);

  double acceptance_probability =
      proposal_kernel->AcceptanceProbability(potential, new_potential, temperature);
  if (acceptance_probability < 1) {
    double random_draw = Random::Uniform(0, 0, 1); // TODO COMM_SPLIT
    if (acceptance_probability < random_draw) { return; }
  }
  potential = new_potential;
  proposal_kernel->Accept();
  observation = new_observation;
}

RVector MetropolisKernel::GetObservation() { return observation; }

ProposalGenerator &MetropolisKernel::GetProposalKernel() { return *proposal_kernel; }

// ------ Measurement ------

Measurement::Measurement(RVector measurement, RMatrix precision) {
  this->mean = measurement;
  this->precision = precision;
}

double Measurement::LogLikelihood(RVector solution) {
  RVector difference;
  double potential;

  difference = solution - mean;
  potential = 0.5 * difference * precision * difference;
  return potential;
}
