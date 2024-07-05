#include "Particle.hpp"

// ------ Particle ------

Particle::Particle(ParticleConfig conf) {
  measurement = std::make_shared<Measurement>(conf.measurement, conf.precision);
  pde_model = CreatePDEModelShared(conf.pde_model_config);
  markov_kernel = std::make_shared<MetropolisKernel>(pde_model, measurement,
                                                     conf.pde_model_config.proposal_config);
}

void Particle::mutate(double temperature) {
  mout.StartBlock("Mutate Particle");
  vout(1) << "Start with temperature " << temperature << endl;
  markov_kernel->Step(temperature);
  mout.EndBlock(verbose > 0);
}

double Particle::GetPotential() { return markov_kernel->potential; }

void Particle::copyState(std::shared_ptr<Particle> other) {
  markov_kernel->potential = other->markov_kernel->potential;
  markov_kernel->GetProposalKernel().SetState(other->markov_kernel->GetProposalKernel().GetState());
  measurement->mean = other->measurement->mean;
  measurement->precision = other->measurement->precision;
  pde_model->GetProblem().GetInput().GetRandomField().SetParameterVector(
      other->pde_model->GetProblem().GetInput().GetRandomField().GetParameterVector());
}

RVector Particle::GetState() { return markov_kernel->GetProposalKernel().GetState(); }

RVector Particle::GetObservation() { return markov_kernel->GetObservation(); }

MetropolisKernel &Particle::GetKernel() { return *markov_kernel; }

// ------ ParticleConfig ------

ParticleConfig::ParticleConfig() {
  std::vector<double> vector_measurement;
  std::vector<double> vector_precision;
  Config::Get("Measurement", vector_measurement);
  Config::Get("Precision", vector_precision);
  if (!vector_measurement.empty()) measurement = RVector(vector_measurement);
  if (!vector_precision.empty()) precision = RMatrix(ToRMatrix(vector_precision));
}

ParticleConfig ParticleConfig::WithMeasurement(std::vector<double> measurement) {
  this->measurement = RVector(measurement);
  return *this;
}

ParticleConfig ParticleConfig::WithPrecision(std::vector<double> flat_precision) {
  this->precision = RMatrix(ToRMatrix(flat_precision));
  return *this;
}

ParticleConfig ParticleConfig::WithPDEModelConfig(PDESolverConfig conf) {
  this->pde_model_config = conf;
  return *this;
}

RMatrix ParticleConfig::ToRMatrix(std::vector<double> flat_matrix) {
  size_t n = static_cast<size_t>(std::sqrt(flat_matrix.size()));
  std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      matrix[i][j] = flat_matrix[i * n + j];
    }
  }
  return RMatrix(matrix);
}
