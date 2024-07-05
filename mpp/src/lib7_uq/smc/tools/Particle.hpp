#ifndef PARTICLE_H
#define PARTICLE_H

#include "Metropolis.hpp"

struct ParticleConfig {
  RVector measurement = RVector({1});
  RMatrix precision = RMatrix({{1}}); // TODO PRECISION_MATRIX
  PDESolverConfig pde_model_config = PDESolverConfig();

  ParticleConfig();
  ParticleConfig WithMeasurement(std::vector<double> measurement);
  ParticleConfig WithPrecision(std::vector<double> precision);
  ParticleConfig WithPDEModelConfig(PDESolverConfig conf);
  RMatrix ToRMatrix(std::vector<double> flat_matrix);
};

class Particle {
private:
  int verbose = 0;
  std::shared_ptr<Measurement> measurement;
  std::shared_ptr<EllipticPDEModel> pde_model;
  std::shared_ptr<MetropolisKernel> markov_kernel;
public:
  Particle(ParticleConfig conf);
  void mutate(double temperature = 1);
  double GetPotential();
  void copyState(std::shared_ptr<Particle> other);
  RVector GetState(); // TODO only valid for double input
  MetropolisKernel &GetKernel();

  EllipticPDEModel &GetModel() { return *pde_model; };

  RVector GetObservation();
};

#endif // PARTICLE_H