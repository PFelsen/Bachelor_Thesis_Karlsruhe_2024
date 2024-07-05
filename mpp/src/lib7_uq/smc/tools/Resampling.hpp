#ifndef RESAMPLINGSCHEME_H
#define RESAMPLINGSCHEME_H

#include "Particle.hpp"
#include "Random.hpp"

class Resampling {
private:
  void prepare_particle_copying(std::vector<int> &indices);
public:
  void resample(std::vector<double> &weights, std::vector<std::shared_ptr<Particle>> &particles);

  virtual std::vector<int> getIndices(const std::vector<double> &weights) = 0;
};

class Systematic : public Resampling {
public:
  std::vector<int> getIndices(const std::vector<double> &weights);
};

class Multinomial : public Resampling {
public:
  std::vector<int> getIndices(const std::vector<double> &weights);
};

Resampling *CreateResampling(std::string resampling_scheme);
std::unique_ptr<Resampling> CreateResamplingUnique(std::string resampling_scheme);

#endif // RESAMPLINGSCHEME_H