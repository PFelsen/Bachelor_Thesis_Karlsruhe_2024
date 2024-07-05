#include "Resampling.hpp"

// ------ Resampling ------

void Resampling::resample(std::vector<double> &weights,
                          std::vector<std::shared_ptr<Particle>> &particles) {
  std::vector<int> indices = getIndices(weights);

  prepare_particle_copying(indices);

  for (int m = 0; m < weights.size(); m++) {
    if (indices[m] != m) { particles[m]->copyState(particles[indices[m]]); }
    weights[m] = 1.0 / weights.size();
  }
}

void Resampling::prepare_particle_copying(std::vector<int> &indices) {
  // Avoids copying conflicts:
  //  if j is in indices, then indices[j] = j, so indices[j] will not get overwritten before
  //  potential use
  int index;
  for (int m = 0; m < indices.size(); m++) {
    index = indices[m];
    if (indices[index] != index) {
      indices[m] = indices[index];
      indices[index] = index;
    }
  }
}

std::vector<int> Systematic::getIndices(const std::vector<double> &weights) {
  std::vector<int> indices(weights.size());
  double v = Random::Uniform(0, 0, 1.0 / weights.size()); // TODO COMM_SPLIT
  double acc_weight = weights[0];

  int index = 0;
  for (int m = 0; m < weights.size(); m++) {
    while (acc_weight < v) {
      acc_weight = acc_weight + weights[index];
      index++;
    }
    indices[m] = index;
    v = v + 1.0 / weights.size();
  }
  return indices;
}

std::vector<int> Multinomial::getIndices(const std::vector<double> &weights) {
  std::vector<int> indices(weights.size());
  std::vector<double> acc_weights(weights.size());
  std::vector<double> rand_thresholds(weights.size());

  for (int m = 0; m < weights.size(); m++) {
    rand_thresholds[m] = Random::Uniform(0, 0, 1); // TODO COMM_SPLIT
  }

  acc_weights[0] = weights[0];
  for (int m = 1; m < weights.size(); m++) {
    acc_weights[m] = acc_weights[m - 1] + weights[m];
  }

  int index;
  for (int m = 0; m < weights.size(); m++) {
    index = 0;
    while (acc_weights[index] < rand_thresholds[m]) {
      index++;
    }
    indices[m] = index;
  }
  return indices;
}

// ------ CreateResampling ------

Resampling *CreateResampling(std::string resampling_scheme) {
  if (resampling_scheme == "Multinomial") { return new Multinomial(); }
  if (resampling_scheme == "Systematic") { return new Systematic(); }

  Exit(resampling_scheme
       + " not found, possible resampling schemes are Multinomial and Systematic.")
};

std::unique_ptr<Resampling> CreateResamplingUnique(std::string resampling_scheme) {
  return std::unique_ptr<Resampling>(CreateResampling(resampling_scheme));
};
