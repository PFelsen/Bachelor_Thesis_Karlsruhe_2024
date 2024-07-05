#include "SamplerSelection.hpp"

RVector NormalSampler::DrawSample(int length) {
  RVector draw = Random::Normal(0, length); // TODO COMM_SPLIT
  draw = draw * sqrt(variance);
  return draw + RVector(mean, length);
}

double NormalSampler::Prob(RVector u) {
  const RVector difference = u - RVector(mean, u.size());
  return std::exp(-0.5 * difference.normSqr() / variance);
}

RVector UniformSampler::DrawSample(int length) {
  if (bound == 0) return RVector(mean, length);
  RVector draw = Random::Uniform(0, length, -bound, bound); // TODO COMM_SPLIT
  return draw + RVector(mean, length);
}

double UniformSampler::Prob(RVector u) {
  for (int k = 0; k < u.size(); k++)
    if (std::fabs(u[k]) > bound) return 0;
  return 1;
}

RVector TestSampler::DrawSample(int length) {
  RVector draw = RVector(0, length);
  draw[0] = 1;
  return draw;
}

double TestSampler::Prob(RVector u) { return 1; }

Sampler *CreateSampler(PriorSamplerConfig conf) {
  if (conf.generator == "Normal") { return new NormalSampler(conf); }
  if (conf.generator == "Uniform") { return new UniformSampler(conf); }
  if (conf.generator == "Test") { return new TestSampler(conf); }

  Exit(conf.generator + " not found")
};

std::unique_ptr<Sampler>
CreateSamplerUnique(PriorSamplerConfig conf) {
  return std::unique_ptr<Sampler>(CreateSampler(conf));
};

PriorSamplerConfig::PriorSamplerConfig() {
  if (!Config::IsInitialized()) return;

  Config::Get("SampleGenerator", generator);
  Config::Get("PriorMean", mean);
  Config::Get("PriorVariance", variance);
}

PriorSamplerConfig PriorSamplerConfig::WithGenerator(std::string generator) {
  this->generator = generator;
  return *this;
}

PriorSamplerConfig PriorSamplerConfig::WithMean(double mean) {
  this->mean = mean;
  return *this;
}

PriorSamplerConfig PriorSamplerConfig::WithVariance(double variance) {
  this->variance = variance;
  return *this;
}