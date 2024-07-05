#ifndef SAMPLERSELECTION_H
#define SAMPLERSELECTION_H

#include "CirculantEmbedding.hpp"
#include "Random.hpp"
#include "SampleGenerator.hpp"

class Sampler {
public:
  virtual RVector DrawSample(int length) = 0;

  virtual double Prob(RVector u) = 0;
};

class NormalSampler : public Sampler {
private:
  double mean, variance;
public:
  NormalSampler(PriorSamplerConfig conf) {
    mean = conf.mean;
    variance = conf.variance;
  }

  RVector DrawSample(int length) override;

  double Prob(RVector u) override;
};

class UniformSampler : public Sampler {
private:
  double mean, bound;
public:
  UniformSampler(PriorSamplerConfig conf) {
    bound = conf.variance;
    mean = conf.mean;
  }

  RVector DrawSample(int length) override;

  double Prob(RVector u) override;
};

class TestSampler : public Sampler {
public:
  TestSampler(PriorSamplerConfig conf) {}

  RVector DrawSample(int length) override;

  double Prob(RVector u) override;
};

Sampler *CreateSampler(PriorSamplerConfig conf);

std::unique_ptr<Sampler> CreateSamplerUnique(PriorSamplerConfig conf);

#endif // SAMPLERSELECTION_H