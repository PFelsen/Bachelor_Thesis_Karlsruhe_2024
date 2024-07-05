#ifndef SAMPLEGENERATOR_HPP
#define SAMPLEGENERATOR_HPP

#include <utility>

#include "Config.hpp"
#include "Sample.hpp"

#include "Meshes.hpp"
#include "RVector.hpp"
#include "CVector.hpp"
#include "RMatrix.hpp"
#include "CMatrix.hpp"
#include "RTensor.hpp"
#include "CTensor.hpp"


typedef std::complex<double> Complex;

static double abs_compl(const Complex &complex) {
  return sqrt(pow(complex.real(), 2) + pow(complex.imag(), 2));
}

template<typename T>
class SampleGenerator {
protected:
  int verbose = 0;

  virtual void drawSample(const SampleID &id) = 0;

public:

  explicit SampleGenerator() {
    Config::Get("GeneratorVerbose", verbose);
  }

  virtual ~SampleGenerator() = default;

  void DrawSample(const SampleID &id) {
    mout.StartBlock(Name());
    vout(3) << id << endl;
    drawSample(id);
    mout.EndBlock(verbose < 3);
  }

  virtual string Name() const = 0;

  virtual T EvalSample() const {
    Exit("Not implemented")
  }

  virtual T EvalSample(const Point &x) const {
    Exit("Not implemented")
  }

  virtual T EvalSample(int n, const Point &x) const {
    Exit("Not implemented")
  }

  virtual T EvalSample(double t, const Point &x) const {
    Exit("Not implemented")
  }

  virtual T EvalSample(const Cell &c) const {
    Exit("Not implemented")
  }

  virtual T EvalSample(int n, const Cell &c) const {
    Exit("Not implemented")
  }

  virtual T EvalSample(double t, const Cell &c) const {
    Exit("Not implemented")
  }

  T DrawAndEvalSample(const SampleID &id) {
    DrawSample(id);
    return EvalSample();
  }

  T DrawAndEvalSample(const SampleID &id, const Point &x) {
    DrawSample(id);
    return EvalSample(x);
  }

  T DrawAndEvalSample(const SampleID &id, int n, const Point &x) {
    DrawSample(id);
    return EvalSample(n, x);
  }

  T DrawAndEvalSample(const SampleID &id, double t, const Point &x) {
    DrawSample(id);
    return EvalSample(t, x);
  }

  T DrawAndEvalSample(const SampleID &id, const Cell &c) {
    DrawSample(id);
    return EvalSample(c);
  }

  T DrawAndEvalSample(const SampleID &id, int n, const Cell &c) {
    DrawSample(id);
    return EvalSample(n, c);
  }

  T DrawAndEvalSample(const SampleID &id, double t, const Cell &c) {
    DrawSample(id);
    return EvalSample(t, c);
  }
};

#endif //SAMPLEGENERATOR_HPP
