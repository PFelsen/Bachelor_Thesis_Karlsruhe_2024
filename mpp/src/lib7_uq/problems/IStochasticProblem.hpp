#ifndef ISTOCHASTICPROBLEM_HPP
#define ISTOCHASTICPROBLEM_HPP

#include <utility>

#include "CirculantEmbedding.hpp"
#include "HybridFluxGenerator.hpp"
#include "IProblem.hpp"
#include "KLExpansionGenerator.hpp"
#include "Random.hpp"
#include "SparseGridGenerator.hpp"


class IStochasticProblem : virtual public IProblem {
protected:
  int verbose = 1;

  SampleID currentID;

  virtual void drawSample(const SampleID &id) = 0;

public:
  explicit IStochasticProblem() {
    Config::Get("ProblemVerbose", verbose);
  }

  void DrawSample(const SampleID &id) {
    currentID = id;
    drawSample(currentID);
  }

  const SampleID &CurrentID() const { return currentID; }

  virtual double SampleWeight() const { return 1.0; }

  virtual double SumOfWeights() { return 1.0; }

  virtual int NumOfSamples() const { return 0; };

  std::string Name() const override { return currentID.IdString(); }

  virtual void LoadRHS(const Vector &u, bool withData) {};

  virtual double FunctionEvaluation() { return 0.0; }
};

#endif
