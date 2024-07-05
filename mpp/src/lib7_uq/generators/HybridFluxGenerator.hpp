#ifndef HYBRIDFLUXGENERATOR_HPP
#define HYBRIDFLUXGENERATOR_HPP

#include "CirculantEmbedding.hpp"
#include "EllipticProblems.hpp"

class EllipticPDESolver;

class HybridFaceNormalFluxGenerator : public SampleGenerator<Scalar> {
private:
  void drawSample(const SampleID &id) override;

public:
  EllipticPDESolver *pdeSolver = nullptr;

  std::shared_ptr<IEllipticProblem> problem = nullptr;

  std::unique_ptr<SampleSolution> solutionFaceFlux;

  std::unique_ptr<SampleSolution> solutionFaceValues;

  HybridFaceNormalFluxGenerator() : SampleGenerator() {}

  ~HybridFaceNormalFluxGenerator() override;

  string Name() const override { return "HybridFaceNormalFluxGenerator"; }

  Scalar EvalSample(int face, const Cell &c) const override;
};

class HybridCellFluxGenerator : public SampleGenerator<VectorField> {
  const HybridFaceNormalFluxGenerator &generator;

  void drawSample(const SampleID &id) override {};

public:
  explicit HybridCellFluxGenerator(HybridFaceNormalFluxGenerator &generator) :
      SampleGenerator(), generator(generator) {}

  VectorField EvalSample(const Cell &c) const override;

  string Name() const override { return "HybridCellFluxGenerator"; }
};

#endif //HYBRIDFLUXGENERATOR_HPP
