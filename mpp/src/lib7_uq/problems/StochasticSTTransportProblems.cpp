#include "StochasticSTTransportProblems.hpp"

class StochasticPollutionRiemann2D : public IStochasticSTTransportProblem {
  HybridFaceNormalFluxGenerator scalarGenerator;

  HybridCellFluxGenerator vectorFieldGenerator;

  double d = 1.0 / 16.0;

protected:
  void drawSample(const SampleID &id) override {};

public:
  StochasticPollutionRiemann2D() : IStochasticSTTransportProblem(2),
        IProblem("STSquare"),
        scalarGenerator(HybridFaceNormalFluxGenerator()),
        vectorFieldGenerator(HybridCellFluxGenerator(scalarGenerator)) {}

  VectorField Flux(const LevelPair &level, const Cell &c, const Point &x) const override {
    return this->vectorFieldGenerator.EvalSample(c);
  }

  Scalar FaceNormalFlux(const LevelPair &level, const Cell &c, int face, const VectorField &N,
                    const Point &x) const override {
    return this->scalarGenerator.EvalSample(face, c);
  }
};

IStochasticSTTransportProblem *
CreateTransportStochasticProblem(const string &problemName) {
  IStochasticSTTransportProblem *problem = nullptr;
  if (problemName == "StochasticPollutionRiemann2D")
    problem = new StochasticPollutionRiemann2D();
  return problem;
}