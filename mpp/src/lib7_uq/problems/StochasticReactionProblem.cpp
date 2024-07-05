#include "StochasticReactionProblems.hpp"


inline Scalar u0(const Point &x) {
  double thickness = 1.0 / 8.0;
  double centerStart = 1 - 1.5 * thickness;
  if ((abs(centerStart - x[1]) < thickness / 2) &&
    (abs(0.5 - x[0]) < 3 * thickness))
    return 1;
  else return 0.0;
}

class StochasticReactionPollutionCosHat2D : public IStochasticReactionProblem {
  HybridFaceNormalFluxGenerator scalarGenerator;

  HybridCellFluxGenerator vectorFieldGenerator;

  double amplitude = 1.0;

  double cc = 6.0;

public:
  explicit StochasticReactionPollutionCosHat2D() :
    IProblem("UnitSquare"),
    scalarGenerator(HybridFaceNormalFluxGenerator()),
    vectorFieldGenerator(HybridCellFluxGenerator(scalarGenerator)) {}

  void drawSample(const SampleID &id) override {
    scalarGenerator.DrawSample(id);
    vectorFieldGenerator.DrawSample(id);
  }

  double Concentration(double t, const Point &x) const override {
    Point midPoint = Point(0.5, 0.8);
    double r = dist(midPoint, x);
    if (r < 1 / cc)
      return amplitude * pow(cos(cc * Pi * r) + 1.0, 2.0);
    return 0.0;
  }

  VectorField CellConvection(const Cell &c, const Point &x) const override {
    return this->vectorFieldGenerator.EvalSample(c);
  }

  Scalar FaceConvection(const Cell &c, int face,
                        const VectorField &N, const Point &x) const override {
    return this->scalarGenerator.EvalSample(face, c);
  }

  Tensor Diffusion(const Point &x) const override {
    return diffusion * Zero;
  }

  virtual Scalar Reaction(const Point &, double c) const override {
    return reaction * 0.0;
  }

  virtual Scalar DerivativeReaction(const Point &, double c) const override {
    return 0.0;
  }

  string Name() const override { return "StochasticReactionPollutionCosHat2D"; }
};

//class LogisticReaction2D : public IStochasticReactionProblem {
//  double r0 = 30.0;
//  double r1 = 10.0;
//  dt = 0.05;
//  T = 0.5;
//}

IStochasticReactionProblem *CreateStochasticReactionProblem(const std::string &problemName) {
//  if (problemName == "StochasticReactionPollutionCosHat2D")
//    return new StochasticReactionPollutionCosHat2D();
//  if (problemName == "StochasticCosHat2D")
//    return new StochasticCosHat2D();
  return nullptr;
}

std::unique_ptr<IStochasticReactionProblem>
CreateStochasticReactionProblemUnique(const std::string &problemName) {
  return std::unique_ptr<IStochasticReactionProblem>(
    CreateStochasticReactionProblem(problemName)
  );
}

std::shared_ptr<IStochasticReactionProblem>
CreateStochasticReactionProblemShared(const std::string &problemName) {
  return std::shared_ptr<IStochasticReactionProblem>(
    CreateStochasticReactionProblem(problemName)
  );
}