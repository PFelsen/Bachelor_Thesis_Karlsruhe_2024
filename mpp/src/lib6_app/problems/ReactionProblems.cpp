#include "ReactionProblems.hpp"
#include "HybridMain.hpp"


Scalar StartConcentration(const Point &x) {
  double thickness = 1.0 / 8.0;
  double centerStart = 1 - 1.5 * thickness;
  if ((abs(centerStart - x[1]) < thickness / 2) &&
    (abs(0.5 - x[0]) < 3 * thickness))
    return 1;
  else return 0.0;
}

class ExponentialReaction2D : public IReactionProblem {
public:
  ExponentialReaction2D() : IProblem("UnitSquare") {};

  Scalar Load(const Point &x) const override { return 0.0; }

  VectorField Flux(const Point &) const override { return zero; }

  Tensor Diffusion(const Point &x) const override { return diffusion * One; }

  VectorField Convection(const Point &x) const override { return {0.0, -1.0}; }

  Scalar Reaction(const Point &, double c) const override { return reaction * c; }

  Scalar Concentration(double t, const Point &x) const override {
    return StartConcentration(x);
  }

  Scalar DerivativeReaction(const Point &, double c) const override {
    return reaction;
  }

  std::string Name() const override { return "ExponentialReaction2D"; }
};

class ExponentialReactionSquare500 : public ExponentialReaction2D {
public:
  ExponentialReactionSquare500() : IProblem("Square500") {};

  std::string Name() const override { return "ExponentialReactionSquare500"; }
};

class LogisticReaction2D : public IReactionProblem {
private:
  double r0 = 30.0;

  double r1 = 10.0;

public:
  LogisticReaction2D() : IProblem("UnitSquare") {
    Config::Get("Reaction0", r0);
    Config::Get("Reaction1", r1);
    dt = 0.05;
    T = 0.5;
  }

  Scalar Load(const Point &x) const override { return 0.0; }

  VectorField Flux(const Point &) const override { return zero; }

  Tensor Diffusion(const Point &x) const override { return diffusion * One; }

  VectorField Convection(const Point &x) const override { return {0.0, -1.0}; }

  Scalar Reaction(const Point &x, double c) const override {
    return r0 * c - r1 * c * c;
  }

  Scalar Concentration(double t, const Point &x) const override {
    return StartConcentration(x);
  }

  Scalar DerivativeReaction(const Point &x, double c) const override {
    return r0 - 2 * r1 * c;
  }

  std::string Name() const override { return "LogisticReaction2D"; }
};

class LogisticReactionSquare500 : public ExponentialReaction2D {
public:
  LogisticReactionSquare500() : IProblem("Square500") {};

  std::string Name() const override { return "LogisticReactionSquare500"; }
};

class ConvectionDiffusionCosHat2D : public IReactionProblem {
private:
  double amplitude = 1.0;

  double frequency = 8.0;

public:
  ConvectionDiffusionCosHat2D() : IProblem("UnitSquare") {}

  Scalar Load(const Point &x) const override { return 0.0; }

  VectorField Flux(const Point &) const override { return zero; }

  Tensor Diffusion(const Point &x) const override { return Zero; }

  static VectorField CircleVectorField(const Point &x) {
    return 2 * Pi * VectorField(0.5 - x[1], x[0] - 0.5, 0.0);
  }

  VectorField Convection(const Point &x) const override {
    return CircleVectorField(x);
  };

  Scalar Reaction(const Point &x, double c) const override { return 0.0; }

  double Concentration(double t, const Point &x) const override {
    Point midPoint = 0.25 * (Point(cos(2.0 * Pi * t), sin(2.0 * Pi * t)))
      + Point(0.5, 0.5);
    double r = dist(midPoint, x);
    if (r < (Pi / 2.0) / frequency)
      return amplitude * pow(cos(frequency * r), 2.0);
    return 0.0;
  }

  Scalar DerivativeReaction(const Point &, double c) const override { return 0.0; }

  Scalar FaceConvection(const Cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return CircleVectorField(x) * N;
  };

  string Name() const override { return "ConvectionDiffusionCosHat2D"; }
};

class PollutionExponentialReaction2D : public ExponentialReaction2D {
private:
  HybridMain hybridMain;

public:
  PollutionExponentialReaction2D() : IProblem("UnitSquare"), hybridMain("Laplace2D") {}

  VectorField CellConvection(const Cell &c, const Point &x) const override {
    return convection * hybridMain.EvaluateCellFlux(c);
  }

  double FaceConvection(const Cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return convection * hybridMain.EvaluateNormalFlux(c, face);
  }

  std::string Name() const override { return "HybridExponentialReaction2D"; }
};

class PollutionExponentialReactionSquare500 : public ExponentialReactionSquare500 {
private:
  HybridMain hybridMain;

public:
  PollutionExponentialReactionSquare500() :
    IProblem("Square500"), hybridMain("LaplaceSquare500") {}

  VectorField CellConvection(const Cell &c, const Point &x) const override {
    return convection * hybridMain.EvaluateCellFlux(c);
  }

  double FaceConvection(const Cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return convection * hybridMain.EvaluateNormalFlux(c, face);
  }

  std::string Name() const override { return "PollutionExponentialReactionSquare500"; }
};

class PollutionLogisticReaction2D : public LogisticReaction2D {
private:
  HybridMain hybridMain;

public:
  PollutionLogisticReaction2D() : IProblem("UnitSquare"), hybridMain("Laplace2D") {}

  VectorField CellConvection(const Cell &c, const Point &x) const override {
    return convection * hybridMain.EvaluateCellFlux(c);
  }

  double FaceConvection(const Cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return convection * hybridMain.EvaluateNormalFlux(c, face);
  }

  std::string Name() const override { return "HybridLogisticReaction2D"; }
};

class PollutionLogisticReactionSquare500 : public LogisticReactionSquare500 {
private:
  HybridMain hybridMain;

public:
  PollutionLogisticReactionSquare500() :
    IProblem("Square500"), hybridMain("LaplaceSquare500") {}

  VectorField CellConvection(const Cell &c, const Point &x) const override {
    return convection * hybridMain.EvaluateCellFlux(c);
  }

  double FaceConvection(const Cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return convection * hybridMain.EvaluateNormalFlux(c, face);
  }

  std::string Name() const override { return "PollutionLogisticReactionSquare500"; }
};

IReactionProblem *CreateReactionProblem(const std::string &problemName) {
  if (problemName == "ConvectionDiffusionCosHat2D")
    return new ConvectionDiffusionCosHat2D();

  if (problemName == "ExponentialReaction2D")
    return new ExponentialReaction2D();
  if (problemName == "ExponentialReactionSquare500")
    return new ExponentialReactionSquare500();

  if (problemName == "LogisticReaction2D")
    return new LogisticReaction2D();
  if (problemName == "LogisticReactionSquare500")
    return new LogisticReactionSquare500();

  if (problemName == "PollutionExponentialReaction2D")
    return new PollutionExponentialReaction2D();
  if (problemName == "PollutionExponentialReactionSquare500")
    return new PollutionExponentialReactionSquare500();

  if (problemName == "PollutionLogisticReaction2D")
    return new PollutionLogisticReaction2D();
  if (problemName == "PollutionLogisticReactionSquare500")
    return new PollutionLogisticReactionSquare500();

  Exit(problemName + " not found")
}

std::unique_ptr<IReactionProblem>
CreateReactionProblemUnique(const std::string &problemName) {
  return std::unique_ptr<IReactionProblem>(
    CreateReactionProblem(problemName)
  );
}

std::shared_ptr<IReactionProblem>
CreateReactionProblemShared(const std::string &problemName) {
  return std::shared_ptr<IReactionProblem>(
    CreateReactionProblem(problemName)
  );
}