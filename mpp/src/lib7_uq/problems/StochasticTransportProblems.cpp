#include "StochasticTransportProblems.hpp"

#if 0
class StochasticSGRiemann1D : public IStochasticTransportProblem {
  SparseGridGenerator sparseGridGen;

  double thickness = 1.0 / 8.0;

  double centerStart = 1.5 * thickness;
public:
  StochasticSGRiemann1D() :
      IProblem("Interval"), sparseGridGen(SparseGridGenerator(GridDomain({{0}, {1}}))) {}

  void drawSample(const SampleID &id) override { sparseGridGen.DrawSample(id); }

  double center(double t) const { return centerStart - t; }

  void InitGenerator(int init) override { sparseGridGen.InitGenerator(meshes, init); }

  int NumOfSamples() const override { return sparseGridGen.GetNumPoints(); }

  double Solution(double t, const Point &x) const override {
    if (abs(center(t) - x[0]) < (thickness * this->sparseGridGen.EvalSample()[0]) / 2.0)
      return 1.0 / (thickness * this->sparseGridGen.EvalSample()[0]);
    else return 0.0;
  }

  static VectorField TransportFlux(const Point &) { return {1.0, 0.0}; }

  VectorField CellFlux(const cell &c, const Point &x) const override { return TransportFlux(x); }

  Scalar FaceNormalFlux(const cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return TransportFlux(x) * N;
  }
};

class StochasticSGGaussHat1D : public IStochasticTransportProblem {
private:
  SparseGridGenerator sparseGridGen;

  double sigma2;

  double factor;
public:
  explicit StochasticSGGaussHat1D(double sigma = 0.05) :
      IProblem("Line"), sigma2(sigma * sigma),
      sparseGridGen(SparseGridGenerator(GridDomain({{0}, {1}}))) {
    T = 0.5;
    factor = 1.0 / sqrt((2.0 * Pi * sigma2));
  }

  void InitGenerator(int init) override { sparseGridGen.InitGenerator(meshes, init); }

  void drawSample(const SampleID &id) override { sparseGridGen.DrawSample(id); }

  int NumOfSamples() const override { return sparseGridGen.GetNumPoints(); }

  static Point mu(double t) { return {0.25, 0.0}; }

  Scalar Solution(double t, const Point &x) const override {
    VectorField diff = (x - mu(t) * this->sparseGridGen.EvalSample()[0]);
    Scalar exponent = -diff * diff / (2 * sigma2);
    return factor * exp(exponent);
  }

  static VectorField TransportFlux(const Point &x) { return {1.0, 0.0}; }

  VectorField CellFlux(const cell &c, const Point &x) const override { return TransportFlux(x); }

  Scalar FaceNormalFlux(const cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return TransportFlux(x) * N;
  }
};

class StochasticSGline1D : public IStochasticTransportProblem {
private:
  SparseGridGenerator sparseGridGen;
public:
  StochasticSGline1D() :
      IProblem("Line"), sparseGridGen(SparseGridGenerator(GridDomain({{0}, {1}}))) {
    T = 0.5;
  }

  void InitGenerator(int init) override { sparseGridGen.InitGenerator(meshes, init); }

  void drawSample(const SampleID &id) override { sparseGridGen.DrawSample(id); }

  int NumOfSamples() const override { return sparseGridGen.GetNumPoints(); }

  Scalar Solution(double t, const Point &x) const override { return x[0] + x[1]; }

  static VectorField TransportFlux(const Point &x) { return {1.0, 0.0}; }

  VectorField CellFlux(const cell &c, const Point &x) const override { return TransportFlux(x); }

  Scalar FaceNormalFlux(const cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return TransportFlux(x) * N;
  }
};

class StochasticSGGaussHat2D : public IStochasticTransportProblem {
private:
  SparseGridGenerator sparseGridGen;

  double factor;

  Tensor sigma;
public:
  StochasticSGGaussHat2D() :
      IProblem("UnitSquare"), sigma(0.005, 0.0, 0.0, 0.005),
      sparseGridGen(SparseGridGenerator(GridDomain({{0, 0}, {0.1, 0.1}}))) {
    factor = 1.0 / sqrt(pow(2.0 * Pi, 2) * sigma.det());
    T = 2.0;
  }

  void InitGenerator(int init) override { sparseGridGen.InitGenerator(meshes, init); }

  void drawSample(const SampleID &id) override { sparseGridGen.DrawSample(id); }

  int NumOfSamples() const override { return sparseGridGen.GetNumPoints(); }

  static Point mu(double t) { return {0.25 + t, 0.25 + t}; }

  Scalar Solution(double t, const Point &x) const override {
    VectorField diff =
        (x - mu(t)
         + Point(this->sparseGridGen.EvalSample()[0], this->sparseGridGen.EvalSample()[1]));
    VectorField temp = Invert(sigma) * diff;
    Scalar exponent = -0.5 * diff * temp;
    return factor * exp(exponent);
  }

  static VectorField TransportFlux(const Point &x) {
    return {cos(Pi * x[0]) / 2, sin(Pi * x[1]) / 2};
  }

  VectorField CellFlux(const cell &c, const Point &x) const override { return TransportFlux(x); }

  Scalar FaceNormalFlux(const cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return TransportFlux(x) * N;
  }
};

class StochasticSGDWave2D : public IStochasticTransportProblem {
private:
  SparseGridGenerator sparseGridGen;
public:
  StochasticSGDWave2D() :
      IProblem("UnitSquare"), sparseGridGen(SparseGridGenerator(GridDomain({{0, 0}, {1, 1}}))) {}

  void InitGenerator(int init) override { sparseGridGen.InitGenerator(meshes, init); }

  void drawSample(const SampleID &id) override { sparseGridGen.DrawSample(id); }

  int NumOfSamples() const override { return sparseGridGen.GetNumPoints(); }

  Scalar Solution(double t, const Point &x) const override {
    if (x[0] == 0) return 1 + this->sparseGridGen.EvalSample()[0];
    else if (x[0] == 1) return 1 + this->sparseGridGen.EvalSample()[1];
    else return 0;
  }

  static VectorField TransportFlux(const Point &x) {
    if (x[0] < 0.5) return {1.0, 0.0};
    else return {-1.0, 0.0};
  }

  VectorField CellFlux(const cell &c, const Point &x) const override { return TransportFlux(x); }

  Scalar FaceNormalFlux(const cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return TransportFlux(x) * N;
  }
};

class StochasticSGCWave2D : public IStochasticTransportProblem {
private:
  SparseGridGenerator sparseGridGen;
public:
  StochasticSGCWave2D() :
      IProblem("UnitSquare"), sparseGridGen(SparseGridGenerator(GridDomain({{0, 0}, {1, 1}}))) {}

  void InitGenerator(int init) override { sparseGridGen.InitGenerator(meshes, init); }

  void drawSample(const SampleID &id) override { sparseGridGen.DrawSample(id); }

  int NumOfSamples() const override { return sparseGridGen.GetNumPoints(); }

  Scalar Solution(double t, const Point &x) const override {
    return exp(-x[0]) * this->sparseGridGen.EvalSample()[0]
           + exp(-x[1]) * this->sparseGridGen.EvalSample()[1];
  }

  static VectorField TransportFlux(const Point &x) { return {x[0] * 0.5, x[1] * 0.5}; }

  VectorField CellFlux(const cell &c, const Point &x) const override { return TransportFlux(x); }

  Scalar FaceNormalFlux(const cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return TransportFlux(x) * N;
  }
};

#endif

class StochasticPollutionRiemann2D : public IStochasticTransportProblem {
  HybridFaceNormalFluxGenerator scalarGenerator;

  HybridCellFluxGenerator vectorFieldGenerator;

  double d = 1.0 / 16.0;
public:
  StochasticPollutionRiemann2D() :
      IProblem("UnitSquare"), scalarGenerator(HybridFaceNormalFluxGenerator()),
      vectorFieldGenerator(HybridCellFluxGenerator(scalarGenerator)) {}

  void drawSample(const SampleID &id) override {
    scalarGenerator.DrawSample(id);
    vectorFieldGenerator.DrawSample(id);
  }

  Scalar Solution(double t, const Cell &c, const Point &x) const override {
    if ((abs(1 - x[1] - 1.5 * d) < d / 2) && (abs(0.5 - x[0]) < 3 * d)) return 1.0;
    return 0.0;
  }

  VectorField Flux(const LevelPair &level, const Cell &c, const Point &x) const override {
    return this->vectorFieldGenerator.EvalSample(c);
  }

  Scalar FaceNormalFlux(const LevelPair &level, const Cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return this->scalarGenerator.EvalSample(face, c);
  }
};

class StochasticPollutionCosHat2D : public IStochasticTransportProblem {
protected:
  HybridFaceNormalFluxGenerator scalarGenerator;

  HybridCellFluxGenerator vectorFieldGenerator;

  double amplitude = 1.0;

  double cc = 6.0;
public:
  explicit StochasticPollutionCosHat2D() :
      IProblem("UnitSquare"), scalarGenerator(HybridFaceNormalFluxGenerator()),
      vectorFieldGenerator(HybridCellFluxGenerator(scalarGenerator)) {}

  void drawSample(const SampleID &id) override {
    scalarGenerator.DrawSample(id);
    vectorFieldGenerator.DrawSample(id);
  }

  double Solution(double t, const Cell &c, const Point &x) const override {
    Point midPoint = Point(0.5, 0.8);
    double r = dist(midPoint, x);
    if (r < 1 / cc) return amplitude * pow(cos(cc * Pi * r) + 1.0, 2.0);
    return 0.0;
  }

  VectorField Flux(const LevelPair &level, const Cell &c, const Point &x) const override {
    return this->vectorFieldGenerator.EvalSample(c);
  }

  Scalar FaceNormalFlux(const LevelPair &level, const Cell &c, int face, const VectorField &N,
                        const Point &x) const override {
    return this->scalarGenerator.EvalSample(face, c);
  }
};

class StochasticPollutionTrippleCosHat2D : public StochasticPollutionCosHat2D {
  RVector sample;
public:
  explicit StochasticPollutionTrippleCosHat2D() :
      IProblem("UnitSquare"), StochasticPollutionCosHat2D() {}

  void drawSample(const SampleID &id) override {
    scalarGenerator.DrawSample(id);
    vectorFieldGenerator.DrawSample(id);
    sample = Random::Uniform(id.commSplit, 3, 0.25, 0.75);
  }

  double Solution(double t, const Cell &c, const Point &x) const override {
    std::vector<Point> midPoints{Point(sample[0], 0.8), Point(sample[1], 0.8),
                                 Point(sample[2], 0.8)};
    double value = 0.0;
    for (auto &point : midPoints) {
      double r = dist(point, x);
      if (r < 1 / cc) value = max(amplitude * pow(cos(cc * Pi * r) + 1.0, 2.0), value);
    }
    return value;
  }
};

IStochasticTransportProblem *CreateStochasticTransportProblem(const string &problemName) {
  IStochasticTransportProblem *problem = nullptr;
#if 0
  if (problemName == "StochasticSGRiemann1D") problem = new StochasticSGRiemann1D();
  if (problemName == "StochasticSGGaussHat1D") problem = new StochasticSGGaussHat1D();
  if (problemName == "StochasticSGline1D") problem = new StochasticSGline1D();
  if (problemName == "StochasticSGGaussHat2D") problem = new StochasticSGGaussHat2D();
  if (problemName == "StochasticSGDWave2D") problem = new StochasticSGDWave2D();
  if (problemName == "StochasticSGCWave2D") problem = new StochasticSGCWave2D();
#endif

  if (problemName == "StochasticPollutionRiemann2D") problem = new StochasticPollutionRiemann2D();
  if (problemName == "StochasticPollutionCosHat2D") problem = new StochasticPollutionCosHat2D();
  if (problemName == "StochasticPollutionTrippleCosHat2D")
    problem = new StochasticPollutionTrippleCosHat2D();

  return problem;
}

std::unique_ptr<IStochasticTransportProblem>
CreateStochasticTransportProblemUnique(const std::string &problemName) {
  return std::unique_ptr<IStochasticTransportProblem>(
      CreateStochasticTransportProblem(problemName));
}

std::shared_ptr<IStochasticTransportProblem>
CreateStochasticTransportProblemShared(const std::string &problemName) {
  return std::shared_ptr<IStochasticTransportProblem>(
      CreateStochasticTransportProblem(problemName));
}