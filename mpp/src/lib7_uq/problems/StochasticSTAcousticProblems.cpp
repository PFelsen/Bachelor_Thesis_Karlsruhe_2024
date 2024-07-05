#include "StochasticSTAcousticProblems.hpp"


class SinCos2D : public IStochasticSTAcousticProblem {
protected:
  void drawSample(const SampleID &id) override {};

public:
  SinCos2D() : IProblem("STSquare"), IStochasticSTAcousticProblem(2, 0) {}

  bool HasExactSolution() const override { return true; }

  bool HasRHS() const override { return false; }

  bool hasInitialData() const override { return true; }

  double ut(double t, const Point &xx, const Cell &c, COMPONENT comp) const override {
    double x = xx[0];
    double y = xx[1];
    switch(comp) {
      case COMPONENT::P0:
        return -sin(M_PI * x) * sin(3 * M_PI * y / 2) * sin(sqrt(13) * M_PI * t / 2);
      case COMPONENT::V_X:
        return 2 * sqrt(13) / 13 * sin(3.0 / 2 * M_PI * y) * cos(M_PI * x) *
          cos(sqrt(13) * M_PI * t / 2);
      case COMPONENT::V_Y:
        return 3 * sqrt(13) / 13 * cos(3.0 / 2 * M_PI * y) * sin(M_PI * x) *
          cos(sqrt(13) * M_PI * t / 2);
      default: Exit("No.");
    }
  };
};

class SphericalWave2D : public IStochasticSTAcousticProblem {
protected:

  void drawSample(const SampleID &id) override {};

  Point midPoint;

  double width;

public:
  SphericalWave2D() :
    IProblem("STSquareCenteredCoarseGeometry"),
    IStochasticSTAcousticProblem(2, 0),
    midPoint(0.5, 0.0, 0.0), width(0.125) {
    Config::Get("MidPoint", midPoint);
    Config::Get("SourceWidth", width);
  }

  double A(const Point &x) const {
    double r = dist(midPoint, x);
    if (r < width) return pow(cos(r * 0.5 * M_PI / width), 6.0);
    else return 0.0;
  }

  double u0(const Point &x, COMPONENT comp) const {
    const double tmp = A(x);
    if (comp == COMPONENT::V_X) return (x[0] - midPoint[0]) * tmp;
    if (comp == COMPONENT::V_Y) return (x[1] - midPoint[1]) * tmp;
    return 0.0;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (x.t() > 0) return 0.0;
    else return u0(x, comp);
  }

  bool hasInitialData() const override { return true; }
};

class StochasticSphericalWave2D : public SphericalWave2D {
  CirculantEmbedding2D gen;

public:
  StochasticSphericalWave2D() : SphericalWave2D(), gen(CirculantEmbedding2D()) {}

  double Rho(const Cell &c, const Point &x) const override {
    return 1.0;
//    return this->gen.EvalSample(x);
  }
};

class Benchmark : public IStochasticSTAcousticProblem {
protected:
  double rho1 = 0.5;

  double rho2 = 2.0;

  double center = -1;

  double width = 1.0;

  void drawSample(const SampleID &id) override {}

public:
  Benchmark() : IProblem("tube_squares"), IStochasticSTAcousticProblem(2, 0) {
    Config::Get("width", width);
    Config::Get("center", center);
  }

  double Rho(const Cell &c, const Point &x) const override {
    if (x[0] < 0) return 1.0;
    else if (x[0] < 1) return rho1;
    else return rho2;
  }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    if (i != 0) { THROW("Kappa_i called with i != 0")}
    return 1 / Rho(c, x);
  }

  double A0(double s) const {
    if (abs(s) < width) return pow(cos(s * 0.5 * M_PI / width), 6.0);
    return 0;
  }

  bool hasInitialData() const override { return true; }

  double A(double s) const { return A0(s - center); }

  double DA0(double s) const {
    if (abs(s) < width)
      return 3.0 * M_PI / width * pow(cos(s * 0.5 * M_PI / width), 5.0)
        * sin(s * 0.5 * M_PI / width);
    else return 0.0;
  }

  double DA(double s) const { return DA0(s - center); }

  double u0(const Point &x, COMPONENT comp) const {
    const double tmp = A(x[0]);
    switch(comp) {
      case COMPONENT::P0:
        return tmp;
      case COMPONENT::V_X:
        return -tmp;
      default:
        return 0.0;
    }
  }

  double ut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    if (x.t() == 0) return u0(x, comp);

    if (x[0] < 0) {
      x[0] -= t;
      return u0(x, comp);
    } else if (x[0] < 1) {
      x[0] = rho1 * x[0] - t;
      return u0(x, comp);
    } else {
      x[0] = rho1 + rho2 * (x[0] - 1) - t;
      return u0(x, comp);
    }
  }

  bool HasExactSolution() const override { return true; }

  double Mdtut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    x[0] -= t;
    const double tmp = DA(x[0]);
    switch(comp) {
      case COMPONENT::P0:
        return -tmp;
      case COMPONENT::V_X:
        return tmp;
      default:
        return 0.0;
    }
  }

  double Aut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    x[0] -= t;
    const double tmp = DA(x[0]);
    switch(comp) {
      case COMPONENT::P0:
        return -tmp;
      case COMPONENT::V_X:
        return tmp;
      default:
        return 0.0;
    }
  }
};

class Benchmark1D : public IStochasticSTAcousticProblem {
protected:
  CirculantEmbedding1D circEmb;

  double center = 0.2;

  double width = 0.3;

  void drawSample(const SampleID &id) override {}

public:
  Benchmark1D() : IProblem("SpaceTimeUnitInterval"), IStochasticSTAcousticProblem(1, 0) {
    Config::Get("width", width);
    Config::Get("center", center);
  }

  double Rho(const Cell &c, const Point &x) const override {
    return 1.0;
//    if (x[0] < 0.25) return 0.5;
//    if (x[0] < 0.75) return 0.5 + this->circEmb.EvalSample(x);
//    return 0.5;
  }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    if (i != 0) { THROW("Kappa_i called with i != 0")}
    return 1 / Rho(c, x);
  }

  double A0(double s) const {
    if (abs(s) < width) return pow(cos(s * 0.5 * M_PI / width), 2.0);
    return 0;
  }

  bool hasInitialData() const override { return true; }

  double A(double s) const { return A0(s - center); }

  double u0(const Point &x, COMPONENT comp) const {
    const double tmp = A(x[0]);
    switch(comp) {
      case COMPONENT::P0:
        return tmp;
      case COMPONENT::V_X:
        return -tmp;
      default:
        return 0.0;
    }
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (t == 0) return u0(x, comp);
    return 0.0;
  }

  bool HasExactSolution() const override { return true; }

  double Mdtut(double t, const Point &xx, const Cell &c, COMPONENT comp) const override { return 0.0; }

  double Aut(double t, const Point &xx, const Cell &c, COMPONENT comp) const override { return 0.0; }
};

class StochasticBenchmark : public Benchmark {
private:
  double factor = exp(-1.0 / 8.0);

  double a = -sqrt(12) / 2.0;

  double b = sqrt(12) / 2.0;

  double a_min = 1.0;

  RVector sample;

  void drawSample(const SampleID &id) override { sample = Random::Uniform(id.commSplit, 6, a, b); }

public:
  StochasticBenchmark() :
    Benchmark() {}

  double Rho(const Cell &c, const Point &x) const override {
    if (0 < x[0] && x[0] < 1)
      return (a_min +
        exp(factor * (cos(2 * M_PI * x[0]) * sample[0]
          + sin(2 * M_PI * x[0]) * sample[1]
          + cos(2 * M_PI * x[1]) * sample[2]
          + sin(2 * M_PI * x[1]) * sample[3]
          + cos(M_PI * x.t()) * sample[4]
          + sin(M_PI * x.t()) * sample[5]
        )));
    else return 1.0;
  }
};

class CosHatAndRicker2D : public IStochasticSTAcousticProblem {
protected:
  void drawSample(const SampleID &id) override {};

  Point shotLocation{0.5, 0.5, 0.0};

  double factor = 10.0;

public:

  CosHatAndRicker2D() : IProblem("STSquare"), IStochasticSTAcousticProblem(2, 0) {}

  double F(double t, const Cell &c, const Point &x, COMPONENT comp) const override {
    if (comp != COMPONENT::P0) return 0.0;
    return factor * Ricker(x.t() - shotLocation.t(), 0.1) * CosHat(dist(shotLocation, x), 0.1);
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override { return 0.0; }

  bool HasRHS() const override { return true; }
};

class GaussHatAndRicker2D : public IStochasticSTAcousticProblem {
protected:
  void drawSample(const SampleID &id) override {};

  double factor = 10.0;

  Point shotLocation;

public:
  GaussHatAndRicker2D(Point shotLocation = {0.5, 0.75, 0.0}) :
    IProblem("STSquare"), shotLocation(shotLocation),
    IStochasticSTAcousticProblem(2, 0) {}

  double F(double t, const Cell &c, const Point &x, COMPONENT comp) const override {
    if (comp != COMPONENT::P0) return 0.0;
    return factor * Ricker(x.t() - shotLocation.t(), 0.1) * GaussHat(dist(shotLocation, x), 0.1);
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override { return 0.0; }

  int BndID(const Point &z) const override { return 1; }

  bool HasRHS() const override { return true; }
};

class StochasticGaussHatAndRicker2D : public GaussHatAndRicker2D {
  CirculantEmbedding2D circEmb;

  double rhoMin = 0.25;

public:
  StochasticGaussHatAndRicker2D() :
    IProblem("STSquare"), circEmb(CirculantEmbedding2D()),
    GaussHatAndRicker2D(Point{0.5, 0.75, 0.0}) {
    Config::Get("RhoMin", rhoMin);
  }

  void drawSample(const SampleID &id) override {}

  double Rho(const Cell &c, const Point &z) const override {
    return 1.0;
//    return max(this->circEmb.EvalSample(z), rhoMin);
  }
};

IStochasticSTAcousticProblem *
CreateAcousticStochasticProblem(const string &problemName) {
  IStochasticSTAcousticProblem *problem = nullptr;
  if (problemName == "Benchmark1D")
    problem = new Benchmark1D();

  if (problemName == "SinCos2D")
    problem = new SinCos2D();

  if (problemName == "SphericalWave2D")
    problem = new SphericalWave2D();
  if (problemName == "StochasticSphericalWave2D")
    problem = new StochasticSphericalWave2D();

  if (problemName == "Benchmark")
    problem = new Benchmark();
  if (problemName == "StochasticBenchmark")
    problem = new StochasticBenchmark();

  if (problemName == "CosHatAndRicker2D")
    problem = new CosHatAndRicker2D();
  if (problemName == "GaussHatAndRicker2D")
    problem = new GaussHatAndRicker2D();
  if (problemName == "StochasticGaussHatAndRicker2D")
    problem = new StochasticGaussHatAndRicker2D();
  return problem;
}

std::unique_ptr<IStochasticSTAcousticProblem>
CreateViscoAcousticStochasticProblemUnique(const std::string &problemName) {
  return std::unique_ptr<IStochasticSTAcousticProblem>(
    CreateAcousticStochasticProblem(problemName)
  );
}

std::shared_ptr<IStochasticSTAcousticProblem>
CreateViscoAcousticStochasticProblemShared(const std::string &problemName) {
  return std::shared_ptr<IStochasticSTAcousticProblem>(
    CreateAcousticStochasticProblem(problemName)
  );
}
