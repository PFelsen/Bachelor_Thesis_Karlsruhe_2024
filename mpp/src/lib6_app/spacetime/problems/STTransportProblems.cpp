#include "STTransportProblems.hpp"
#include "ParameterImage.hpp"
#include "VectorField.hpp"

class LinearTransport : public TransportProblem {

public:
  VectorField q;

  LinearTransport() : TransportProblem(2) {
    q[0] = 1;
    q[1] = 1;
  }

  std::string Name() const override { return "LinearTransport"; }

  VectorField FluxVector(const LevelPair &level, const cell c, const Point &x) const override { return q; }

  double divFluxVector(const LevelPair &level, const cell c, const Point &x) const override { return 0; }

  double ut(const cell c, const Point &x) const override {
    double t = x.t();
    if (t == 0.0) return x[0] + 2 * x[1];
    Point y = x;
    y[0] -= t * q[0];
    y[1] -= t * q[1];
    y[2] = 0.0;
    return ut(c, y);
  }

  double InflowData(const LevelPair &level, const cell c, const Point &x, VectorField &N) const override {
    return Flux(level, c, x, ut(c, x)) * N;
  }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }
};

class ConstantTransport : public TransportProblem {

public:
  VectorField q;

  ConstantTransport() : TransportProblem(2) {
    q[0] = 1;
    q[1] = 1;
  }

  std::string Name() const override { return "ConstantTransport"; }

  VectorField FluxVector(const LevelPair &level, const cell c,const Point &x) const override { return q; }

  double divFluxVector(const LevelPair &level, const cell c,const Point &x) const override { return 0; }

  double ut(const cell c,const Point &x) const override {
    double t = x.t();
    if (t == 0.0) return 1;
    Point y = x;
    y[0] -= t * q[0];
    y[1] -= t * q[1];
    y[2] = 0.0;
    return ut(c, y);
  }

  double InflowData(const LevelPair &level, const cell c,const Point &x, VectorField &N) const override {
    return Flux(level, c, x, ut(c, x)) * N;
  }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }
};

class RiemannTransport : public TransportProblem {

public:
  VectorField q;

  RiemannTransport() : TransportProblem(2) {
    Point Q(1.0, 2.0);
    Config::Get("FluxVector", Q);
    q = Q;
  }

  std::string Name() const override { return "RiemannTransport"; }

  VectorField FluxVector(const LevelPair &level, const cell c, const Point &x) const override { return q; }

  double divFluxVector(const LevelPair &level, const cell c, const Point &x) const override { return 0; }

  double ut(const cell c, const Point &x) const override {
    double t = x.t();
    if (t == 0.0) {
      if (x[0] + 2 * x[1] < 2 / 3.0) return 1;
      return 0;
    }
    Point y = x;
    y[0] -= t * q[0];
    y[1] -= t * q[1];
    y[2] = 0.0;
    return ut(c, y);
  }

  double InflowData(const LevelPair &level, const cell c, const Point &x, VectorField &N) const override {
    return Flux(level, c, x, ut(c, x)) * N;
  }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }
};

class BubbleTransport : public TransportProblem {

public:
  VectorField q;

  BubbleTransport() : TransportProblem(2) {
    Point Q;
    Q[0] = 0.2;
    Q[1] = 0.3;
    Config::Get("FluxVector", Q);
    q[0] = Q[0];
    q[1] = Q[1];
  }

  std::string Name() const override { return "BubbleTransport"; }

  VectorField FluxVector(const LevelPair &level, const cell c, const Point &x) const override {
    VectorField qq(0.1 * x[1] * x[1], 0.0);
    return qq;
  }

  double divFluxVector(const LevelPair &level, const cell c, const Point &x) const override { return 0; }

  double Amplitude(const cell c, const Point &x) const {
    const double d = 0.125;
    double nrm = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (nrm >= d) return 0;
    return exp(nrm / (nrm - d));
  }

  double ut(const cell c, const Point &x) const override {
    double t = x.t();
    Point mid(0.25, 0.25);
    if (t == 0.0) return Amplitude(c, x - mid);
    Point y = x;
    y[0] = x[0] - 0.1 * t * x[1] * x[1];
    y[1] = x[1];
    return Amplitude(c, y - mid);
  }

  double InflowData(const LevelPair &level, const cell c, const Point &x, VectorField &N) const override {
    return 0;
  }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }
};

class RotatingTransport : public TransportProblem {
  const double d = 0.125;
  const Point mid = Point(0.5, 0.5);
public:
  RotatingTransport() : TransportProblem(2) {}

  std::string Name() const override { return "RotatingTransport"; }

  VectorField FluxVector(const LevelPair &level, const cell c, const Point &x) const override {
    return {-2 * Pi * (x[1] - mid[1]), 2 * Pi * (x[0] - mid[0])};
  }

  double divFluxVector(const LevelPair &level, const cell c, const Point &x) const override { return 0; }

  double Amplitude(const cell c, double &r) const {
    if (r >= d) return 0;
    return exp(r / (r - d));
  }

  double ut(const cell c, const Point &x) const override {
    double t = x.t();
    Point Pmid = mid + 0.25 * Point(cos(2 * Pi * t), sin(2 * Pi * t));
    double r = sqrt(pow(Pmid[0] - x[0], 2.) + pow(Pmid[1] - x[1], 2.));
    return Amplitude(c, r);
  }

  double InflowData(const LevelPair &level, const cell c, const Point &x, VectorField &N) const override {
    return 0;
  }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }
};

class Rotating10Transport : public TransportProblem {

public:
  Rotating10Transport() : TransportProblem(2) {}

  std::string Name() const override { return "RotatingTransport"; }

  VectorField FluxVector(const LevelPair &level, const cell c, const Point &x) const override {
    return {-2 * Pi * x[1], 2 * Pi * x[0]};
  }

  double divFluxVector(const LevelPair &level, const cell c, const Point &x) const override { return 0; }

  double Amplitude(const cell c, const Point &x) const {
    const double d = 0.125;
    double nrm = sqrt(x[0] * x[0] + x[1] * x[1]);
    if (nrm >= d) return 0;
    return exp(nrm / (nrm - d));
  }

  double ut(const cell c, const Point &x) const override {
    double t = x.t() / 20;
    Point Pmid = 5.0 * Point(cos(2 * Pi * t), sin(2 * Pi * t));
    double r = sqrt(pow(Pmid[0] - x[0], 2.) + pow(Pmid[1] - x[1], 2.));
    if (r < 2) return exp(r / (r - 2));
    return 0;
  }

  double InflowData(const LevelPair &level, const cell c, const Point &x, VectorField &N) const override {
    return 0;
  }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }
};

TransportProblem *CreateTransportProblemST(const string &name) {
  if (name == "ConstantTransport") return new ConstantTransport();
  if (name == "LinearTransport") return new LinearTransport();
  if (name == "RiemannTransport") return new RiemannTransport();
  if (name == "BubbleTransport") return new BubbleTransport();
  if (name == "RotatingTransport") return new RotatingTransport();

  Exit(name + " Error: no TransportProblem found!");
}

std::unique_ptr<TransportProblem> CreateTransportProblemUniqueST(const string &name) {
  return std::unique_ptr<TransportProblem>(CreateTransportProblemST(name));
}

std::shared_ptr<TransportProblem> CreateTransportProblemSharedST(const string &name) {
  return std::shared_ptr<TransportProblem>(CreateTransportProblemST(name));
}
