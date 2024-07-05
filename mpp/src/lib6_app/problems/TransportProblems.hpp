#ifndef _DGPROBLEM_H_
#define _DGPROBLEM_H_

#include "IProblem.hpp"
#include "Algebra.hpp"


class ITransportProblem : virtual public IProblem {
protected:
  double cfl = 0.25;

  double t0 = 0.0;

  double T = 1.0;

public:
  ITransportProblem() {
    Config::Get("CFL", cfl);
    Config::Get("t0", t0);
    Config::Get("T", T);
  }

  virtual double GetStepSize(int level) const { return cfl * pow(2, -level); }

  double GetStartTime() const { return t0; }

  double GetEndTime() const { return T; }

  double GetCFL() const { return cfl; }

  virtual bool RHS() const { return false; };

  virtual double Solution(double t, const Cell &c, const Point &x) const {
    THROW("Solution for Problem " + Name() + " not implemented.")
  }

  virtual double Inflow(double t, const Cell &c, const Point &x) const {
    THROW("Solution for Problem " + Name() + " not implemented.")
  }

  virtual VectorField Flux(const LevelPair &level, const Cell &c, const Point &x) const { return zero; }

  virtual double FaceNormalFlux(const LevelPair &level, const Cell &c, int f, const VectorField &N,
                                const Point &x) const {
    return Flux(level, c, x) * N;
  }

  virtual double divFluxVector(const LevelPair &level, const Cell &c, const Point &x) const { return 0; }

  double DivFlux(const LevelPair &level, const Cell &c, const Point &x, double &u, VectorField &Du) const {
    return u * divFluxVector(level, c, x) + Flux(level, c, x) * Du;
  }

  bool HasInflow(const LevelPair &level,
                 const Cell &c,
                 const Point &x,
                 VectorField N) const {
    return Flux(level, c, x) * N < 0;
  }

  virtual double InflowData(const LevelPair &level, const Cell &c, double time, const Point &x, VectorField &N) const {
    return 0.0;
  }

  bool HasOutflow(const LevelPair &level, const Cell &c, const Point &x, VectorField N) const {
    return Flux(level, c, x) * N >= 0;
  }

  virtual IProblem& GetFluxProblem() {
    THROW("FluxProblem not available for " + Name());
    return *this;
  }

  virtual bool HasFluxProblem() const {
    return false;
  }

  virtual bool InternalInflow() const { return false; }
};

class TransportPDEProblem : public ITransportProblem {
private:
  friend class TransportPDEProblemBuilder;

  using FluxFunction = std::function<VectorField(const LevelPair &level, const Cell &c, const Point &x)>;
  using SolutionFunction = std::function<Scalar(double t, const Cell &C, const Point &x)>;
  using FaceNormalFluxFunction = std::function<double(const LevelPair &level, const Cell &c, int f, const VectorField &N, const Point &x)>;
  using NameFunction = std::function<std::string()>;
  using InflowFunction = std::function<Scalar(double t, const Cell &c, const Point &x)>;

  bool rhs = false;
  bool hasExactSolution = false;
  bool internalInflow = false;

  FluxFunction fluxFunction = nullptr;

  SolutionFunction solutionFunction = nullptr;

  FaceNormalFluxFunction faceNormalFluxFunction = nullptr;

  NameFunction nameFunction = nullptr;

  InflowFunction inflowFunction = nullptr;
public:
  explicit TransportPDEProblem(std::shared_ptr<Meshes> meshes) : IProblem(std::move(meshes)) {}

  explicit TransportPDEProblem(std::string meshName) : IProblem(std::move(meshName)) {}

  bool RHS() const override { return rhs; }

  bool InternalInflow() const override { return internalInflow; }

  bool HasExactSolution() const override { return hasExactSolution; }

  Scalar Solution(double t, const Cell &c, const Point &x) const override {
    if (solutionFunction) return solutionFunction(t, c, x);
    return 0.;
  }

  Scalar Inflow(double t, const Cell &c, const Point &x) const override {
    if (inflowFunction) return inflowFunction(t, c, x);
    return 0.0;
  }

  VectorField Flux(const LevelPair &level, const Cell &c, const Point &x) const override {
    if (fluxFunction) return fluxFunction(level, c, x);
    return zero;
  }

  double FaceNormalFlux(const LevelPair &level, const Cell &c, int f, const VectorField &N,
                                const Point &x) const override {
    if (faceNormalFluxFunction) return faceNormalFluxFunction(level, c, f, N, x);
    return fluxFunction(level, c, x) * N;
  }

  std::string Name() const override {
    if (nameFunction) return nameFunction();
    return "TransportPDEProblem";
  }
};

class TransportPDEProblemBuilder {
private:
  TransportPDEProblem *problem;

public:
  explicit TransportPDEProblemBuilder(std::shared_ptr<Meshes> meshes) :
      problem(new TransportPDEProblem(std::move(meshes))) {}

  explicit TransportPDEProblemBuilder(std::string meshName) :
      problem(new TransportPDEProblem(meshName)) {}

  TransportPDEProblemBuilder &
  WithRHS(const bool rhs) {
    problem->rhs = rhs;
    return *this;
  }

  TransportPDEProblemBuilder &WithInternalInflow(const bool internalInflow) {
    problem->internalInflow = internalInflow;
    return *this;
  }

  TransportPDEProblemBuilder &WithInflow(TransportPDEProblem::InflowFunction func) {
    problem->inflowFunction = std::move(func);
    return *this;
  }

  TransportPDEProblemBuilder &
  WithHasExactSolution(const bool hasExactSolution) {
    problem->hasExactSolution = hasExactSolution;
    return *this;
  }

  TransportPDEProblemBuilder &
  WithSolution(TransportPDEProblem::SolutionFunction func) {
    problem->solutionFunction = std::move(func);
    return *this;
  }

  TransportPDEProblemBuilder &
  WithFlux(TransportPDEProblem::FluxFunction func) {
    problem->fluxFunction = std::move(func);
    return *this;
  }

  TransportPDEProblemBuilder &
  WithFaceNormalFlux(TransportPDEProblem::FaceNormalFluxFunction func) {
    problem->faceNormalFluxFunction = std::move(func);
    return *this;
  }

  TransportPDEProblemBuilder &
  WithName(TransportPDEProblem::NameFunction func) {
    problem->nameFunction = std::move(func);
    return *this;
  }

  TransportPDEProblem *Build() {
    return problem;
  }

  static TransportPDEProblemBuilder RiemannTransport1D();
  static TransportPDEProblemBuilder InflowTest();
  static TransportPDEProblemBuilder Inflow();
  static TransportPDEProblemBuilder Hat();
  static TransportPDEProblemBuilder FatHat();
  static TransportPDEProblemBuilder TravelingWave();
  static TransportPDEProblemBuilder SphericalWave2D();
  static TransportPDEProblemBuilder CircleWave2D();
  static TransportPDEProblemBuilder SineWave2D();
  static TransportPDEProblemBuilder CirclePacman();
  static TransportPDEProblemBuilder Riemann1D();
  static TransportPDEProblemBuilder Riemann2D();
  static TransportPDEProblemBuilder GaussHat2D();
  static TransportPDEProblemBuilder CosHat1D();
  static TransportPDEProblemBuilder CosHat2D();
  static TransportPDEProblemBuilder Inflow2D();
};


ITransportProblem *CreateTransportProblem(const string &problemName);

std::unique_ptr<ITransportProblem>
CreateTransportProblemUnique(const std::string &problemName);

std::shared_ptr<ITransportProblem>
CreateTransportProblemShared(const std::string &problemName);

/*#if USE_SPACETIME
#include "STTransportProblems.hpp"
std::unique_ptr<TransportProblem> CreateSTTransportProblemForTutorial(const string &problemName);
#endif*/

#endif
