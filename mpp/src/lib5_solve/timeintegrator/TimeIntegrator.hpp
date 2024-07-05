#ifndef TIMEINTEGRATOR_HPP
#define TIMEINTEGRATOR_HPP

#include "Assemble.hpp"
#include "LinearSolver.hpp"
#include "Newton.hpp"
#include "RMatrix.hpp"

enum TIMEINTEGRATOR {
  EXPLICIT_EULER = 1,
  HEUN = 2,
  RUNGE = 3,
  RUNGE_KUTTA = 4,
  IMPLICIT_EULER = -1,
  IMPLICIT_MIDPOINT = -2,
  CRANK_NICOLSON = -102,
  CRANK_NICOLSON_HARDCODED = -302,
  DI_RK_ORDER2 = -202,
  DI_RK_ORDER3 = -3,
  EXPONENTIAL_EULER = 1001,
  EXPONENTIAL_MIDPOINT = 1002,
  EXPONENTIAL_SHIFT_EULER = -1001,
  EXPONENTIAL_SHIFT_MIDPOINT = -1002,
  EXPONENTIAL = 0
};

bool NotDiverged(Vector &u, bool divergenceCheck);

bool Diverged(Vector &u, bool divergenceCheck);

class TimeIntegrator {
protected:
  int verbose = 1;

  bool divergenceCheck = true;

  std::vector<LinearSolver *> solvers{};

  int stagesSum = 0;
  std::unique_ptr<Matrix> massMatrix = nullptr;

  std::unique_ptr<Matrix> systemMatrix = nullptr;

  std::unique_ptr<Vector> du = nullptr;

  virtual void initializeSpecific(ILinearTimeAssemble *assemble, Vector &u){};
  virtual void destructSpecific(){};
public:
  TimeIntegrator(LinearSolver *solver1) {
    Config::Get("TimeIntegratorVerbose", verbose);
    solvers.push_back(solver1);
  }

  TimeIntegrator(LinearSolver *solver1, LinearSolver *solver2) {
    Config::Get("TimeIntegratorVerbose", verbose);
    solvers.push_back(solver1);
    solvers.push_back(solver2);
  }

  TimeIntegrator(LinearSolver *solver1, LinearSolver *solver2, LinearSolver *solver3,
                 LinearSolver *solver4) {
    Config::Get("TimeIntegratorVerbose", verbose);
    solvers.push_back(solver1);
    solvers.push_back(solver2);
    solvers.push_back(solver3);
    solvers.push_back(solver4);
  }

  virtual ~TimeIntegrator() {
    for (auto &solver : solvers)
      delete solver;
    solvers.clear();
  }

  void SetDivergenceCheck(bool);
  virtual int GetOrder() const = 0;

  virtual const char *Name() const = 0;

  virtual void StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                              const Matrix &systemMatrix, Vector &du) = 0;

  virtual void AssembleRHS(ILinearTimeAssemble *assemble, const Vector &u) = 0;

  void Initialize(ILinearTimeAssemble *assemble, Vector &u);


  bool Method(ILinearTimeAssemble *assemble, Vector &u);

  bool Method(ILinearTimeAssemble *assemble, Vectors &u);

  bool OneTimeStep(ILinearTimeAssemble *assemble, Vector &u);

  void Destruct();

  virtual void PrintInfo() const;
};

class NonLinearTimeIntegrator {
protected:
  int verbose = 1;

  std::unique_ptr<NonLinearSolver> nonLinSolver{};
public:
  NonLinearTimeIntegrator() { Config::Get("TimeIntegratorVerbose", verbose); }

  explicit NonLinearTimeIntegrator(NonLinearSolver *nonLinSolver) :
      nonLinSolver(std::unique_ptr<NonLinearSolver>(nonLinSolver)) {
    Config::Get("TimeIntegratorVerbose", verbose);
  }

  virtual bool Method(INonLinearTimeAssemble &assemble, Vector &) = 0;

  virtual int GetOrder() const = 0;

  virtual const char *Name() const = 0;

  virtual void PrintInfo() const;
};

class TimeIntegratorCreator {
private:
  TIMEINTEGRATOR type;

  std::vector<LinearSolver *> linearSolvers;

  std::vector<NonLinearSolver *> nonLinearSolvers;
public:
  TimeIntegratorCreator() = default;

  TimeIntegratorCreator(TIMEINTEGRATOR type) : type(type) {}

  TimeIntegratorCreator WithConfigTypeEntry() {
    int _type;
    Config::Get("rkorder", _type);
    this->type = static_cast<TIMEINTEGRATOR>(_type);
    return *this;
  }

  TimeIntegratorCreator WithRkOrder(int rkorder) {
    this->type = static_cast<TIMEINTEGRATOR>(rkorder);
    return *this;
  }

  TimeIntegratorCreator WithLinearSolver(LinearSolver *linearSolver) {
    linearSolvers.push_back(std::move(linearSolver));
    return *this;
  }

  TimeIntegratorCreator WithNonLinearSolver(NonLinearSolver *nonLinearSolver) {
    nonLinearSolvers.push_back(std::move(nonLinearSolver));
    return *this;
  }

  TimeIntegrator *Create();

  std::unique_ptr<TimeIntegrator> CreateUnique();

  NonLinearTimeIntegrator *CreateNonLinearTimeIntegrator();

  std::unique_ptr<NonLinearTimeIntegrator> CreateNonLinearTimeIntegratorUnique();
};

#endif // TIMEINTEGRATOR_HPP
