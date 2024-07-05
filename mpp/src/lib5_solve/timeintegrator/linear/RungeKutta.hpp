#ifndef RUNGEKUTTA_HPP
#define RUNGEKUTTA_HPP

#include <memory>
#include <utility>

#include "TimeIntegrator.hpp"

struct ButcherTableau {
  RVector c;
  RVector b;
  RMatrix A;
};

class CrankNicolsonHardcoded : public TimeIntegrator {
  std::unique_ptr<Vector> RHS;
  bool rhsInitialized = false;
  bool assembled = false;
  double dtOld = -infty;
  std::unique_ptr<Matrix> stageMatrix = nullptr;
public:
  CrankNicolsonHardcoded(LinearSolver *solver1, LinearSolver *solver2, LinearSolver *solver3,
                         LinearSolver *solver4) :
      TimeIntegrator(solver1, solver2, solver3, solver4) {}

  int GetOrder() const override { return 2; };

  const char *Name() const override { return "CrankNicolsonHardcoded"; };

  void StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                      const Matrix &systemMatrix, Vector &du) override;

  void AssembleRHS(ILinearTimeAssemble *assemble, const Vector &u) override;
};

class RungeKutta : public TimeIntegrator {
protected:
  ButcherTableau tableau;

  const RVector &c;

  const RVector &b;

  const RMatrix &A;

  std::vector<Vector> rhs;

  bool rhsInitialized = false;
public:
  RungeKutta(LinearSolver *solver, ButcherTableau tableau) :
      TimeIntegrator(solver), tableau(std::move(tableau)), c(this->tableau.c), b(this->tableau.b),
      A(this->tableau.A) {}

  RungeKutta(LinearSolver *solver1, LinearSolver *solver2, LinearSolver *solver3,
             LinearSolver *solver4, ButcherTableau tableau) :
      TimeIntegrator(solver1, solver2, solver3, solver4), tableau(std::move(tableau)),
      c(this->tableau.c), b(this->tableau.b), A(this->tableau.A) {}

  void AssembleRHS(ILinearTimeAssemble *assemble, const Vector &u) override {
    if (!rhsInitialized) {
      rhs = std::vector<Vector>(c.size(), u);
      rhsInitialized = false;
    }
    for (int i = 0; i < c.size(); ++i) {
      assemble->RHS(assemble->Time() - (1 - c[i]) * assemble->StepSize(), rhs[i]);
    }
  }

  int GetOrder() const override;

  void CheckButcherTableau() const;

  ButcherTableau GetButcherTableau() { return tableau; }
};

/*
 * Explicit methods
 */

class ExplicitRungeKutta : public RungeKutta {
protected:
  bool mInvAssembled = false;
public:
  ExplicitRungeKutta(LinearSolver *solver, ButcherTableau tableau) :
      RungeKutta(solver, std::move(tableau)) {}

  void StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                      const Matrix &systemMatrix, Vector &du) override;

  void initializeSpecific(ILinearTimeAssemble *assemble, Vector &u) override {
    (*solvers.at(0))(*massMatrix, 1);
  }
};

class ExplicitEuler : public ExplicitRungeKutta {
public:
  explicit ExplicitEuler(LinearSolver *solver) :
      ExplicitRungeKutta(solver, ButcherTableau{{0.0}, {1.0}, {{0.0}}}) {}

  const char *Name() const override { return "ExplicitEuler"; }
};

class Heun : public ExplicitRungeKutta {
public:
  explicit Heun(LinearSolver *solver) :
      ExplicitRungeKutta(solver, ButcherTableau{{0.0, 1.0}, {0.5, 0.5}, {{0.0, 0.0}, {1.0, 0.0}}}) {
  }

  const char *Name() const override { return "Heun"; }
};

class Runge : public ExplicitRungeKutta {
public:
  explicit Runge(LinearSolver *solver) :
      ExplicitRungeKutta(solver,
                         ButcherTableau{{0.0, 0.5, 1.0},
                                        {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
                                        {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, {-1.0, 2.0, 0.0}}}) {}

  const char *Name() const override { return "Runge"; }
};

class RungeKutta4 : public ExplicitRungeKutta {
public:
  explicit RungeKutta4(LinearSolver *solver) :
      ExplicitRungeKutta(solver, ButcherTableau{{0.0, 0.5, 0.5, 1.0},
                                                {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0},
                                                {{0.0, 0.0, 0.0, 0.0},
                                                 {0.5, 0.0, 0.0, 0.0},
                                                 {0.0, 0.5, 0.0, 0.0},
                                                 {0.0, 0.0, 1.0, 0.0}}}) {}

  const char *Name() const override { return "RungeKutta4"; }
};

/*
 * Diagonal implicit methods
 */

class DiagonalImplicitRungeKutta : public RungeKutta {
  bool assembled = false;
  double dtOld = -infty;
  std::vector<std::shared_ptr<Matrix>> stageMatrix;
protected:
  void initializeStageMatrices(double dt);

  std::shared_ptr<Matrix> getStageMatrix(double dt, int stage);
public:
  DiagonalImplicitRungeKutta(LinearSolver *solver1, LinearSolver *solver2, LinearSolver *solver3,
                             LinearSolver *solver4, ButcherTableau tableau) :
      RungeKutta(solver1, solver2, solver3, solver4, std::move(tableau)) {}

  void destructSpecific() override {
    dtOld = -infty;
    assembled = false;
    for (auto M : stageMatrix) {
      M.reset();
    }
    stageMatrix.clear();
  }

  void StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                      const Matrix &systemMatrix, Vector &du) override;

  void initializeSpecific(ILinearTimeAssemble *assemble, Vector &u) override {}
};

class GenericImplicitEuler : public DiagonalImplicitRungeKutta {
public:
  GenericImplicitEuler(LinearSolver *solver1, LinearSolver *solver2, LinearSolver *solver3,
                       LinearSolver *solver4) :
      DiagonalImplicitRungeKutta(solver1, solver2, solver3, solver4,
                                 ButcherTableau{{1.0}, {1.0}, {{1.0}}}) {}

  const char *Name() const override { return "ImplicitEuler"; }
};

class ImplicitMidpointRule : public DiagonalImplicitRungeKutta {
public:
  ImplicitMidpointRule(LinearSolver *solver1, LinearSolver *solver2, LinearSolver *solver3,
                       LinearSolver *solver4) :
      DiagonalImplicitRungeKutta(solver1, solver2, solver3, solver4,
                                 ButcherTableau{{0.5}, {1.0}, {{0.5}}}) {}

  const char *Name() const override { return "ImplicitMidpointRule"; }
};

class CrankNicolson : public DiagonalImplicitRungeKutta {
public:
  CrankNicolson(LinearSolver *solver1, LinearSolver *solver2, LinearSolver *solver3,
                LinearSolver *solver4) :
      DiagonalImplicitRungeKutta(solver1, solver2, solver3, solver4,
                                 ButcherTableau{{0.0, 1.0}, {0.5, 0.5}, {{0.0, 0.0}, {0.5, 0.5}}}) {
  }

  const char *Name() const override { return "CrankNicolson"; }
};

class DiagonalImplicitRungeKutta2 : public DiagonalImplicitRungeKutta {
public:
  DiagonalImplicitRungeKutta2(LinearSolver *solver1, LinearSolver *solver2, LinearSolver *solver3,
                              LinearSolver *solver4) :
      DiagonalImplicitRungeKutta(solver1, solver2, solver3, solver4,
                                 ButcherTableau{{0.25, 0.75},
                                                {0.5, 0.5},
                                                {{0.25, 0.0}, {0.5, 0.25}}}) {}

  const char *Name() const override { return "DiagonalImplicitRungeKutta2"; }
};

class DiagonalImplicitRungeKutta3 : public DiagonalImplicitRungeKutta {
public:
  DiagonalImplicitRungeKutta3(LinearSolver *solver1, LinearSolver *solver2, LinearSolver *solver3,
                              LinearSolver *solver4) :
      DiagonalImplicitRungeKutta(solver1, solver2, solver3, solver4,
                                 ButcherTableau{{0.5 + sqrt(3) / 6.0, 0.5 - sqrt(3) / 6.0},
                                                {0.5, 0.5},
                                                {{0.5 + sqrt(3) / 6.0, 0.0},
                                                 {-sqrt(3.0) / 3.0, 0.5 + sqrt(3) / 6.0, 0.0}}}) {}

  const char *Name() const override { return "DiagonalImplicitRungeKutta3"; }
};

#endif // DIRK_HPP
