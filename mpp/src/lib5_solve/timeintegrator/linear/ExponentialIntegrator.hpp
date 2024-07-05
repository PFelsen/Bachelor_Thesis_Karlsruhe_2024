#ifndef ARNOLDI_HPP
#define ARNOLDI_HPP

#include <Jacobi.hpp>
#include "TimeIntegrator.hpp"

/*
 * Todo find theory and test appropriately
 */

class ExponentialIntegrator : public TimeIntegrator {
protected:
  bool assembled = false;
  double dtOld = -infty;
  double gamma = 5;
  std::unique_ptr<Matrix> stageMatrix = nullptr;
  bool keepLazy = false;
  int KMax = 25;
  double Keps = 0.00001;
  vector<Vector> rhs;
  std::shared_ptr<LazyVectors> vInt = nullptr;
  std::shared_ptr<LazyVectors> wInt = nullptr;
  double betaEps = Eps;
  double KepsRel = 10e-8;

  static void combineSolution(const RVector &s, LazyVectors &v, Vector &du, int k, double dt);

  void destructSpecific() override {
    dtOld = -infty;
    assembled = false;
  }

  void initializeSpecific(ILinearTimeAssemble *assemble, Vector &u) override {
    rhs = std::vector<Vector>(1, u);
  };

  void assembleRationalMatrix(double dt, const Matrix &MassMatrix, const Matrix &SystemMatrix);
  void arnoldiPolynomial(double dt, const Vector &u, const Matrix &massMatrix,
                         const Matrix &systemMatrix, Vector &du, int function);
  void arnoldiRational(double dt, const Vector &u, const Matrix &massMatrix,
                       const Matrix &systemMatrix, Vector &du, int function);

  void fillMat(const RMatrix &h, RMatrix &H_k);

  void polynomialOneStep(const Matrix &systemMatrix, const Matrix &massMatrix, RMatrix &h,
                         Vector &w, LazyVectors &v, int k);

  void rationalOneStep(const Matrix &systemMatrix, const Matrix &massMatrix, RMatrix &h,
                       LazyVectors &w, LazyVectors &v, int k);
public:
  ExponentialIntegrator(LinearSolver *solver, LinearSolver *solver2) :
      TimeIntegrator(solver, solver2) {
    Config::Get("Kmax", KMax);
    Config::Get("Keps", Keps);
    Config::Get("KeepLazy", keepLazy);
    Config::Get("KepsRel", KepsRel);
    Config::Get("BetaEps", betaEps);
  }

  double getEta(RVector &s_old, RVector &s);

  void evaluateExp(RMatrix &H_k, RVector &s, int function);
};

class ExponentialEuler : public ExponentialIntegrator {
public:
  ExponentialEuler(LinearSolver *solver) : ExponentialIntegrator(solver, nullptr){};

  void AssembleRHS(ILinearTimeAssemble *assemble, const Vector &u) override;
  void StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                      const Matrix &systemMatrix, Vector &du) override;

  int GetOrder() const override { return 1; }

  const char *Name() const override { return "ExponentialEuler"; }
};

class ExponentialMidpoint : public ExponentialIntegrator {
public:
  ExponentialMidpoint(LinearSolver *solver) : ExponentialIntegrator(solver, nullptr){};

  void AssembleRHS(ILinearTimeAssemble *assemble, const Vector &u) override;
  void StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                      const Matrix &systemMatrix, Vector &du) override;

  int GetOrder() const override { return 2; }

  const char *Name() const override { return "ExponentialMidpoint"; }
};

class ExponentialMidpointWithShift : public ExponentialIntegrator {
public:
  ExponentialMidpointWithShift(LinearSolver *solver, LinearSolver *solver2) :
      ExponentialIntegrator(solver, solver2) {
    Config::Get("Gamma", gamma);
  };

  void AssembleRHS(ILinearTimeAssemble *assemble, const Vector &u) override;
  void StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                      const Matrix &systemMatrix, Vector &du) override;

  int GetOrder() const override { return 2; }

  const char *Name() const override { return "Rational ExponentialMidpoint"; }
};

class ExponentialEulerWithShift : public ExponentialIntegrator {
public:
  ExponentialEulerWithShift(LinearSolver *solver, LinearSolver *solver2) :
      ExponentialIntegrator(solver, solver2){};

  void AssembleRHS(ILinearTimeAssemble *assemble, const Vector &u) override {
    rhs = std::vector<Vector>(1, u);
    assemble->RHS(assemble->Time() - assemble->StepSize(), rhs[0]);
  }

  void StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                      const Matrix &systemMatrix, Vector &du) override;

  int GetOrder() const override { return 1; }

  const char *Name() const override { return "ExponentialEuler"; }
};

#endif // ARNOLDI_HPP
