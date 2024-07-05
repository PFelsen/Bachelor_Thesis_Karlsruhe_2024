#ifndef IELLIPTICASSEMBLE_HPP
#define IELLIPTICASSEMBLE_HPP

#include "EllipticProblems.hpp"
#include "Assemble.hpp"

typedef std::pair<double, double> FluxPair;

class IEllipticAssemble : public IAssemble {
protected:
  const IEllipticProblem &problem;

public:
  explicit IEllipticAssemble(const IEllipticProblem &problem) : problem(problem) {};

  virtual const IDiscretization &GetDisc() const = 0;

  virtual std::shared_ptr<const IDiscretization> GetSharedDisc() const = 0;

  const IEllipticProblem &GetProblem() const noexcept { return problem; };

  const char *Name() const override = 0;

  void Initialize(Vector &u) const override = 0;

  double Energy(const Vector &u) const override = 0;

  void Residual(const cell &c, const Vector &u, Vector &r) const override = 0;

  void Jacobi(const cell &c, const Vector &u, Matrix &A) const override = 0;

  virtual double EnergyError(const Vector &u) const = 0;

  virtual double L2(const Vector &u) const = 0;

  virtual double H1(const Vector &u) const = 0;

  virtual double L2Error(const Vector &u) const = 0;

  virtual double L2CellAvgError(const Vector &u) const = 0;

  virtual double MaxError(const Vector &u) const = 0;

  virtual double FluxError(const Vector &u) const = 0;

  virtual double FaceError(const Vector &u) const = 0;

  virtual FluxPair InflowOutflow(const Vector &u) const = 0;

  virtual FluxPair PrescribedInflowOutflow(const Vector &u) const = 0;

  virtual FluxPair OutflowLeftRight(const Vector &u) const = 0;

  virtual double GoalFunctional(const Vector &u) const = 0;

  virtual void SetExactSolution(Vector &u_ex) const = 0;

  virtual void SetPressure(const Vector &u, Vector &p) const { p = u; };

  virtual void SetFlux(const Vector &u, Vector &flux) = 0;

  virtual void SetNormalFlux(const Vector &u, Vector &flux) {
    Exit("Not implemented for " + string(typeid(*this).name()))
  }

  virtual double DualPrimalError(const Vector &u) const {
    Exit("Not implemented for " + string(typeid(*this).name()))
  }

  void PrintInfo() const {
    mout.PrintInfo("Assemble", verbose,
      PrintInfoEntry("Assemble Name", Name()),
      PrintInfoEntry("Problem", problem.Name()));
  }

};

struct PDESolverConfig;

IEllipticAssemble *
CreateEllipticAssemble(const IEllipticProblem &problem, const PDESolverConfig &conf);

std::unique_ptr<IEllipticAssemble>
CreateEllipticAssembleUnique(const IEllipticProblem &problem, const PDESolverConfig &conf);

std::shared_ptr<IEllipticAssemble>
CreateEllipticAssembleShared(const IEllipticProblem &problem, const PDESolverConfig &conf);

#endif //IELLIPTICASSEMBLE_HPP
