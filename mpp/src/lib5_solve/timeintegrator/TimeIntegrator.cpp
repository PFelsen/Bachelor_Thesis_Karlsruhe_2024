#include "TimeIntegrator.hpp"
#include "ExponentialIntegrator.hpp"
#include "ImplicitEuler.hpp"
#include "RungeKutta.hpp"

bool NotDiverged(Vector &u, bool divergenceCheck) {
  if (u.norm() <= 1e10 || !divergenceCheck) {
    return true;
  } else {
    u = 0.0;
    Warning("Time integrator diverged") return false;
  }
}

bool Diverged(Vector &u, bool divergenceCheck) { return !NotDiverged(u, divergenceCheck); }

void TimeIntegrator::SetDivergenceCheck(bool diverged) { divergenceCheck = diverged; }

bool TimeIntegrator::Method(ILinearTimeAssemble *assemble, Vector &u) {
  mout.StartBlock("TI");
  vout(1) << "Starting " << Name() << " (Order=" << GetOrder() << ")" << endl;
  Initialize(assemble, u);
  while (!assemble->IsFinished()) {
    if (!OneTimeStep(assemble, u)) return false;
  }
  Destruct();
  mout.EndBlock(verbose == 0);
  return true;
}

void TimeIntegrator::PrintInfo() const {
  mout.PrintInfo("Time Integrator", verbose, PrintInfoEntry("Name", Name(), 1),
                 PrintInfoEntry("Order", GetOrder(), 1));
}

void TimeIntegrator::Destruct() {
  massMatrix = nullptr;
  systemMatrix = nullptr;
  du = nullptr;
  destructSpecific();
}

bool TimeIntegrator::OneTimeStep(ILinearTimeAssemble *assemble, Vector &u) {
  assemble->NextTimeStep();
  double dt = assemble->StepSize();
  AssembleRHS(assemble, u);
  StageFunctions(dt, u, *massMatrix, *systemMatrix, *du);
  u += *du;
  if (NotDiverged(u, divergenceCheck)) {
    assemble->FinishTimeStep(u);
    return true;
  } else {
    mout.EndBlock(verbose == 0);
    return false;
  }
}

void TimeIntegrator::Initialize(ILinearTimeAssemble *assemble, Vector &u) {
  massMatrix = std::make_unique<Matrix>(u);
  systemMatrix = std::make_unique<Matrix>(u);
  du = std::make_unique<Vector>(0.0, u);

  assemble->Initialize(u);
  assemble->SystemMatrix(*systemMatrix);
  assemble->MassMatrix(*massMatrix);
  assemble->SetInitialValue(u);
  assemble->FinishTimeStep(u);
  initializeSpecific(assemble, u);
  u.Accumulate();
}

bool TimeIntegrator::Method(ILinearTimeAssemble *assemble, Vectors &u) {
  mout.StartBlock("TI");
  vout(1) << "Starting " << Name() << " (Order=" << GetOrder() << ")" << endl;
  Initialize(assemble, u[0]);
  int step = 1;
  while (!assemble->IsFinished()) {
    u[step] = u[step - 1];
    if (!OneTimeStep(assemble, u[step])) return false;
    step++;
  }
  Destruct();
  mout.EndBlock(verbose == 0);
  return false;
}

/*
 * NonLinearTimeIntegrator
 */

void NonLinearTimeIntegrator::PrintInfo() const {
  mout.PrintInfo("Time Integrator", verbose, PrintInfoEntry("Name", Name(), 1),
                 PrintInfoEntry("Order", GetOrder(), 1));
}

/*
 * TimeIntegratorCreator
 */

TimeIntegrator *TimeIntegratorCreator::Create() {
  try {
    switch (type) {
    case EXPLICIT_EULER:
      return new ExplicitEuler(linearSolvers.at(0));
    case HEUN:
      return new Heun(linearSolvers.at(0));
    case RUNGE:
      return new Runge(linearSolvers.at(0));
    case RUNGE_KUTTA:
      return new RungeKutta4(linearSolvers.at(0));
    case IMPLICIT_EULER: // Todo figure out how many are really needed
      return new GenericImplicitEuler(linearSolvers.at(0), linearSolvers.at(1), linearSolvers.at(2),
                                      linearSolvers.at(3));
    case IMPLICIT_MIDPOINT:
      return new ImplicitMidpointRule(linearSolvers.at(0), linearSolvers.at(1), linearSolvers.at(2),
                                      linearSolvers.at(3));
    case CRANK_NICOLSON:
      return new CrankNicolson(linearSolvers.at(0), linearSolvers.at(1), linearSolvers.at(2),
                               linearSolvers.at(3));
    case CRANK_NICOLSON_HARDCODED:
      return new CrankNicolsonHardcoded(linearSolvers.at(0), linearSolvers.at(1),
                                        linearSolvers.at(2), linearSolvers.at(3));

    case DI_RK_ORDER2:
      return new DiagonalImplicitRungeKutta2(linearSolvers.at(0), linearSolvers.at(1),
                                             linearSolvers.at(2), linearSolvers.at(3));
    case DI_RK_ORDER3:
      return new DiagonalImplicitRungeKutta3(linearSolvers.at(0), linearSolvers.at(1),
                                             linearSolvers.at(2), linearSolvers.at(3));
    case EXPONENTIAL_EULER:
      return new ExponentialEuler(linearSolvers.at(0));
    case EXPONENTIAL_MIDPOINT:
      return new ExponentialMidpoint(linearSolvers.at(0));
    case EXPONENTIAL_SHIFT_MIDPOINT:
      return new ExponentialMidpointWithShift(linearSolvers.at(0), linearSolvers.at(1));
    case EXPONENTIAL_SHIFT_EULER:
      return new ExponentialEulerWithShift(linearSolvers.at(0), linearSolvers.at(1));
    default:
      Exit(to_string(type) + " TimeIntegrator not implemented")
    }
  } catch (std::out_of_range) { THROW("Number of solvers does not match selected time integrator") }
}

std::unique_ptr<TimeIntegrator> TimeIntegratorCreator::CreateUnique() {
  return std::unique_ptr<TimeIntegrator>(Create());
}

NonLinearTimeIntegrator *TimeIntegratorCreator::CreateNonLinearTimeIntegrator() {
  try {
    switch (type) {
    case IMPLICIT_EULER:
      return new ImplicitEuler(nonLinearSolvers.at(0));
    default:
      Exit(to_string(type) + " TimeIntegrator not implemented")
    }
  } catch (std::out_of_range) { THROW("Number of solvers does not match selected time integrator") }
}

std::unique_ptr<NonLinearTimeIntegrator>
TimeIntegratorCreator::CreateNonLinearTimeIntegratorUnique() {
  return std::unique_ptr<NonLinearTimeIntegrator>(CreateNonLinearTimeIntegrator());
}
