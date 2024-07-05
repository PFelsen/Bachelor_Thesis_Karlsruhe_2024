#include "ImplicitEuler.hpp"

bool ImplicitEuler::Method(INonLinearTimeAssemble &assemble, Vector &u) {
  mout.StartBlock("Euler");
  vout(1) << "Starting implicit Euler for non-linear problems" << endl;

  double t = assemble.FirstTStep();
  assemble.Initialize(u);
  assemble.SetInitialValue(u);
  assemble.FinishTimeStep(u);
  Vector u_old(u);
  Vector u_new(u);

  while (true) {
    assemble.InitTimeStep(u);

    if (assemble.Step() >= 2) {
      double factor = assemble.StepSize() / assemble.OldStepSize();
      u_new -= 1 * factor * u_old;
      u_new += factor * u;
    }
    nonLinSolver->operator()(assemble, u_new);

    if (NotDiverged(u, false)) {
      u_old = u;
      u = u_new;
      if (t == assemble.LastTStep()) { break; }
      t = assemble.NextTimeStep(t);
    } else {
      t = assemble.NextHalfTimeStep(t);
      u_new = u;
      mout << "Repeating time step with " << assemble.StepSize() << endl;
      continue;
    }

    if (NotDiverged(u, false)) {
      assemble.FinishTimeStep(u);
    } else {
      mout.EndBlock();
      return false;
    }
  }
  mout.EndBlock();
  return true;
}
