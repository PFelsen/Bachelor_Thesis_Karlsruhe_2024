#include "STAdaptiveConvergenceMethod.hpp"
#include "ErrorEstimator.hpp"

STAdaptiveConvergenceMethod::STAdaptiveConvergenceMethod() : STAdaptiveConvergenceMethod(
    STMainBuilder().ReadConfig()) {};

STAdaptiveConvergenceMethod::STAdaptiveConvergenceMethod(STMainBuilder builder) : STMethod(builder) {
  refinement = builder.refinement_steps;
}

void STAdaptiveConvergenceMethod::run() {
  LevelPair endlevel = levels.Next();
  for (LevelPair level = meshes->PLevel(); level != endlevel; level = level.Next()) {
    for (int refstep = 0; refstep <= refinement; ++refstep, ++currentRun) {
      mout << "Starting on level "  << level
           << " adaptive  step " << refstep << "/" << refinement << endl;
      Vector U(0.0, assemble->GetDisc(), level);
      U.PrintInfo();
      U.GetMatrixGraph().PrintInfo();
      U.GetMatrixGraph().PrintMatrixMemoryInfo();
      U.GetMesh().PrintInfo();
      U.SetAccumulateFlag(true);
      Matrix B(U);
      Vector RHS(0.0, U);
      RHS.Clear();
      checkForHighDegree(level);
      Date Start;
      assemble->System(B, RHS);
      Time t_ass = Date() - Start;
      SetStartingValue(B, U, RHS);
      Start = Date();
      (*solver)(B, true);
      Time t_prep = Date() - Start;
      Start = Date();
      B.applications = 0;
      U += (*solver) * RHS;
      Time t_solve = Date() - Start;
      int steps = solver->GetIteration().Steps();

      results["LinearSolverSteps"].append(steps);
      results["LinearSolverMatrixMultiplications"].append(B.applications);
      results["Time_Assemble"].append(t_ass.t);
      results["Time_LinearSolver_prep"].append(t_prep.t);
      results["Time_LinearSolver_relative"].append(t_solve.t/t_ass.t);
      results["Time_LinearSolver"].append(t_solve.t);




      oldSolution = std::make_unique<Vector>(U);
      printVTK(U, *assemble, currentRun);
      computeNorms(U, B, RHS);
      Vector Eta = ee->EstimateError(B, U, currentRun);
      results["ErrorEstimator_2"].append(Eta.norm());
      if (refstep < refinement) {
        handleAdaptivity(Eta);
      }
    }
  }
  results.PrintInfo();

}

void STAdaptiveConvergenceMethod::handleAdaptivity(const Vector &Eta) {
  mout.StartBlock("Adaptivity");
  double theta = 0;
  double theta_min = 0;
  double theta_factor = 1.0;
  std::string refine_by;

  Config::Get("theta", theta);
  Config::Get("theta_factor", theta_factor);
  Config::Get("theta_min", theta_min);
  Config::Get("refine_by", refine_by);

  theta *= pow(theta_factor, currentRun);


  assemble->AdaptPolynomialDegrees(Eta, theta, theta_min, refine_by);
  assemble->AdaptQuadrature(Eta.Level(), -1, true);
  mout.EndBlock();
}
