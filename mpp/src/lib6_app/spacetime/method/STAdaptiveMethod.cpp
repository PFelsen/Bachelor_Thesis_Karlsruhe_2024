#include "STAdaptiveMethod.hpp"
#include "ErrorEstimator.hpp"

STAdaptiveMethod::STAdaptiveMethod() : STAdaptiveMethod(STMainBuilder().ReadConfig()) {};

STAdaptiveMethod::STAdaptiveMethod(STMainBuilder builder) : STMethod(builder) {
  refinement = builder.refinement_steps;
}

void STAdaptiveMethod::run() {

  bool checkSolution = false;
  Config::Get("CheckSolution", checkSolution);

  for (LevelPair level = levels;
       level.adaptivityLevel <= refinement;
       level = level.NextInAdaptivity()) {

    if (level.adaptivityLevel > 0 && Config::Exists("Linear2Preconditioner")){
      solver = std::unique_ptr<LinearSolver>(GetLinearSolverByPrefix("Linear2"));
    }

    mout << "Adaptivity on level " << level << endl;
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
    Vector RHS_orig = RHS;
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
    results["Time_LinearSolver_relative"].append(t_solve.t / t_ass.t);
    results["Time_LinearSolver"].append(t_solve.t);

    oldSolution = std::make_unique<Vector>(U);

    printVTK(U, *assemble, level.adaptivityLevel);
    computeNorms(U, B, RHS);

    Vector Eta = ee->EstimateError(B, U, level.adaptivityLevel);
    results["ErrorEstimator_2"].append(Eta.norm());

    if (checkSolution) {
      double old_reduction = solver->GetReduction();
      double old_epsilob = solver->GetEpsilon();
      double reduction = 1e-50;
      double epsilon = solver->GetIteration().GetLastDefect() * 0.1;
      Config::Get("Linear2SolverReduction", reduction);
      Config::Get("Linear2SolverEpsilon", epsilon);

      solver->SetReduction(reduction);
      solver->SetEpsilon(epsilon);

      RHS = RHS_orig - B * U;
      B.applications = 0;
      Start = Date();
      U += (*solver) * RHS;
      t_solve = Date() - Start;
      solver->SetReduction(old_reduction);
      solver->SetEpsilon(old_epsilob);

      computeNorms(U, B, RHS_orig);
      results["LinearSolverSteps"].append(solver->GetIteration().Steps());
      results["LinearSolverMatrixMultiplications"].append(B.applications);
      results["Time_Assemble"].append(-1.0);
      results["Time_LinearSolver_prep"].append(-1.0);
      results["Time_LinearSolver_relative"].append(-1.0);
      results["Time_LinearSolver"].append(t_solve.t);

      Eta = ee->EstimateError(B, U, level.adaptivityLevel);
      results["ErrorEstimator_2"].append(Eta.norm());

    }

    if (level.adaptivityLevel < refinement) {
      handleAdaptivity(Eta);
    }


  }
  results.PrintInfo();

}
