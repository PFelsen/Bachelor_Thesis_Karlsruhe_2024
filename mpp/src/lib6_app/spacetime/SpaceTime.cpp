#include "SpaceTime.hpp"
#include "PrintUtil.hpp"

#include "SpaceTimeMeasurement.hpp"
#include "SpaceTimePlotting.hpp"
#include "Results.hpp"
#include "DebuggingTools.hpp"
#include "LinearSolver.hpp"
#include <memory>
#include <ostream>
#include "STGLViscoAcousticAssemble.hpp"
#include "STPathManager.hpp"
#include "MemoryLogger.hpp"


/*
 * Todo Split up this main ! ! !
 * level loop could also be handled in Test or in python file
 *
 */

void SpaceTimeMain::run() {
  int run1 = 0;

  mout << endl << Level << " global refinement of " << Level << endl << endl;

  std::string shapes = "EQ";
  Config::Get("nodalPoints", shapes);

  bool previous_solution = false;
  Config::Get("PreviousSolution", previous_solution);

  std::shared_ptr<STDiscretization> prev_disc;
  Vector *prev_v;

  for (int level = meshes->pLevel(); level <= Level; ++level) {
    LevelPair levels = {level, level};
    int L = level - meshes->pLevel();
    mout << endl << endl << " ----- Set level " << level << " of " << Level
         << " ----- " << endl << endl;

    assemble = CreateSTAssemble(modelName, *meshes, sDeg, tDeg, probName);
    exact = assemble->GetProblem().HasExactSolution();

    for (int run = 0; run <= refinement; ++run, ++run1) {
      //checkForRestore(levels, run);

      mout << endl
           << run << " adaptive refinement of " << refinement << " begins at "
           << Date() << endl;

      if (run > 0 && load_balancing && level != meshes->pLevel()) {
        assemble->RedistributeCells(levels);
      }
      std::string pathChoice = "none";
      Config::Get("PathChoice", pathChoice);
      std::unique_ptr<Preconditioner> PC0 = createPreconditioner(pathChoice, levels, *assemble);

      Preconditioner *PC = PC0.get();
      LinearSolver *S = GetLinearSolver(PC);

      Vector U(0.0, assemble->GetDisc(), level, level);
      U.SetAccumulateFlag(true);


      results["DoFs"].append(U.pSize());

      U.PrintInfo();

      // Todo adapt verbose such that desired info is logged
      U.GetMesh().PrintInfo();
      results["Cells"].append(U.GetMesh().CellCountGeometry());

      Vector RHS(0.0, U);
      Matrix B(RHS);
      B = 0;
      assemble->System(B, RHS);

      if (useSparseMatrix) B.CreateSparse();
      Date Start;
      (*S)(B, true);
      mout.PrintInfo("Solver", true,PrintInfoEntry<Time>("init time", Date() - Start));

      bool set_solution = false;
      Config::Get<bool>("set_exact_solution", set_solution);
      if (set_solution) {
        if (assemble->GetProblem().HasExactSolution()) {
          assemble->get_exact_solution(U);
          RHS -= B * U;
        } else {
          mout << "Problem has not an exact solution." << endl;
        }
      } else if (previous_solution && prev_disc) {
        RHS -= B * U;
      }
      RHS.SetAccumulateFlag(false);
      U += (*S) * RHS;

      results["VNormSolution"].append(assemble->VNorm(B, U));

      if (previous_solution) {
        prev_v = new Vector(U);
        mout << "previous Norm: " << assemble->L2Norm(*prev_v) << endl;
      }

      mout << "Total SystemMemory: " << MemoryLogger::TotalMemoryUsage(false) << " MB" << endl;
      handleExactSolution(U, B, run1, results);

      if (doMeasure) {
        measurePkt(U, assemble->GetProblem(), assemble->GetDisc(), run);
      }

      if (vtkplot) {
        printVTK(U, *assemble,  run1);
      }

      double E_U = 0;
      if (functional_string == "linear")
        E_U = assemble->Goal_Functional(U, roi_min, roi_max);
      else if (functional_string == "quadratic")
        E_U = assemble->Energy(U, roi_min, roi_max);

      if (functional_string != "none") {
        results["GoalFunctional"].append(E_U);
        mout << std::defaultfloat;
        mout << "GoalFunctional " << functional_string;
        mout << " on " << roi_min << " x " << roi_max << endl;
        mout << "E(U_h) = " << std::scientific << E_U << endl;
        PrintValues("GoalFunctional E(U)", results["GoalFunctional"].asType<double>(), run);
      }
      //mout << endl;
      //mout.printEntryRaw("space-time Cells", to_string(errors.NCells));
      //mout.printEntryRaw("space-time DoFs", to_string(errors.NDoFs));
      //mout << std::defaultfloat;
      if (theta < 0) {
        mout << "uniform degree increase." << endl;
        assemble->AdaptPolynomialDegreesUniformly(levels);
        assemble->AdaptQuadrature(levels, -1, true);
      } else if (run < refinement) {
        Vector Dual_U(0.0, U);
        Vector Dual_RHS(0.0, U);

        assemble->MakeAdjointMatrix(B);
        if (useSparseMatrix) B.CreateSparse();

        if (functional_string == "linear")
          assemble->DualRHS_linear(Dual_RHS, U, roi_min, roi_max);
        if (functional_string == "quadratic")
          assemble->DualRHS_quadratic(Dual_RHS, U, roi_min, roi_max);
        Dual_RHS.Collect();


        PC->transpose();

        LinearSolver *T = GetLinearSolver(PC);
        (*T)(B);
        (*T).SetEpsilon(dual_eps);
        (*T).SetReduction(dual_red);
        mout << "Solve dual problem for p-adaptivity" << endl;
        Dual_U = (*T) * Dual_RHS;

        Vector Eta(0.0, Dual_U);
        double eta_max = assemble->DualErrorEstimateJump(U, Dual_U, Eta);


        printVTK_dual(Dual_U, Eta, *assemble, run1);
        assemble->AdaptPolynomialDegrees(Eta, theta, theta_min, refine_by);
        assemble->AdaptQuadrature(Eta.Level(), -1, true);
      }
    }
  }

  if (!(exact)) return;

  results.PrintInfo();

  /*PrintValues("W_Error", errors.W);
  PrintValues("V_Error", errors.V);
  PrintValues("L1_Error", errors.L1);
  PrintValues("L2_Error", errors.L2);
  PrintValues("LInf_Error", errors.LInf);
  PrintValues("L1_int_Error", errors.L1_int);
  PrintValues("L2_int_Error", errors.L2_int);
  PrintValues("LInf_int_Error", errors.LInf_int);
  PrintValues("GN_Error", errors.GN);
  PrintValues("DGSemiNormExact", errors.DGSemiNormExactSolution);

  for (auto &[key, value]: errors.errors) {
    RVector err = RVector(value.size());
    for (int i = 0; i < value.size(); i++) {
      err[i] = value[i];
    }
    PrintValues(key, err);
  }*/

  //PrintValues("DGSemiNorm           ", errors.DGSemiNormSolution);



}

void SpaceTimeMain::handleExactSolution(Vector &U,
                                        Matrix &B,
                                        const int run1,
                                        Results &results) {
  if(!assemble->GetProblem().HasExactSolution()){
    return;
  }

  mout.StartBlock("CalcErrors");

  assemble->AdaptQuadrature(U.Level(), 6, true);

  double l1error = assemble->L1Error(U);
  if(l1error > 0){
    results["L1_Error"].append(l1error);
    double l1norm = assemble->L1Norm(U);
    if(l1norm > 0){
      results["L1_Error_relative"].append(l1error / l1norm);
    }
  }

  double l2error = assemble->L2Error(U);
  if(l2error > 0){
    results["L2_Error"].append(l2error);
    double l2norm = assemble->L2Norm(U);
    if(l2norm > 0){
      results["L2_Error_relative"].append(l2error / l2norm);
    }
  }
  double linferror = assemble->LInfError(U);
  if(linferror > 0){
    results["LInf_Error"].append(linferror);
  }

  //errors.GN[run1] = assemble->GNError(U);

  Vector exactSol(0.0, U);
  assemble->get_exact_solution(exactSol);
  Vector error = U - exactSol;


  double t = 0.91;
  std::string t_str = std::to_string(t);
  results["L2SpaceAtT"+t_str+"_error"].append(assemble->L2SpaceNormAtTimeError(U, t));
  results["L2SpaceAtT"+t_str+"_int_error"].append(assemble->L2SpaceNormAtTime(error, t));

  results["L2SpaceAtT1_error"].append(assemble->L2SpaceNormAtTimeError(U, 1));
  results["L2SpaceAtT1_int_error"].append(assemble->L2SpaceNormAtTime(error, 1));

  //results["VNormIntExact"].append(assemble->VNorm(B, U));
  //results["VNormError"].append(assemble->VNorm(B, error));


  results["L1_int_Error"].append(assemble->L1Error(exactSol));
  results["L2_int_Error"].append(assemble->L2Error(exactSol));
  results["LInf_int_Error"].append(assemble->LInfError(exactSol));

  mout.EndBlock();

  assemble->AdaptQuadrature(U.Level(), -1, true);
}
/*
void SpaceTimeMain::checkForRestore(const LevelPair levels, int &run) {
  if (restore) {
    restore = false;
    mout.StartBlock("Restore");
    restoreAdapt(meshes, *assemble->GetDiscSharedPtr(), levels, run);

    assemble->GetDisc().communicate(levels);
    assemble->GetDiscSharedPtr()->adaptQuadrature(-1);
    mout.EndBlock();
  } else if (run > 0 && theta >= 0) {
    saveAdapt(meshes, *assemble->GetDiscSharedPtr(), levels, run);
  }
}
*/



