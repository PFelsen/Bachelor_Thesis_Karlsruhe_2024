#include "STMethod.hpp"

#include "ErrorEstimator.hpp"
#include "Plotting.hpp"
#include "STDGDGViscoAcousticElement.hpp"
#include "STPathManager.hpp"
#include "SpaceTimeTransfer.hpp"

STMethod::STMethod(STMainBuilder builder) {

  Config::Get("MainVerbose", verbose);
  Config::Get("vtkplot", plotting);

  mout.StartBlock("CreateMesh");
  meshes = MeshesCreator(builder.meshName)
      .WithLevel(builder.level)
      .WithPLevel(builder.plevel)
      .CreateUnique();
  mout.EndBlock(verbose == 0);
  if(meshes->dim() < 3){
    globalPlotBackendCreator = std::make_unique<DefaultSpaceTimeBackend>;
  }else {
    globalPlotBackendCreator = std::make_unique<SpaceTimeBackend>;
  }


  assemble = CreateSTAssemble(builder.modelName, *meshes, builder.space_deg,
                              builder.time_deg, builder.problemName);
  assemble->PrintInfo();

  levels = {meshes->Level(), meshes->Level()};

  solver = GetLinearSolverUnique(createPC());

  std::string ee_name = builder.errorEstimator;
  if (!ee_name.empty()) {
    ee = CreateErrorEstimator(ee_name, *assemble, *solver);
  }
}

std::unique_ptr<Preconditioner> STMethod::createPC() {
  std::string pathChoice = "none";
  Config::Get("PathChoice", pathChoice);

  return createPreconditioner(pathChoice, *assemble);
}

std::string CalcDegreeDistribution(const Vector &u) {
  std::vector<int> numberOfCells(7, 0);
  for (cell c = u.cells(); c != u.cells_end(); c++) {
    DegreePair deg = u.GetDoF().GetDegree(c());
    numberOfCells[max(deg.space, deg.time)]++;
  }
  std::stringstream ss;
  ss << "(";
  for (int i = 0; i < numberOfCells.size(); i++) {
    numberOfCells[i] = PPM->SumOnCommSplit(numberOfCells[i], u.CommSplit());
    ss << numberOfCells[i];
    if (i < numberOfCells.size() - 1) {
      ss << ",";
    }
  }
  ss << ")";
  return ss.str();
}

void STMethod::computeNorms(const Vector &U, const Matrix &M,
                            const Vector &RHS) {
  mout.StartBlock("ComputeNorms");

  vector<string> excludedResults;
  Config::Get("ExcludedResults", excludedResults);

  auto includeResult = [&excludedResults](const std::string &result) {
    return std::find(excludedResults.begin(), excludedResults.end(), result) ==
           excludedResults.end();
  };

  double l2norm = 0;
  if (includeResult("L2Norm")) {
    l2norm = assemble->L2Norm(U);
    results["L2_Norm"].append(l2norm);
  }

  results["level"].append(U.Level().str());
  results["DoFs"].append(U.pSize());
  results["Cells"].append(U.GetMesh().CellCountGeometry());
  results["CellsWithDegree"].append(CalcDegreeDistribution(U));

  if (assemble->GetAbstractProblem().HasExactSolution()) {
    if (includeResult("L2Error")) {
      double l2error = assemble->L2Error(U);
      results["L2_Error"].append(l2error);
      double l1error = assemble->L1Error(U);
      results["L1_Error"].append(l1error);
      if (contains(assemble->Name(), "Transport") ||
          assemble->GetAbstractProblem().nL() == 0) {
        double dgerror = assemble->DGError(U);
        results["DG_Error"].append(dgerror);
      }
      if (l2norm > 0) {
        results["L2_Error_relative"].append(l2error / l2norm);
      }
    }
  }

  if (!includeResult("all")) {
    mout.EndBlock();
    return;
  }

  assemble->AdaptQuadrature(U.Level(), 6, true);

  if (includeResult("DGNorm")) {
    double dgnorm = assemble->DGNorm(U, M);
    results["DGNorm"].append(dgnorm);
  }
  if (includeResult("Conf") && assemble->HasConformingLagrangeInterpolation()) {
    Vector U_conf(0.0, U);
    double T = U.GetMesh().GetEndTime();
    assemble->ConformingLagrangeInterpolation(U, U_conf);

    results["||Lu_h - f||_Q"].append(
        assemble->MhalfLuMinusF(U, "Residual_u_h"));
    results["||Lu_conf - f||_Q"].append(
        assemble->MhalfLuMinusF(U_conf, "Residual_u_conf"));

    if (assemble->GetAbstractProblem().HasExactSolution()) {
      results["||Lu_exact -f||_Q"].append(
          assemble->MhalfLu_exactMinusF(U.Level()));
      results["||u_conf - u  ||_Q"].append(assemble->L2Error(U_conf));
      results["||u_conf(T)-u(T)||_Omega"].append(
          assemble->L2SpaceNormAtTimeError(U_conf, T));
      results["||u_conf(0)-u_h(0)||_Omega"].append(
          assemble->L2SpaceNormAtTime(U_conf - U, 0));
      results["||u_conf(T)-u_h(T)||_Omega"].append(
          assemble->L2SpaceNormAtTime(U_conf - U, T));
    }
  }

  if (includeResult("GoalFunctional")) {
    std::string goal_functional = "none";

    Config::Get("GoalFunctional", goal_functional);
    if (goal_functional != "none") {
      Point roi_min = Origin;
      Point roi_max = Origin;
      Config::Get("roi_min", roi_min);
      Config::Get("roi_max", roi_max);
      if (goal_functional == "linear") {
        results["GoalFunctional"].append(
            assemble->Goal_Functional(U, roi_min, roi_max));
      } else if (goal_functional == "quadratic") {
        results["GoalFunctional"].append(assemble->Energy(U, roi_min, roi_max));
      }
    }
  }

  double l1norm = 0;
  if (includeResult("L1Norm")) {
    l1norm = assemble->L1Norm(U);
    if (l1norm >= 0) {
      results["L1_Norm"].append(l1norm);
    }
  }

  if (!assemble->GetAbstractProblem().HasExactSolution()) {
    assemble->AdaptQuadrature(U.Level(), -1, true);
    mout.EndBlock();
    return;
  }

  if (includeResult("L1Error")) {
    double l1error =
        results["L1_Error"].last<double>();  // assemble->L1Error(U);
    // results["L1_Error"].append(l1error);
    if (l1norm > 0) {
      results["L1_Error_relative"].append(l1error / l1norm);
    }
  }
  double linferror = assemble->LInfError(U);
  results["LInf_Error"].append(linferror);

  if (includeResult("DGError")) {
    double dgerror = assemble->DGError(U);
    if (dgerror > 0) {
      results["DG_Error"].append(dgerror);
    }
  }

  Vector exactSol_int(0.0, U);
  Vector exactSol_proj(0.0, U);
  assemble->get_exact_solution(exactSol_int);
  assemble->get_projected_exact_solution(exactSol_proj);

  if (includeResult("exact")) {
    results["||Lu_int - f||_Q"].append(
        assemble->MhalfLuMinusF(exactSol_int, "Residual_u_int"));
    results["||Lu_proj - f||_Q"].append(
        assemble->MhalfLuMinusF(exactSol_proj, "Residual_u_proj"));

    if (includeResult("DG_Int_Error")) {
      results["DG_Int_Error"].append(assemble->DGError(exactSol_int));
    }

    if (includeResult("||u_int-u_h||_DG")) {
      results["||u_int-u_h||_DG"].append(assemble->DGNorm(exactSol_int - U, M));
    }

    if (includeResult("||u_proj-u_h||_DG")) {
      results["||u_proj-u_h||_DG"].append(
          assemble->DGNorm(exactSol_proj - U, M));
    }
  }

  Vector error = U - exactSol_int;
  Vector error_proj = U - exactSol_proj;

  if (includeResult("DiscNorm")) {
    auto discNorms = assemble->DiscNorm(M, error);
    double normW = discNorms.first;
    double normV = discNorms.second;

    results["ErrorW"].append(normW);
    results["ErrorV"].append(normV);
  }

  if (includeResult("L2SpaceAtT")) {
    bool includeL2SpaceAtT_int = includeResult("L2SpaceAtT_int");
    for (double t: (*meshes).coarse().GetTimesteps()) {
      std::stringstream s;
      s << "L2SpaceAtT" << std::to_string(t) /*.substr(0,4)*/ << "_error";
      results[s.str()].append(assemble->L2SpaceNormAtTimeError(U, t));
      if (includeL2SpaceAtT_int)
        results[s.str() + "_int"].append(
            assemble->L2SpaceNormAtTimeError(exactSol_int, t));
    }
  }

  if (includeResult("projError")) {
    results["L1_proj_Error"].append(assemble->L1Error(exactSol_proj));
    results["L2_proj_Error"].append(assemble->L2Error(exactSol_proj));
    results["LInf_proj_Error"].append(assemble->LInfError(exactSol_proj));
  }
  if (includeResult("intError")) {
    results["L1_int_Error"].append(assemble->L1Error(exactSol_int));
    results["L2_int_Error"].append(assemble->L2Error(exactSol_int));
    results["LInf_int_Error"].append(assemble->LInfError(exactSol_int));
  }
  if (includeResult("EE")) {
    ResidualErrorEstimator ee_res(*assemble, *solver);
    Vector Eta_res = ee_res.EstimateError(M, U, -1);
    results["eta_(res,h)"].append(norm(Eta_res));

    if (includeResult("Conf") &&
        assemble->HasConformingLagrangeInterpolation()) {
      L2ReliableResidualErrorEstimator ee_l2(*assemble, *solver);
      DGErrorEstimator ee_dg(*assemble, *solver);
      Vector Eta_l2 = ee_l2.EstimateError(M, U, -1);
      results["eta_(L2,h)"].append(norm(Eta_l2));
      Vector Eta_dg = ee_dg.EstimateError(M, U, -1, Eta_res, Eta_l2);
      results["eta_(dg,h)"].append(norm(Eta_dg));
    }
  }

  assemble->AdaptQuadrature(U.Level(), -1, true);
  mout.EndBlock(verbose == 0);
}

Results STMethod::GetResults() {
  return results;
}

void STMethod::plotParams(const Vector &U) {
  STDiscretizationT_DGDG paramDisc(*meshes, 0, 0, 2);
  Vector rhoAndKappa(0.0, paramDisc, U.Level());
  for (row r = rhoAndKappa.rows(); r != rhoAndKappa.rows_end(); ++r) {
    rhoAndKappa(r)[0] = assemble->GetAbstractProblem().Rho(r());
    rhoAndKappa(r)[1] = assemble->GetAbstractProblem().Kappa(r());
  }
  VtuPlot plot("Parameters");
  plot.AddData("Rho", rhoAndKappa, 0);
  plot.AddData("Kappa", rhoAndKappa, 1);
  plot.PlotFile();
}

void STMethod::plotDegrees(const Vector &U) {
  STDiscretizationT_DGDG degreeDisc(*meshes, 0, 0, 2);
  Vector degrees(0.0, degreeDisc, U.Level());
  for (row r = degrees.rows(); r != degrees.rows_end(); ++r) {
    degrees(r)[0] = U.GetDoF().get_time_deg(r());
    degrees(r)[1] = U.GetDoF().get_space_deg(r());
  }
  VtuPlot plot{"Degrees"};
  plot.AddData("TimeDegree", degrees, 0);
  plot.AddData("SpaceDegree", degrees, 1);
  plot.PlotFile();
}

void STMethod::plotSolution(const Vector &U) {
  STDiscretizationT_DGDG vtuDisc(U.GetDisc().GetMeshes(), 0, 0, assemble->GetAbstractProblem().Dim());
  Vector U_vtu(0.0, vtuDisc);

  for (cell c = U.cells(); c != U.cells_end(); ++c) {
    SpaceTimeViscoAcousticDGTElement elem(U, c, assemble->GetAbstractProblem().nL());
    row r = U_vtu.find_row(c());
    for (size_t i = 0; i < elem.GetComponentCount(); i++){
      U_vtu(r, i) = elem.EvaluateComponentGlobal(c(), U, i);
    }
  }

  VtuPlot plot{"NumericalSolution"};

  vector<int> velocity_indices(SpaceDimension);
  std::iota(begin(velocity_indices), end(velocity_indices), 0);
  plot.AddData("Velocity", U_vtu, velocity_indices);
  plot.AddData("Pressure", U_vtu, SpaceDimension);
  for (size_t i = 0; i < assemble->GetAbstractProblem().nL(); i++) {
    plot.AddData("PDamping"+std::to_string(i), U_vtu, SpaceDimension + i);
  }

  plot.PlotFile();

  if (assemble->GetAbstractProblem().HasExactSolution()) {
    U_vtu = 0.0;
    Vector exactSolution(0.0, U_vtu);
    assemble->get_exact_solution(exactSolution);
    // TODO: Plot PDamping of exact Solution
    VtuPlot exactSolutionPlot{"ExactSolution"};
    exactSolutionPlot.AddData("Velocity", exactSolution, {0, 1});
    exactSolutionPlot.AddData("Pressure", exactSolution, 2);
    exactSolutionPlot.PlotFile();
  }
}

void STMethod::plotVtu(const Vector &U) {
  if (typeid(U.GetDisc()) != typeid(STDiscretization_DGDG))
    return;
  if (plotting == 0)
    return;
  mout.StartBlock("vtuplot");
  plotParams(U);
  plotDegrees(U);
  plotSolution(U);
  mout.EndBlock();
}

void STMethod::checkForHighDegree(LevelPair levels) {
  bool highQuadDeg = false;
  Config::Get("highFaceQuadDeg", highQuadDeg);
  if (highQuadDeg) {
    bool adaptCellQuad = false;
    Config::Get("adaptCellQuad", adaptCellQuad);
    assemble->AdaptQuadrature(levels, 6, adaptCellQuad);
  }
}

void STMethod::SetStartingValue(Matrix &M, Vector &U, Vector &RHS) {
  bool use_random_initial = false;
  Config::Get("use_random_initial", use_random_initial, true);
  if (use_random_initial) {
    for (cell c = M.cells(); c != M.cells_end(); c++) {
      row r = M.find_row(c());
      for (int i = 0; i < r.NumberOfEntries(); i++) {
        U(c(), i) = 2 * (rand() / double(RAND_MAX)) - 1;
      }
    }
    RHS -= M * U;
    return;
  }
  bool set_solution = false;
  Config::Get<bool>("set_exact_solution", set_solution);
  if (set_solution) {
    if (assemble->GetAbstractProblem().HasExactSolution()) {
      U.Clear();
      assemble->get_exact_solution(U);
      RHS -= M * U;
    } else {
      mout << "Problem has not an exact solution." << endl;
    }
    return;
  }

  set_solution = false;
  Config::Get<bool>("set_projected_exact_solution", set_solution);
  if (set_solution) {
    if (assemble->GetAbstractProblem().HasExactSolution()) {
      U.Clear();
      assemble->get_projected_exact_solution(U);
      RHS -= M * U;
    } else {
      mout << "Problem has not an exact solution." << endl;
    }
    return;
  }

  bool usePrevious = false;
  Config::Get("UsePrevious", usePrevious);
  if (usePrevious && oldSolution) {
    std::string transferName = "Projection";
    Config::Get("Transfer", transferName);
    std::unique_ptr<SpaceTimeTransfer> transfer =
        GetSpaceTimeTransfer(assemble->GetDisc(), transferName);
    transfer->Prolongate(U, *oldSolution);
    Vector mu = M * U;
    RHS -= mu;
  }
}
