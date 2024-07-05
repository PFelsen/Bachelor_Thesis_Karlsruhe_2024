#ifndef SPACETIME_STACOUSTICMAIN_HPP
#define SPACETIME_STACOUSTICMAIN_HPP

#include "STMethod.hpp"
#include "STDGViscoAcousticAssemble.hpp"
#include "AcousticProblems.hpp"
#include "PDESolver.hpp"
#include "STPathManager.hpp"
#include "STDGMatrixFreeViscoAcousticAssemble.hpp"
#include "STGLViscoAcousticAssemble.hpp"
#include "STPGViscoAcousticAssemble.hpp"
#include "AdaptiveConfig.hpp"
#include "STCallbackStrategy.hpp"

template<class AcousticAssemble>
class STAcousticMainT : public PDESolver<AcousticProblem> {
protected:
  std::unique_ptr<AcousticAssemble> assemble = nullptr;

  std::unique_ptr<LinearSolver> solver = nullptr;

  std::unique_ptr<ObservationSpecification> specification = nullptr;

  bool measure = false;

  mutable int count = 0;

  mutable std::unique_ptr<Matrix> M = nullptr;

  // Todo: discuss if this as a state is really needed.
  // UQOC has same problem; i.e. control has to be passed too
  // Proposal: Rather give Run method possibility to take
  // Problem and old solution
  mutable std::unique_ptr<Vector> oldSolution = nullptr;

public:
  using GenericPDESolver::Run;

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return assemble->GetSharedDisc();
  }

  explicit STAcousticMainT(PDESolverConfig _mc) : GenericPDESolver(_mc) { }

  void createAssemble(std::shared_ptr<AcousticProblem> problem) override;

  void createAdjointAssemble(std::shared_ptr<AcousticProblem> problem) override {
    solver->GetPreconditioner().transpose();
    assemble->SetProblem(problem);
  };

  STAcousticMainT() : STAcousticMainT(PDESolverConfig{}) {}

  std::string Name() const override {
    return "STAcoustic";
  }

  void plotVtu(Solution &solution) const override {
    PlotVtu(solution.vector);
  }

  void computeValues(Solution &solution) const override {
    Vector &U = solution.vector;
    mout.StartBlock("Measure");
    if (measure && specification) {
      try {
        SeismogramData data = assemble->measure(*specification, U);
        solution.seismograms = std::make_shared<SeismogramData>(data);
        data.WriteToFile("data/seismo/seismogram"+ std::to_string(count));
      } catch (NotImplementedException &e) {
        mout << "Measuring not implemented.";
      }
    }

    mout.EndBlock(true);

    std::map<std::string, double> &values = solution.values;

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
      values["L2_Norm"] = l2norm;
    }

    values["SpaceLevel"] = U.Level().space;
    values["TimeLevel"] = U.Level().time;
    values["AdaptiveLevel"] = U.Level().adaptivityLevel;

    values["DoFs"] = U.pSize();
    values["MatrixAllocationSize"] = PPM->Sum(U.Size());

    values["Cells"] = U.GetMesh().CellCountGeometry();
    //value["CellsWithDegree"] = CalcDegreeDistribution(U));

    if (assemble->GetProblem().HasExactSolution()) {
      if (includeResult("L2Error")) {
        double l2error = assemble->L2Error(U);
        values["L2_Error"] = l2error;
        double l1error = assemble->L1Error(U);
        values["L1_Error"] = l1error;
        if (contains(assemble->Name(), "Transport") ||
            assemble->GetProblem().nL() == 0) {
          double dgerror = assemble->DGError(U);
          values["DG_Error"] = dgerror;
        }
        if (l2norm > 0) {
          values["L2_Error_relative"] = l2error / l2norm;
        }
      }
    }

    if (!includeResult("all")) {
      mout.EndBlock();
      return;
    }

    assemble->AdaptQuadrature(U.Level(), 6, true);

    if (includeResult("DGNorm")) {
      double dgnorm = assemble->DGNorm(U, *M);
      values["DGNorm"] = dgnorm;
    }
    if (includeResult("Conf") && assemble->HasConformingLagrangeInterpolation()) {
      Vector U_conf(0.0, U);
      double T = U.GetMesh().GetEndTime();
      assemble->ConformingLagrangeInterpolation(U, U_conf);

      values["||Lu_h - f||_Q"] = assemble->MhalfLuMinusF(U, "Residual_u_h");
      values["||Lu_conf - f||_Q"] = assemble->MhalfLuMinusF(U_conf, "Residual_u_conf");

      if (assemble->GetProblem().HasExactSolution()) {
        values["||Lu_exact -f||_Q"] = assemble->MhalfLu_exactMinusF(U.Level());
        values["||u_conf - u  ||_Q"] = assemble->L2Error(U_conf);
        values["||u_conf(T)-u(T)||_Omega"] = assemble->L2SpaceNormAtTimeError(U_conf, T);
        values["||u_conf(0)-u_h(0)||_Omega"] = assemble->L2SpaceNormAtTime(U_conf - U, 0);
        values["||u_conf(T)-u_h(T)||_Omega"] = assemble->L2SpaceNormAtTime(U_conf - U, T);
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
          values["GoalFunctional"] = assemble->Goal_Functional(U, roi_min, roi_max);
        } else if (goal_functional == "quadratic") {
          values["GoalFunctional"] = assemble->Energy(U, roi_min, roi_max);
        }
      }
    }

    double l1norm = 0;
    if (includeResult("L1Norm")) {
      l1norm = assemble->L1Norm(U);
      if (l1norm >= 0) {
        values["L1_Norm"] = l1norm;
      }
    }

    if (!assemble->GetProblem().HasExactSolution()) {
      assemble->AdaptQuadrature(U.Level(), -1, true);
      mout.EndBlock();
      return;
    }

    if (includeResult("L1Error")) {
      double l1error = values["L1_Error"];
      if (l1norm > 0) {
        values["L1_Error_relative"] = l1error / l1norm;
      }
    }
    double linferror = assemble->LInfError(U);
    values["LInf_Error"] = linferror;

    if (includeResult("DGError")) {
      double dgerror = assemble->DGError(U);
      if (dgerror > 0) {
        values["DG_Error"] = dgerror;
      }
    }

    Vector exactSol_int(0.0, U);
    Vector exactSol_proj(0.0, U);
    assemble->get_exact_solution(exactSol_int);
    assemble->get_projected_exact_solution(exactSol_proj);

    if (includeResult("exact")) {
      values["||Lu_int - f||_Q"] = assemble->MhalfLuMinusF(exactSol_int, "Residual_u_int");
      values["||Lu_proj - f||_Q"] = assemble->MhalfLuMinusF(exactSol_proj, "Residual_u_proj");

      if (includeResult("DG_Int_Error")) {
        values["DG_Int_Error"] = assemble->DGError(exactSol_int);
      }

      if (includeResult("||u_int-u_h||_DG")) {
        values["||u_int-u_h||_DG"] = assemble->DGNorm(exactSol_int - U, *M);
      }

      if (includeResult("||u_proj-u_h||_DG")) {
        values["||u_proj-u_h||_DG"] = assemble->DGNorm(exactSol_proj - U, *M);
      }
    }

    Vector error = U - exactSol_int;
    Vector error_proj = U - exactSol_proj;


    if (includeResult("L2SpaceAtT")) {
      bool includeL2SpaceAtT_int = includeResult("L2SpaceAtT_int");
      for (double t: U.GetMeshes().coarse().GetTimesteps()) {
        std::stringstream s;
        s << "L2SpaceAtT" << std::to_string(t) /*.substr(0,4)*/ << "_error";
        values[s.str()] = assemble->L2SpaceNormAtTimeError(U, t);
        if (includeL2SpaceAtT_int)
          values[s.str() + "_int"] = assemble->L2SpaceNormAtTimeError(exactSol_int, t);
      }
    }

    if (includeResult("projError")) {
      values["L1_proj_Error"] = assemble->L1Error(exactSol_proj);
      values["L2_proj_Error"] = assemble->L2Error(exactSol_proj);
      values["LInf_proj_Error"] = assemble->LInfError(exactSol_proj);
    }
    if (includeResult("intError")) {
      values["L1_int_Error"] = assemble->L1Error(exactSol_int);
      values["L2_int_Error"] = assemble->L2Error(exactSol_int);
      values["LInf_int_Error"] = assemble->LInfError(exactSol_int);
    }
    if (includeResult("EE")) {
      ResidualErrorEstimator ee_res(*assemble, *solver);
      Vector Eta_res = ee_res.EstimateError(*M, U, -1);
      values["eta_(res,h)"] = norm(Eta_res);

      if (includeResult("Conf") && assemble->HasConformingLagrangeInterpolation()) {
        L2ReliableResidualErrorEstimator ee_l2(*assemble, *solver);
        DGErrorEstimator ee_dg(*assemble, *solver);
        Vector Eta_l2 = ee_l2.EstimateError(*M, U, -1);
        values["eta_(L2,h)"] = norm(Eta_l2);
        Vector Eta_dg = ee_dg.EstimateError(*M, U, -1, Eta_res, Eta_l2);
        values["eta_(dg,h)"] = norm(Eta_dg);
      }
    }

    assemble->AdaptQuadrature(U.Level(), -1, true);
    mout.EndBlock(verbose == 0);
  }

  void PlotVtu(const Vector &u) const {
    bool vtuPlot = false;
    Config::Get("VtuPlot", vtuPlot);
    if (!vtuPlot) return;
    globalPlotBackendCreator = std::make_unique<DefaultSpaceTimeBackend>;
    plotParams(u);
    plotDegrees(u);
    plotSolution(u);
    vtuPlot = false;
    Config::Get("VtuPlotTimeSeries", vtuPlot);
    if (vtuPlot) {
      globalPlotBackendCreator = std::make_unique<SpaceTimeBackend>;
      plotDegrees(u);
      plotSolution(u);
    }
    count++;
  }

  void plotParams(const Vector &U) const {
    assemble->PlotParameters(U, std::to_string(count));
  }

  void plotDegrees(const Vector &U) const {
    auto degreeDisc = std::make_shared<STDiscretizationT_DGDG<>>(U.GetMeshes(), DegreePair{0, 0}, 2);
    Vector degrees(0.0, degreeDisc, U.Level());
    for (row r = degrees.rows(); r != degrees.rows_end(); ++r) {
      degrees(r)[0] = U.GetDoF().get_time_deg(r());
      degrees(r)[1] = U.GetDoF().get_space_deg(r());
    }
    VtuPlot plot{"Degrees"+ std::to_string(count)};
    plot.AddData("TimeDegree", degrees, 0);
    plot.AddData("SpaceDegree", degrees, 1);
    plot.PlotFile();
  }

  void plotSolution(const Vector &U) const {
    assemble->PlotSingleSolution(U, "NumericalSolution" + std::to_string(count));

    if (assemble->GetAbstractProblem().HasExactSolution()) {
      Vector exactSolution(0.0, U);
      assemble->get_exact_solution(exactSolution);
      assemble->PlotSingleSolution(exactSolution, "ExactSolution" + std::to_string(count));
    }
  }

  void PrintValues(const Solution &solution) override {

    const std::map<std::string, std::vector<double>> &allValues = solution.multiValues;

    /*for (int n = 0; n < allValues.begin()->second.size(); ++n) {
      auto timestep = solution.vector.GetMesh().GetTimesteps()[n];
      std::vector<PrintIterEntry<double>>
          entries{{"n", double(n), 5, 1},
                  {"t", timestep, 11, 1},
                  //  {"Mass",   allValues["Mass"][n],           11, 1},
                  {"Energy", allValues.at("Energy")[n], 11, 1}};

      mout.PrintIteration(0, entries);

    }*/

    auto &values = solution.values;


    for (auto &[key, value]: values) {
      mout << key << ":";
      for (int i = 0; i < 40 - key.size(); i++) {
        mout << ".";
      }
      mout << value << endl;
    }
  }

  void clearOldSolution() {
    oldSolution = nullptr;
  }

  void run(Solution &solution) const override {
    Vector &u = solution.vector;
    if(oldSolution && conf.usePrevious){
      u = *oldSolution;
    }

    mout.StartBlock("STAcousticMain");
    mout << "Run" << endl;
    Vector RHS(0.0, u.GetSharedDisc(), u.Level());
    RHS.PrintInfo();

    M = std::make_unique<Matrix>(RHS);

    assemble->System(*M, RHS);
    solver->operator()(*M, true);
    RHS -= (*M) * u;
    u += (*solver) * RHS;
    if (conf.usePrevious) {
      oldSolution = std::make_unique<Vector>(u);
    }
    mout.EndBlock();
  }

  void runAdjoint(Solution &solution) const override {
    Vector &u = solution.vector;
    mout.StartBlock("STAcousticMain");
    Vector RHS(0.0, u.GetSharedDisc(), u.Level());
    if (count == 0)RHS.PrintInfo();

    if (!M){
      M = std::make_unique<Matrix>(RHS);
      assemble->System(*M, RHS);
      RHS.Clear();
    }
    assemble->RHS(RHS);

    assemble->MakeAdjointMatrix(*M);

    solver->operator()(*M, true);
    RHS -= (*M) * u;
    u += (*solver) * RHS;

    mout.EndBlock();
  }

  void EstimateErrorAndApplyAdaptivity(const Vector &U, AdaptiveConfig conf) {
    mout.StartBlock("Adaptivity");
    auto errorEstimator = CreateErrorEstimator(conf.data.errorEstimator,
                                               *assemble,
                                               *solver);
    Vector Eta = errorEstimator->EstimateError(*M, U, U.Level().adaptivityLevel);
    mout << "U.Level in EstimateErrorAndApplyAdaptivity" << U.Level() << endl;
    mout << DOUT(Eta.Level()) << endl;
    double theta = conf.data.theta;
    double theta_min = conf.data.theta_min;
    double theta_factor = conf.data.theta_factor;
    std::string refine_by = std::ref(conf.data.refine_by);
    theta *= pow(theta_factor, Eta.Level().adaptivityLevel);
    assemble->AdaptPolynomialDegrees(Eta, theta, theta_min, refine_by);
    assemble->AdaptQuadrature(Eta.Level().NextInAdaptivity(), -1, false);
    mout.EndBlock();
  }

  Vector CalculateMaterialUpdate(const Vector &material, const Vector &forwardSolution, const Vector &backwardSolution){
    return assemble->CalculateMaterialUpdate(material, forwardSolution, backwardSolution);
  }

};

template<class AcousticAssemble>
void STAcousticMainT<AcousticAssemble>::createAssemble(std::shared_ptr<AcousticProblem> problem) {
  auto acousticProblem = std::dynamic_pointer_cast<AcousticProblem>(problem);

  if(!problem->IsMeshInitialized()) {
    if (!conf.mesh.empty()) {
      problem->CreateMeshes(
          MeshesCreator(conf.mesh).WithPLevel(conf.pLevel).WithLevel(conf.level));
    } else {
      problem->CreateMeshes(MeshesCreator().WithPLevel(conf.pLevel).WithLevel(conf.level));
    }
    problem->GetMeshes().PrintInfo();
  }

  if(!assemble){
    assemble = std::make_unique<AcousticAssemble>(acousticProblem->GetMeshes(), DegreePair(conf.degree, conf.timeDegree), acousticProblem);

  }else {
    assemble->SetProblem(acousticProblem);
  }
  std::string pathChoice = "none";
  Config::Get("PathChoice", pathChoice);
  solver = GetLinearSolverUnique(createPreconditioner(pathChoice, *assemble));

  //assemble->PrintInfo();

  Config::Get("doMeasure", measure);


  if (measure && !specification) {
    const Mesh &fineMesh = problem->GetMeshes()[problem->GetMeshes().FineLevel()];
    double dt = fineMesh.GetTimesteps()[1] / 4;
    specification = std::move(ObservationSpecificationBuilder().
        WithReceiverCounts({20}).WithT(fineMesh.GetEndTime()).Withdt(dt).
        WithWaypoints({Point(750, -350), Point(750, 1650)}).CreateUnique());
  }
}

using STAcousticMain = STAcousticMainT<STDGViscoAcousticAssemble>;
using STAcousticMainGL = STAcousticMainT<STGLViscoAcousticAssemble>;
using STAcousticMainMF = STAcousticMainT<STDGMatrixFreeViscoAcousticAssemble>;
using STAcousticMainPG = STAcousticMainT<STPGViscoAcousticAssemble>;


#endif //SPACETIME_STACOUSTICMAIN_HPP
