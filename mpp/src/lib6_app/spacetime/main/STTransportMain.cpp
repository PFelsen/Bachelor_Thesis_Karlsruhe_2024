#include "STTransportMain.hpp"
#include "GMRES.hpp"

#include "STDGDGTransportElement.hpp"
#include "StringUtil.hpp"

#include "VtuPlot.hpp"

void STTransportMain::createAssemble(std::shared_ptr<ITransportProblem> problem) {

  std::shared_ptr<Distribution> stDistribution;
  if (problem->HasFluxProblem()) {
    auto &fluxProblem = problem->GetFluxProblem();
    auto spaceMeshCreator = MeshesCreator(fluxProblem.Domain().Name()).
        WithPLevel(conf.pLevel).
        WithLevel(conf.level);
    fluxProblem.CreateMeshes(spaceMeshCreator);
    const Meshes &spaceMeshes = fluxProblem.GetMeshes();
    if (PPM->Size() > 1) {
      auto dist = spaceMeshes.Settings().distributions.at(0);
      stDistribution = std::make_shared<Distribution>("LikeSpace",dist);
    }

  } else {
    stDistribution = std::make_shared<Distribution>(conf.distribution);
  }
  MeshesCreator meshesCreator = MeshesCreator(conf.mesh).
      WithTimesteps(CreateTimestepsFromConfig(conf));

  if (PPM->Size() > 1){
    meshesCreator.WithDistribution(stDistribution);
  }

  problem->CreateMeshes(meshesCreator);


  problem->GetMeshes().PrintInfo();

  assemble = std::make_unique<STDGTransportAssemble>(
      DegreePair(conf.degree, conf.timeDegree),
      problem);
  assemble->PrintInfo();

  solver = GetLinearSolverUnique("GMRES", GetPC("PointBlockGaussSeidel"), "Linear");
}


bool STTransportMain::AdaptiveStep(Vector &u, int step) {
  int degree = 0;
  Config::Get("degree", degree);
  if (step > degree) return true;

  mout.StartBlock("Adaptivity");
  double theta = 0;
  Config::Get("theta", theta);
  if (theta == -1) return true;
  double theta_min = 0;
  double theta_factor = 1.0;
  std::string refine_by("abs_value");
  Config::Get("theta_factor", theta_factor);
  Config::Get("theta_min", theta_min);
  Config::Get("refine_by", refine_by);

  Vector Eta(u);
  theta *= pow(theta_factor, Eta.Level().adaptivityLevel);

  assemble->AdaptPolynomialDegrees(Eta, theta, theta_min, refine_by);
  assemble->AdaptQuadrature(Eta.Level().NextInAdaptivity(), -1, true);
  mout.EndBlock();

  Vector u_ref(0.0, u.GetSharedDisc(), u.Level().NextInAdaptivity());
  u_ref.Clear();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    STDGDGTransportElement elem(u, c);
    vector<Point> z = u_ref.GetDoF().GetNodalPoints(*c);
    for (int i = 0; i < z.size(); ++i)
      u_ref(c(), i) = elem.Density(z[i], u);
  }
  u_ref.Accumulate();

  Vector RHS_ref(0.0, u_ref);
  RHS_ref.PrintInfo();
  Matrix M_ref(RHS_ref);
  assemble->System(M_ref, RHS_ref);
  solver->operator()(M_ref, true);
  RHS_ref -= M_ref * u_ref;
  u_ref += (*solver) * RHS_ref;

  PlotVtu(u_ref);

  Solution sol(u_ref);
  computeValues(sol);
  PrintValues(sol);
  AdaptiveStep(u_ref, step + 1);
  return true;
}

void STTransportMain::run(Solution &solution) const {
  mout.StartBlock("STTransportMain");
  Vector RHS(0.0, solution.vector.GetSharedDisc());
  RHS.PrintInfo();
  Matrix M(RHS);
  assemble->System(M, RHS);
  solver->operator()(M, true);
  if (assemble->GetProblem().HasExactSolution()) {
    bool set_exact_solution = false;
    Config::Get("set_exact_solution", set_exact_solution);
    if (set_exact_solution)
      assemble->get_exact_solution(solution.vector);
  }
  RHS -= M * solution.vector;
  solution.vector += (*solver) * RHS;
}

void STTransportMain::computeValues(Solution &solution) const {
  std::map<std::string, std::vector<double>> &values = solution.multiValues;
  Vector &u = solution.vector;  
  double T = u.GetMesh().GetEndTime();
  double dt = 0;
  int level = u.Level().space;
  Config::Get("dt", dt);
  dt *= pow(2, -level);
  size_t N = u.GetMesh().GetTimesteps().size();
  vector<double> Energy(N);
  vector<double> Mass(N);
  vector<double> inflow(N);
  vector<double> outflow(N);
  vector<double> outflow_bnd(N);
  
  assemble->Energy(u, Energy);
  assemble->Mass(u, T, Mass, inflow, outflow, outflow_bnd);
  
  double Isum = 0;
  double Osum = 0;
  for (int n = 0; n < N; ++n) {
    Isum += inflow[n];
    Osum += outflow[n];
    values["Energy"].push_back(Energy[n]);
    values["Mass"].push_back(Mass[n]);
    values["Inflow"].push_back(inflow[n] / dt);
    values["Outflow"].push_back(outflow[n] / dt);
    values["IF_sum"].push_back(Isum);
    values["OF_sum"].push_back(Osum);    
  }
  std::map<std::string, double> &singleValues = solution.values;
  singleValues["norm"] = u.norm();
  singleValues["L2_Norm"] = assemble->L2Norm(u);
  if (assemble->GetProblem().HasExactSolution()){
    singleValues["L1_Error"] = assemble->L1Error(u);
    singleValues["L2_Error"] = assemble->L2Error(u);
    singleValues["GN_Error"] = assemble->GNError(u);
    singleValues["LInf_Error"] = assemble->LInfError(u);
  }  
}

void STTransportMain::PrintValues(const Solution &solution) {
  if (verbose == 0) return;
  const std::map<std::string, std::vector<double>> &allValues = solution.multiValues;
  const Vector &u = solution.vector;
  if (assemble->GetProblem().HasExactSolution()) {
    for (int n = 0; n < allValues.begin()->second.size(); ++n) {
      double t = solution.vector.GetMesh().GetTimesteps()[n];
      std::vector<PrintIterEntry<double>> entries{
        {"n",      double(n),                5,  1},
        {"t",      t, 11, 1},
        {"Mass",   allValues.at("Mass")[n],    11, 1},
        {"IFR",    allValues.at("Inflow")[n],  11, 1},
        {"IF_sum", allValues.at("IF_sum")[n],  11, 1},
        {"OFR",    allValues.at("Outflow")[n], 11, 1},
        {"OF_sum", allValues.at("OF_sum")[n],  11, 1},
        {"Error",  assemble->L2SpaceNormAtTimeError(solution.vector,t),  11, 1}    
      };
      mout.PrintIteration(verbose, entries);	
    }
  } else {
    for (int n = 0; n < allValues.begin()->second.size(); ++n) {
      std::vector<PrintIterEntry<double>> entries{
        {"n",      double(n),                5,  1},
        {"t",      solution.vector.GetMesh().GetTimesteps()[n], 11, 1},
        {"Mass",   allValues.at("Mass")[n],    11, 1},
        {"IFR",    allValues.at("Inflow")[n],  11, 1},
        {"IF_sum", allValues.at("IF_sum")[n],  11, 1},
        {"OFR",    allValues.at("Outflow")[n], 11, 1},
        {"OF_sum", allValues.at("OF_sum")[n],  11, 1}
    };
    mout.PrintIteration(verbose, entries);      
    }
  }
  auto values = solution.values.empty() ? ComputeValues(solution.vector) : solution.values;
  mout.PrintInfo("STTransport", verbose,
                 PrintInfoEntry("norm", values["norm"]),
                 PrintInfoEntry("L2_Norm", values["L2_Norm"]));
  if (assemble->GetProblem().HasExactSolution())
    mout.PrintInfo("STTransport", verbose,
		   PrintInfoEntry("L1_Error", values["L1_Error"]),
		   PrintInfoEntry("L2_Error", values["L2_Error"]),
		   PrintInfoEntry("GN_Error", values["GN_Error"]),
		   PrintInfoEntry("LInf_Error", values["LInf_Error"]));
  mout << "Problem " << assemble->GetProblem().Name() << endl;      
}

void STTransportMain::PlotVtu(const Vector &U) const {
  bool vtuPlot = false;
  Config::Get("VtuPlot", vtuPlot);
  if (!vtuPlot) return;

  auto vtuDisc = std::make_shared<STDiscretizationT_DGDG<>>(U.GetDisc().GetMeshes(), DegreePair{0, 0}, 1);
  Vector U_vtu(0.0, vtuDisc, U.Level());

  for (cell c = U.cells(); c != U.cells_end(); ++c) {
    STDGDGTransportElement elem(U, c);
    U_vtu(elem.r(), 0) = elem.Density(c(), U);
  }

  std::string name = "NumericalSolutionST_l" + std::to_string(U.Level().space) +
                     "_a" + std::to_string(U.Level().adaptivityLevel);

  VtuPlot plot(name, {.parallelPlotting=false,
      .plotBackendCreator=std::make_unique<DefaultSpaceTimeBackend>});
  plot.AddData("Density", U_vtu, 0);
  plot.PlotFile();

  vtuPlot = false;
  Config::Get("VtuPlotTimeSeries", vtuPlot);
  if (vtuPlot) {
    VtuPlot plot2(name, {.parallelPlotting=false,
        .plotBackendCreator=std::make_unique<SpaceTimeBackend>});
    plot2.AddData("Density", U_vtu, 0);
    plot2.PlotFile();
  }
}
