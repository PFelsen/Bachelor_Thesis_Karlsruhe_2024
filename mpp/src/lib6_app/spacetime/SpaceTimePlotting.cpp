#include "SpaceTimePlotting.hpp"


void printSingleSolutionVTK(const Vector &U,
                           const STAssemble &assemble,
                           const string &name) {
  bool vtkplot = false;
  Config::Get("vtkplot", vtkplot);
  if(!vtkplot){
    return;
  }
  mout.StartBlock("vtkplot");
  const Mesh &mesh = U.GetMesh();
  bool dgInTime = assemble.GetDisc().isDgInTime();
  int probDim = assemble.GetAbProblem().Dim();

  VtuPlot plot(mesh);

  Vector U_vtk(0.0, U);

  for (cell c = U.cells(); c != U.cells_end(); ++c) {
    row r = U.find_row(c());
    DegreePair deg = U.GetDoF().GetDegree(*c);
    int time_deg = deg.time;
    int shift = r.n() * time_deg / (time_deg + 1);
    if (!dgInTime) shift = r.n() * (time_deg - 1) / time_deg;

    const Shape &shape(assemble.GetDisc().GetCellShape(deg.space));

    Point localPT;
    transformPointGlobalToLocal(c(), localPT, c);
    int NX = 0;
    for (int p = 0; p < probDim; ++p) {
      double tmp = 0.0;
      for (int n = 0; n < U.GetDoF().NodalPointsLocal(deg.space, p); n++) {
        tmp += shape(localPT, n) * U(r)[shift + NX + n];
      }
      U_vtk(r)[p] = tmp;
      NX += U.GetDoF().NodalPointsLocal(deg.space, p);
    }
  }
  for (int i = 0; i < probDim; ++i) {
    COMPONENT comp = assemble.GetProblem().GetComponents()[i];
    plot.AddData(to_string(comp)+ "_" + name, U_vtk, i);
  }
  plot.PlotFile("U.vtu");
  mout.EndBlock();
}

void printVTK(const Vector &U, const STAssemble &assemble, int run) {
  bool vtkplot = false;
  Config::Get("vtkplot", vtkplot);
  if(!vtkplot){
    return;
  }
  mout.StartBlock("vtkplot");
  const Mesh &mesh = U.GetMesh();
  bool dgInTime = assemble.GetDisc().isDgInTime();
  int probDim = assemble.GetProblem().Dim();
  bool exact = assemble.GetProblem().HasExactSolution();

  VtuPlot plot(mesh);

  Vector ProcDistribution(0.0, U);
  Vector TimeDegreeDistribution(0.0, U);
  Vector SpaceDegreeDistribution(0.0, U);
  Vector U_vtk(0.0, U);

  for (cell c = U.cells(); c != U.cells_end(); ++c) {
    row r = U.find_row(c());
    DegreePair deg = U.GetDoF().GetDegree(*c);
    ProcDistribution(r)[0] = PPM->Proc(U.CommSplit());
    TimeDegreeDistribution(r)[0] = deg.time;
    SpaceDegreeDistribution(r)[0] = deg.space;

    int time_deg = deg.time;
    int shift = r.n() * time_deg / (time_deg + 1);
    if (!dgInTime) shift = r.n() * (time_deg - 1) / time_deg;

    const Shape &shape(assemble.GetDisc().GetCellShape(deg.space));

    Point localPT;
    transformPointGlobalToLocal(c(), localPT, c);
    int NX = 0;
    for (int p = 0; p < probDim; ++p) {
      double tmp = 0.0;
      for (int n = 0; n < U.GetDoF().NodalPointsLocal(deg.space, p); n++) {
        tmp += shape(localPT, n) * U(r)[shift + NX + n];
      }
      U_vtk(r)[p] = tmp;
      NX += U.GetDoF().NodalPointsLocal(deg.space, p);
    }
  }

  for (int i = 0; i < probDim; ++i) {
    COMPONENT comp = assemble.GetProblem().GetComponents()[i];
    plot.AddData(to_string(comp)+ "_" + to_string(run), U_vtk, i);
  }
  plot.AddData("Procs_" + to_string(run), ProcDistribution, 0);
  plot.AddData("Time_Refinements_" + to_string(run), TimeDegreeDistribution, 0);
  plot.AddData("Space_Refinements_" + to_string(run), SpaceDegreeDistribution, 0);
  plot.PlotFile("U.vtu");

  if (exact) {
    Vector Exact_Solution(0.0, U);
    assemble.get_exact_solution(Exact_Solution);
    Vector U_exact(0.0, U);
    Vector U_error(0.0, U);

    for (cell c = U.cells(); c != U.cells_end(); ++c) {
      row r = U.find_row(c());
      DegreePair deg = U.GetDoF().GetDegree(*c);

      int time_deg = deg.time;
      int shift = r.n() * time_deg / (time_deg + 1);
      if (!dgInTime) shift = r.n() * (time_deg - 1) / time_deg;

      const Shape &shape(assemble.GetDisc().GetCellShape(deg.space));

      Point localPT;
      transformPointGlobalToLocal(c(), localPT, c);
      int NX = 0;
      for (int p = 0; p < probDim; ++p) {
        double tmp2 = 0.0;
        for (int n = 0; n < U.GetDoF().NodalPointsLocal(deg.space, p); n++) {
          tmp2 += shape(localPT, n) * Exact_Solution(r)[shift + NX + n];
        }
        U_error(r)[p] = U_vtk(r)[p] - tmp2;
        U_exact(r)[p] = tmp2;
        NX += U.GetDoF().NodalPointsLocal(deg.space, p);
      }
    }
    VtuPlot error_plot(mesh);
    for (int i = 0; i < probDim; ++i) {
      COMPONENT comp = assemble.GetProblem().GetComponents()[i];
      error_plot.AddData(to_string(comp) + "_error",U_error, i);
    }
    error_plot.PlotFile("Error.vtu");

    VtuPlot exact_plot(mesh);
    for (int i = 0; i < probDim; ++i) {
      COMPONENT comp = assemble.GetProblem().GetComponents()[i];
      exact_plot.AddData(to_string(comp) + "_exact", U_exact, i);
    }
    error_plot.PlotFile("Exact.vtu");
  }


  if (contains(assemble.Name(), "Acoustic")) {
    Vector Rho(0.0, U);
    Vector Kappa(0.0, U);

    for (cell c = U.cells(); c != U.cells_end(); ++c) {
      row r = U.find_row(c());
      Rho(r)[0] = assemble.GetProblem().Rho(c());
      Kappa(r)[0] = assemble.GetProblem().Kappa(c());
    }
    VtuPlot material_plot(mesh);
    material_plot.AddData("Rho_" + std::to_string(run), Rho,  0);
    material_plot.AddData("Kappa_" + std::to_string(run), Kappa, 1);
    material_plot.PlotFile("Material.vtu");
  }
  mout.EndBlock();
}

void printVTK_Eta(const Vector &Eta,
                  const STAssemble &assemble,
                  const std::string &filename) {
  bool vtkplot = false;
  Config::Get("vtkplot", vtkplot);
  if(!vtkplot){
    return;
  }
  VtuPlot plot(Eta.GetMesh());
  plot.AddData("Eta", Eta, 0);
  plot.PlotFile(filename);
}



void printVTK_dual(const Vector &U_dual,
                   const Vector &Eta,
                   const STAssemble &assemble,
                   int run) {
  bool vtkplot = false;
  Config::Get("vtkplot", vtkplot);
  if(!vtkplot){
    return;
  }

  const ProblemBase &problem = assemble.GetProblem();
  const STDiscretization &disc = dynamic_cast<const STDiscretization&>(U_dual.GetDisc());
  const Mesh &STMesh = U_dual.GetMesh();
  bool dgInTime = disc.isDgInTime();

  int probDim = problem.Dim();
  Vector U_vtk_dual(0.0, U_dual);

  for (cell c = U_dual.cells(); c != U_dual.cells_end(); ++c) {
    row r = U_dual.find_row(c());
    DegreePair deg = U_dual.GetDoF().GetDegree(*c);
    int shift = r.n() * deg.time / (deg.time + 1);
    if (!dgInTime) shift = (deg.time - 1) * r.n() / deg.time;
    Point localPT;
    transformPointGlobalToLocal(c(), localPT, c);
    int NX = 0;
    double tmp = 0.0;
    for (int p = 0; p < probDim; ++p) {
      const Shape &shape(disc.GetCellShape(deg.space));
      for (int n = 0; n < U_dual.GetDoF().NodalPointsLocal(deg.space, p); n++) {
        tmp += shape(localPT, n) * U_dual(r)[shift + NX + n];
      }
      U_vtk_dual(r)[p] = tmp;
      NX += U_dual.GetDoF().NodalPointsLocal(deg.space, p);
    }
  }
  VtuPlot plot(STMesh);
  for (int i = 0; i < probDim; ++i) {
    COMPONENT comp = assemble.GetProblem().GetComponents()[i];
    string name = "Dual_" + to_string(comp) +"_" + to_string(run);
    plot.AddData(name, U_vtk_dual, i);
  }
  plot.PlotFile("Dual");
  VtuPlot plot_eta(STMesh);
  string name = "Est_Error_" + to_string(run);
  plot_eta.AddData(name, Eta, 0);
  plot_eta.PlotFile(name);
}
