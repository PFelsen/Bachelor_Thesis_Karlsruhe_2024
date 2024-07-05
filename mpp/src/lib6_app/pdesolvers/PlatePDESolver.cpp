#include "PlatePDESolver.hpp"

PlatePDESolver::PlatePDESolver(const PDESolverConfig &conf) :
    newton(CreateNewton(conf.linearSolver, conf.preconditioner)), GenericPDESolver(conf) {}

void PlatePDESolver::run(Solution &solution) const {
  solution.vector.SetAccumulateFlag(true);
  newton->operator()(*assemble, solution.vector);
  solution.converged = newton->converged();
  int N = 10;
  Config::Get("eigenvalues", N);
  Vectors U(N, assemble->GetSharedDisc());
  for (int i = 0; i < N; ++i)
    assemble->Initialize(U[i]);
  Matrix A(solution.vector);
  Matrix B(solution.vector, false);
  assemble->Matrices(A, B);
  IEigenSolver *esolver = EigenSolverCreator("LOBPCG").Create();
  Eigenvalues lambda;
  (*esolver)(U, lambda, A, B);
  delete esolver;
  for (int i = 0; i < N; ++i)
    mout << "eigenvalue[" << i << "] = " << lambda[i] << endl;
  for (int i = 0; i < N; ++i) {
    std::string filename = "U." + to_string(i);
    mpp::plot(filename) << U[i] << mpp::endp;
  }
}

void PlatePDESolver::computeValues(Solution &solution) const {
  solution.values = ComputeValues(solution.vector);
}

void PlatePDESolver::plotVtu(Solution &solution) const {
  if (!plotting) return;
  mpp::plot("u") << solution.vector << mpp::endp;
}

void PlatePDESolver::createAssemble(std::shared_ptr<PlateProblem> problem) {
  try {
    assemble = CreatePlateAssembleUnique(*std::dynamic_pointer_cast<PlateProblem>(problem), conf);
  } catch (const std::bad_cast &ex) { THROW("PDESolver can't solve problems of given type") }
}

std::map<std::string, double> PlatePDESolver::ComputeValues(const Vector &u) const {
  std::map<std::string, double> values{};

  values["DoFCount"] = u.size();
  values["LinearSteps"] = newton->GetLinearSolver().GetIteration().Steps();

  values["Energy"] = assemble->Energy(u);

  return values;
}

void PlatePDESolver::PrintValues(const Solution &solution) {
  if (verbose == 0) return;
  else mout << "max(u) = " << solution.vector.Max() << endl;
}