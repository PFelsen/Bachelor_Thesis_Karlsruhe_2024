#include "LinearSolver.hpp"
#include "Vector.hpp"

#include "BiCGStab.hpp"
#include "CG.hpp"
// #include "BiCGStab2.hpp"
#include "GMRES.hpp"
// #include "GMRES2.hpp"
// #include "GMRES3.hpp"
#include "FGMRES.hpp"
#include "FullMultigrid.hpp"
#include "MINRES.hpp"

void LinearSolver::multiply(Vector &u, const Vector &b) const {
  u = 0.0;
  multiply_plus(u, b);
}

void LinearSolver::multiply(Vectors &u, const Vectors &b) const {
  u = 0;
  multiply_plus(u, b);
}

void LinearSolver::operator()(const ILinearAssemble &assemble, Vector &u) {
  assemble.Initialize(u);
  Vector rhs(u);
  Matrix systemMatrix(u);
  assemble.AssembleSystem(systemMatrix, rhs);
  u = (*this)(systemMatrix)*rhs;
}

void LinearSolver::operator()(const ILinearAssemble &assemble, Vectors &u) {
  assemble.Initialize(u);
  Vectors rhs(u);
  Matrix systemMatrix(u[0]);
  assemble.AssembleSystem(systemMatrix, rhs);
  (*this)(systemMatrix);
  for (int i = 0; i < u.size(); ++i)
    u[i] = (*this) * rhs[i];
}

LinearSolver &LinearSolver::operator()(Preconditioner *pc) {
  preCond = std::unique_ptr<Preconditioner>(pc);
  return *this;
}

LinearSolver &LinearSolver::operator()(const Matrix &_A, const IAssemble &assemble,
                                       int reassemblePC) {
  A = &_A;
  if (reassemblePC) {
    if (preCond) preCond->Destruct();
    if (preCond) preCond->Construct(*A, assemble);
  }
  return *this;
}

LinearSolver &LinearSolver::operator()(const Matrix &_A, int reassemblePC) {
  A = &_A;
  if (reassemblePC) {
    if (preCond) preCond->Destruct();
    if (preCond) preCond->Construct(*A);
  }
  return *this;
}

LinearSolver &LinearSolver::operator()(const Matrix &AA, const Matrix &WW, int reassemblePC) {
  W = &WW;
  return (*this)(AA, reassemblePC);
}

int LinearSolver::Iter() const {
  Warning("This function is deprecated. Information about the iterations are stored in"
          "an Iteration struct, which is currently only locally available in the solve "
          "method.") return 0;
}

void LinearSolver::PrintInfo() const {
  mout.PrintInfo("Solver", verbose, PrintInfoEntry(prefix + "Solver", this->Name()),
                 PrintInfoEntry(prefix + "Preconditioner", this->preCond->Name()),
                 PrintInfoEntry(prefix + "Steps", max_iter),
                 PrintInfoEntry(prefix + "Epsilon", defaultEpsilon),
                 PrintInfoEntry(prefix + "Reduction", defaultReduction));
}

void LinearSolver::solve(const Operator &A, const Operator &B, Vector &u, Vector &r, double d0,
                         double epsilon) const {
  vout(11) << "u " << u << endl;
  double d = d0;
  Vector c(0.0, u);
  int iter = 0;
  // Todo bool iter(): checks if current d smaller than epsilon and count smaller than max_iter
  for (; iter < max_iter; ++iter) {
    if (d < epsilon) break;
    if (iter % printSteps == 0) { vout(1) << iteration; }
    vout(10) << "r " << r << endl;
    c = B * r;
    vout(10) << "c " << c << endl;
    u += c;
    vout(10) << "u " << u << endl;
    r -= A * c;
    d = r.norm();

    iteration.push_back(d);

    {
      bool isAnyActive = false;
      for (auto &cb : callbacks) {
        isAnyActive |= cb->isActive(iter, iter, d);
      };
      if (isAnyActive) {
        bool isAnyFinished = false;
        for (auto &cb : callbacks) {
          if (cb->isActive(iter, iter, d)) {
            double error = cb->checkCurrentSolution((const Matrix &)A, u);
            if (error != -1) { vout(1) << "E(" << iter << ")= " << error << "\n"; }
          }
          isAnyFinished |= cb->isFinished(d);
        }
        if (isAnyFinished) { return; }
      }
    }
  }
}

void LinearSolver::solve(const Operator &A, const Operator &B, Vectors &u, Vectors &r, double d0,
                         double epsilon) const {
  long unsigned int size = u.size();
  vector<double> d2(size);
  d2 = r.norm();
  double d = *max_element(d2.begin(), d2.begin() + size);
  double d_0 = d;
  double d_min = *min_element(d2.begin(), d2.begin() + size);
  double eps = defaultEpsilon + defaultReduction * d;
  Vectors c(u);
  int iter = 0;
  for (; iter < max_iter; ++iter) {
    if (d < eps) break;
    vout(1) << "d(" << iter << ")= " << d << "  (min: " << d_min << ")" << endl;
    c = B * r;
    u += c;
    r -= A * c;
    d2 = r.norm();
    d = *max_element(d2.begin(), d2.begin() + size);
    d_min = *min_element(d2.begin(), d2.begin() + size);
  }
}

LinearSolver::LinearSolver(Preconditioner *preCond, const std::string &prefix) :
    LinearSolver(std::unique_ptr<Preconditioner>(preCond), prefix) {}

LinearSolver::LinearSolver(std::unique_ptr<Preconditioner> preCond, const std::string &prefix) :
    preCond(std::move(preCond)), prefix(prefix), printSteps(max_iter) {
  Config::Get(prefix + "Verbose", verbose);
  Config::Get(prefix + "Steps", max_iter);
  Config::Get(prefix + "Epsilon", defaultEpsilon);
  Config::Get(prefix + "Reduction", defaultReduction);
  Config::Get(prefix + "MinimalStepNumber", min_iter);
  Config::Get(prefix + "PrintSteps", printSteps);
}

void LinearSolver::multiply_plus(Vector &u, const Vector &b) const {
  Vector r = b;
  Solve(u, *A, *preCond, r);
}

void LinearSolver::multiply_plus(Vectors &u, const Vectors &b) const {
  Vectors r = b;
  Solve(u, *A, *preCond, r);
}

void LinearSolver::Solve(Vector &u, const Operator &A, const Operator &B, Vector &r) const {
  mout.StartBlock(Name());
  double d0 = norm(r);
  iteration.Init(d0);
  double epsilon = max(defaultEpsilon, defaultReduction * d0);
  solve(A, B, u, r, d0, epsilon);
  vout(1) << iteration;
  mout.EndBlock(verbose <= 0);
}

void LinearSolver::Solve(Vector &u, const Operator &A, const Operator &B, Vector &&r) const {
  mout.StartBlock(Name());
  double d0 = norm(r);
  iteration.Init(d0);
  double epsilon = max(defaultEpsilon, defaultReduction * d0);
  solve(A, B, u, r, d0, epsilon);
  vout(1) << iteration;
  mout.EndBlock(verbose <= 0);
}

void LinearSolver::Solve(Vectors &u, const Operator &A, const Operator &B, Vectors &r) const {
  mout.StartBlock(Name());
  Warning("TODO: What to do with norm(Vectors)") double d0 = norm(r[0]);
  double epsilon = max(defaultEpsilon, defaultReduction * d0);
  solve(A, B, u, r, d0, epsilon);
  vout(0) << iteration << endl;
  mout.EndBlock(verbose <= 0);
}

void LinearSolver::SolveWithInitialValue(Vector &u, const Vector &b) const {
  //    Vector r(b);
  //    r -= (*A) * u;
  Solve(u, *A, *preCond, b - (*A) * u);
}

void LinearSolver::SolveWithInitialValue(Vectors &u, const Vectors &b) const {
  Vectors r(b);

  r -= (*A) * u;
  Solve(u, *A, *preCond, r);
}

std::string SolverName(const std::string &prefix) {
  std::string name = "GMRES";
  Config::Get(prefix + "Solver", name);
  return name;
}

LinearSolver *GetLinearSolver() { return GetLinearSolver(SolverName("Linear"), GetPC()); }

std::unique_ptr<LinearSolver> GetLinearSolverUnique() {
  return std::unique_ptr<LinearSolver>(GetLinearSolver());
}

LinearSolver *GetLinearSolver(const string &name) { return GetLinearSolver(name, GetPC()); }

LinearSolver *GetLinearSolver(Preconditioner *pc) {
  return GetLinearSolver(SolverName("Linear"), pc);
}

LinearSolver *GetLinearSolver(const string &name, Preconditioner *pc) {
  return GetLinearSolver(name, pc, "Linear");
}

/*
 *
 *
 *
 */
LinearSolver *GetLinearSolver(const string &name, std::unique_ptr<Preconditioner> pc,
                              const string &prefix) {
  if (name == "LS") return new LS(std::move(pc), prefix);
  if (name == "CG") return new CG(std::move(pc), prefix); // Todo Tests urgently needed
  if (name == "CGX")
    return new CGX(std::move(pc), prefix); // Todo what is that? => Routine that can show how the
                                           // condition number affects the speed of convergence
  if (name == "pcg") return new CG(std::move(pc), prefix); // Todo remove
  if (name == "CGNE")
    return new CGNE(std::move(pc), prefix); // Todo what is that? => CG for normal equation
  if (name == "BiCGStab") return new BiCGStab(std::move(pc), prefix);
  if (name == "BiCGStab2")
    return new BiCGStab2(std::move(pc), prefix); // Todo Tests urgently needed
  if (name == "GMRES") return new GMRES(std::move(pc), prefix);

  //  if (name == "GMRES3") return new GMRES2(pc); // Todo Tests urgently needed
  //  if (name == "GMRES2") return new GMRES3(pc); // Todo Tests urgently needed
  if (name == "FGMRES") return new FGMRES(std::move(pc), prefix);
  if (name == "MINRES") return new MINRES(std::move(pc), prefix);
  if (name == "gmres") return new GMRES(std::move(pc), prefix);
  if (lowerCase(name) == "fullmultigrid") return new FullMultigrid();

  THROW("No linear solver " + name + " implemented")
}

std::unique_ptr<LinearSolver> GetLinearSolverUnique(const string &name, Preconditioner *pc,
                                                    const string &prefix) {
  return std::unique_ptr<LinearSolver>(GetLinearSolver(name, pc, prefix));
}

LinearSolver *GetLinearSolverByPrefix(const string &prefix) {
  return GetLinearSolver(SolverName(prefix), GetPCByPrefix(prefix), prefix);
}

LinearSolver *GetLinearSolver(std::unique_ptr<Preconditioner> pc) {
  return GetLinearSolver(SolverName("Linear"), std::move(pc));
}

LinearSolver *GetLinearSolver(const string &name, std::unique_ptr<Preconditioner> pc) {
  return GetLinearSolver(name, std::move(pc), "Linear");
}

LinearSolver *GetLinearSolver(const string &name, Preconditioner *pc, const string &prefix) {
  return GetLinearSolver(name, std::unique_ptr<Preconditioner>(pc), prefix);
}

std::unique_ptr<LinearSolver> GetLinearSolverUnique(std::unique_ptr<Preconditioner> pc) {
  std::string prefix = "Linear";
  return GetLinearSolverUnique(SolverName(prefix), std::move(pc), prefix);
}

std::unique_ptr<LinearSolver> GetLinearSolverUnique(const string &name,
                                                    std::unique_ptr<Preconditioner> pc,
                                                    const string &prefix) {
  return std::unique_ptr<LinearSolver>(GetLinearSolver(name, std::move(pc), prefix));
}

void applyLinearSolver(std::unique_ptr<LinearSolver> solver, ILinearAssemble &assemble, Vector &u,
                       bool mute) {
  if (!mute) solver->PrintInfo();
  if (!mute) assemble.PrintInfo(u);
  (*solver)(assemble, u);
}

void ApplyLinearSolver(ILinearAssemble &assemble, Vector &u, bool mute) {
  applyLinearSolver(std::unique_ptr<LinearSolver>(GetLinearSolver()), assemble, u, mute);
}

void ApplyLinearSolver(ILinearAssemble &assemble, Vector &u, const std::string &name, bool mute) {
  applyLinearSolver(std::unique_ptr<LinearSolver>(GetLinearSolver(name)), assemble, u, mute);
}

void applyLinearSolver(std::unique_ptr<LinearSolver> solver, ILinearAssemble &assemble, Vectors &u,
                       bool mute) {
  if (!mute) solver->PrintInfo();
  if (!mute) assemble.PrintInfo(u);
  (*solver)(assemble, u);
}

void ApplyLinearSolver(ILinearAssemble &assemble, Vectors &u, bool mute) {
  applyLinearSolver(std::unique_ptr<LinearSolver>(GetLinearSolver()), assemble, u, mute);
}

void ApplyLinearSolver(ILinearAssemble &assemble, Vectors &u, const std::string &name, bool mute) {
  applyLinearSolver(std::unique_ptr<LinearSolver>(GetLinearSolver(name)), assemble, u, mute);
}
