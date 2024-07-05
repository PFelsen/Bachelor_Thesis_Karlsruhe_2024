#ifndef _NEWTON_H_
#define _NEWTON_H_

#include "Assemble.hpp"
#include "LinearSolver.hpp"

template<typename AssembleType>
class NonLinearSolverT {
protected:
  int verbose = 1;

  int LS_iter = 3;

  int min_iter = 0;

  int max_iter = 10;

  double defaultEpsilon = 1e-10;

  double defaultReduction = 1e-5;

  std::string prefix = "NonLinear";

  int iter{0}; // todo remove
  int line_iter{0};

  double d, d_0; // todo remove

  std::unique_ptr<LinearSolver> solver;
public:
  NonLinearSolverT(std::unique_ptr<LinearSolver> &&solver, const string &prefix = "NonLinear") :
      solver(std::move(solver)) {
    Config::Get(prefix + "Verbose", verbose);
    Config::Get(prefix + "Steps", max_iter);
    Config::Get(prefix + "LineSearchSteps", LS_iter);
    Config::Get(prefix + "Epsilon", defaultEpsilon);
    Config::Get(prefix + "Reduction", defaultReduction);
    Config::Get(prefix + "MinimalStepNumber", min_iter);
  }

  void PrintInfo() const {
    mout.PrintInfo("Solver", verbose, PrintInfoEntry(prefix + " Solver", this->Name(), 0),
                   PrintInfoEntry(prefix + " Steps", max_iter, 1),
                   PrintInfoEntry(prefix + " Epsilon", defaultEpsilon, 1),
                   PrintInfoEntry(prefix + " Reduction", defaultReduction, 1),
                   PrintInfoEntry(this->solver->GetLinearPrefix() + " Solver", this->solver->Name(),
                                  0),
                   PrintInfoEntry(this->solver->GetLinearPrefix() + " Preconditioner",
                                  this->solver->GetPreconditioner().Name(), 0),
                   PrintInfoEntry(this->solver->GetLinearPrefix() + " Steps",
                                  this->solver->GetLinearMaxIter(), 1),
                   PrintInfoEntry(this->solver->GetLinearPrefix() + " Epsilon",
                                  this->solver->GetLinearDefaultEps(), 1),
                   PrintInfoEntry(this->solver->GetLinearPrefix() + " Reduction",
                                  this->solver->GetLinearDefaultReduction(), 1));
  }

  void PrintIterationInfo() const {
    mout.PrintInfo("Newton iterations", verbose, PrintInfoEntry("Steps", iter, 0),
                   PrintInfoEntry("line search steps", line_iter, 0));
  }

  virtual void Initialize(const AssembleType &, Vector &){};

  virtual void operator()(const AssembleType &, Vector &) = 0;

  virtual std::string Name() const = 0;

  ~NonLinearSolverT() {}

  // Todo make converged return value of operator
  bool converged() const { return (iter < max_iter); }

  double rate() const {
    if (abs(d) > VeryLarge) THROW("no convergence in iteration")
    if (iter)
      if (d_0 != 0) return pow(d / d_0, 1.0 / iter);
    return 0;
  }

  int Iter() const { return iter; }

  double defect() const { return d; }

  void endOut(int v) {
    vout(v) << "d(" << iter << ")= " << d << " rate " << rate() << endl;
    mout.EndBlock(verbose, v);
  }

  const LinearSolver &GetLinearSolver() const { return *solver; }

  friend std::ostream &operator<<(std::ostream &s, const NonLinearSolverT &S) {
    return s << "d(" << S.iter << ")= " << S.d << " rate " << S.rate() << endl;
  }

  virtual void FillInformation(const AssembleType &assemble,
                               std::unordered_map<std::string, double> &info) const {}

  virtual void GetDynamicVectors(std::vector<Vector> &values) const {};
};

using NonLinearSolver = NonLinearSolverT<IAssemble>;

template<typename AssembleType>
class NewtonT : public NonLinearSolverT<AssembleType> {
protected:
  using NonLinearSolverT<AssembleType>::verbose;
  using NonLinearSolverT<AssembleType>::LS_iter;
  using NonLinearSolverT<AssembleType>::min_iter;
  using NonLinearSolverT<AssembleType>::max_iter;
  using NonLinearSolverT<AssembleType>::defaultEpsilon;
  using NonLinearSolverT<AssembleType>::defaultReduction;
  using NonLinearSolverT<AssembleType>::prefix;
  using NonLinearSolverT<AssembleType>::iter; // todo remove
  using NonLinearSolverT<AssembleType>::line_iter;
  using NonLinearSolverT<AssembleType>::d;
  using NonLinearSolverT<AssembleType>::d_0; // todo remove
  using NonLinearSolverT<AssembleType>::solver;

  int suppressLS = 0;
  int JacobiUpdate = 1;

  //  double energy;
  /// Damping parameter of the newton method. d=1 means no damping
  double damping{1.0};

  bool isDamped() const { return damping < 1.0; }

  virtual double calculateResidualUpdate(const AssembleType &assemble, const Vector &u,
                                         Vector &r) const {
    vout(10) << "Step in Residual" << endl;
    double d = assemble.Residual(u, r);
    vout(10) << "Stepped out of Residual" << endl;
    return d;
  }

  virtual void calculateJacobiUpdate(const AssembleType &assemble, const Vector &u,
                                     Matrix &J) const {
    vout(10) << "Step in Jacobi" << endl;
    assemble.Jacobi(u, J);
    vout(10) << "Stepped out of Jacobi" << endl;
  }
public:
  /*
   *  Gives Newton's method with GMRES as linear solver
   *  preconditioned with an incomplete LU decomposition,
   *  if "LinearSolver" or "Preconditioner" are not found in config.
   */
  NewtonT() : NewtonT(std::unique_ptr<LinearSolver>(GetLinearSolver())) {}

  explicit NewtonT(std::unique_ptr<LinearSolver> &&linearSolver) :
      NonLinearSolverT<AssembleType>(std::move(linearSolver), "Newton") {
    Config::Get("NewtonSuppressFirstLineSearch", suppressLS);
    Config::Get("NewtonJacobiUpdate", JacobiUpdate);
    Config::Get("NewtonDamping", damping);
    double reduction = solver->GetReduction();
    Config::Get("NewtonLinearizationReduction", reduction);
    solver->SetReduction(reduction);
  }

  void operator()(const AssembleType &assemble, Vector &u) override {
    mout.StartBlock("Newton");

    vout(10) << "Step in Initialize" << endl;
    assemble.Initialize(u);
    vout(10) << "Stepped out of Initialize" << endl;

    vout(90) << "u(0)= " << endl << u << endl;

    Vector r(u);
    Matrix J(u);

    d = d_0 = calculateResidualUpdate(assemble, u, r);

    double eps = defaultEpsilon + defaultReduction * d;
    double d_previous = d;
    int LS_cnt = 0;

    int JU_cnt = 0; // counts iteration steps without Jacobian update
    Vector c(r);

    line_iter = 0;
    for (iter = 0; iter < max_iter; ++iter) {
      vout(3) << "r(" << iter << ")= " << endl << r << endl;
      if (d < eps) { break; }
      vout(1) << "d(" << iter << ")= " << d << endl;

      // Determines whether Jacobi Matrix is updated
      if (((iter - JU_cnt) >= JacobiUpdate) || (iter == 0)) {
        calculateJacobiUpdate(assemble, u, J);
        JU_cnt = iter;
        c = (*solver)(J, assemble) * r;
      } else c = (*solver)(J, assemble, 0) * r;

      vout(3) << "c(" << iter << ")= " << endl << c << endl;
      if (isDamped()) { c *= damping; }
      u -= c;
      vout(5) << "u-c " << endl << u << endl;
      d = calculateResidualUpdate(assemble, u, r);
      if (d > d_previous) {
        for (int l = 1; l <= LS_iter; ++l) {
          if (iter == 0)
            if (suppressLS) {
              vout(2) << "line search suppressed" << endl;
              break;
            }
          vout(1) << "line search " << l << ": d(" << iter << ")= " << d << endl;
          c *= 0.5;
          u += c;
          d = calculateResidualUpdate(assemble, u, r);
          line_iter++;
          if (d < d_previous) break;
        }
      }
      if (d > d_previous) {
        vout(5) << "line search unsuccessful." << endl;
        ++LS_cnt;
        if (LS_cnt == 3) {
          vout(1) << "too many line searches unsuccessful." << endl;
          iter = max_iter;
        }
      }
      d_previous = d;
    }
    this->endOut(1);
  };

  std::string Name() const override { return "Newton"; }
};

using Newton = NewtonT<IAssemble>;

bool NewtonMethod(IAssemble &assemble, Vector &u, bool mute = false);

/// Note that the assemble is not deleted here!
bool NewtonMethod(IAssemble *assemble, Vector &u, bool mute = false);

class ENewton : public NewtonT<IAssemble> {
  double sigma = 0.25;
public:
  ENewton(LinearSolver &linearSolver) :
      NewtonT<IAssemble>(std::unique_ptr<LinearSolver>(std::move(&linearSolver))) {
    Config::Get("NewtonSigma", sigma);
  }

  void operator()(const IAssemble &, Vector &) override;

  std::string Name() const override { return "ENewton"; }
};

std::shared_ptr<Newton> CreateNewton(const std::string &solverName,
                                     const std::string &preconditionerName);


#endif // of #ifndef _NEWTON_H_
