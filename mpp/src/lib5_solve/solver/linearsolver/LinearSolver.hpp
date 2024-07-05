#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_

#include "Algebra.hpp"
#include "Preconditioner.hpp"
#include "RMatrix.hpp"
#include "RVector.hpp"

class CallbackStrategy {
public:
  virtual bool isActive(int restartIter, int totalIter, double residual) = 0;
  virtual double checkCurrentSolution(const Matrix &M, const Vector &u) = 0;
  virtual bool isFinished(double residual) const = 0;
  virtual ~CallbackStrategy() = default;
};

struct Iteration {
  std::vector<double> defect;

  void Init(double d0) {
    defect.clear();
    defect.push_back(d0);
  }

  Iteration() {}

  void push_back(double d) { defect.push_back(d); }

  void replace_back(double d) {
    defect.back() = d;
    //        defect.pop_back();
  }

  double GetLastDefect() { return defect.at(Steps()); }

  int Steps() const { return defect.size() - 1; }

  double Rate() const {
    if (abs(defect.back()) > VeryLarge) THROW("No convergence in iteration")
    if (defect.at(0) != 0) return pow(defect.back() / defect.at(0), 1.0 / (defect.size() - 1));
    return 0;
  }

  template<typename T>
  friend LogTextStream<T> &operator<<(LogTextStream<T> &s, const Iteration &iteration) {
    s << "d(" << iteration.defect.size() - 1 << ")"
      << "= " << iteration.defect.back() << " ";
    if (iteration.defect.size() > 1) s << "rate=" << iteration.Rate();
    s << endl;
    return s;
  }
};

class LinearSolver : public Operator {
protected:
  mutable Iteration iteration;
  int verbose = 1;

  int min_iter = 0;

  int max_iter = 800; // 200;

  double defaultEpsilon = 1e-10;

  double defaultReduction = 1e-12; // 1e-5;

  int printSteps = 1;

  std::string prefix = "Linear";

  std::unique_ptr<Preconditioner> preCond;

  std::vector<std::unique_ptr<CallbackStrategy>> callbacks;

  // Todo how to manage these guys properly
  std::shared_ptr<Vector> x0 = nullptr;
  const Matrix *A = nullptr;
  const Matrix *W = nullptr;

  virtual void solve(const Operator &A, const Operator &B, Vector &u, Vector &r, double d0,
                     double epsilon) const;

  virtual void solve(const Operator &A, const Operator &B, Vectors &u, Vectors &r, double d0,
                     double epsilon) const;
public:
  explicit LinearSolver(Preconditioner *preCond, const std::string &prefix = "Linear");

  explicit LinearSolver(std::unique_ptr<Preconditioner> preCond,
                        const std::string &prefix = "Linear");

  virtual ~LinearSolver() override = default;

  virtual std::string Name() const { return "LinearSolver"; };

  virtual void Solve(Vector &u, const Operator &A, const Operator &B, Vector &r) const;

  virtual void Solve(Vectors &u, const Operator &A, const Operator &B, Vectors &r) const;

  virtual void Solve(Vector &u, const Operator &A, const Operator &B, Vector &&r) const;

  virtual void SolveWithInitialValue(Vector &u, const Vector &b) const;

  virtual void SolveWithInitialValue(Vectors &u, const Vectors &b) const;

  virtual const Preconditioner &GetPreconditioner() { return *preCond; }

  virtual const int &GetLinearVerbose() { return verbose; }

  virtual const int &GetLinearMinIter() { return min_iter; }

  virtual const int &GetLinearMaxIter() { return max_iter; }

  virtual const double &GetLinearDefaultEps() { return defaultEpsilon; }

  virtual const double &GetLinearDefaultReduction() { return defaultReduction; }

  virtual const std::string &GetLinearPrefix() { return prefix; }

  virtual void operator()(const ILinearAssemble &assemble, Vector &u);

  virtual void operator()(const ILinearAssemble &assemble, Vectors &u);

  virtual LinearSolver &operator()(Preconditioner *pc);

  virtual LinearSolver &operator()(const Matrix &_A, const IAssemble &assemble,
                                   int reassemblePC = 1);

  virtual LinearSolver &operator()(const Matrix &_A, int reassemblePC = 1);

  virtual LinearSolver &operator()(const Matrix &AA, const Matrix &WW, int reassemblePC = 1);

  void multiply(Vector &u, const Vector &v) const override;

  void multiply(Vectors &u, const Vectors &v) const override;

  virtual void multiply_plus(Vector &u, const Vector &b) const override;

  virtual void multiply_plus(Vectors &u, const Vectors &b) const override;

  virtual int Iter() const;

  virtual Iteration GetIteration() const { return iteration; }

  virtual void PrintInfo() const;

  void AddCallback(std::unique_ptr<CallbackStrategy> cb) { callbacks.push_back(std::move(cb)); }

  // Todo discuss purpose -> needed e.g. in EigenSolver to disable output of LinearSolver without
  // config entry!!!
  virtual void SetVerbose(int v) {
    verbose = v;
    if (preCond) preCond->SetVerbose(v);
  }

  // Todo discuss purpose
  virtual void SetSteps(int s) { max_iter = s; }

  // Todo discuss purpose
  virtual void SetEpsilon(double s) { defaultEpsilon = s; }

  // Todo discuss purpose
  virtual void SetReduction(double s) { defaultReduction = s; }

  // Todo discuss purpose
  double GetEpsilon() const { return defaultEpsilon; }

  // Todo discuss purpose
  double GetReduction() const { return defaultReduction; }

  // Todo discuss purpose
  virtual void SetStartingValue(const Vector &v) { x0 = std::make_shared<Vector>(v); }
};

typedef LinearSolver LS;

std::string SolverName(const std::string &prefix);

LinearSolver *GetLinearSolver();

LinearSolver *GetLinearSolver(const string &name);

LinearSolver *GetLinearSolver(Preconditioner *pc);

LinearSolver *GetLinearSolver(const string &name, Preconditioner *pc);

LinearSolver *GetLinearSolver(const string &name, Preconditioner *pc, const string &prefix);

std::unique_ptr<LinearSolver> GetLinearSolverUnique(const string &name, Preconditioner *pc,
                                                    const string &prefix);

LinearSolver *GetLinearSolverByPrefix(const string &prefix);

LinearSolver *GetLinearSolver(std::unique_ptr<Preconditioner> pc);

LinearSolver *GetLinearSolver(const string &name, std::unique_ptr<Preconditioner> pc);

LinearSolver *GetLinearSolver(const string &name, std::unique_ptr<Preconditioner> pc,
                              const string &prefix);

std::unique_ptr<LinearSolver> GetLinearSolverUnique(std::unique_ptr<Preconditioner> pc);

std::unique_ptr<LinearSolver>
GetLinearSolverUnique(const string &name, std::unique_ptr<Preconditioner> pc, const string &prefix);

void ApplyLinearSolver(ILinearAssemble &assemble, Vector &u, bool mute = false);

void ApplyLinearSolver(ILinearAssemble &assemble, Vector &u, const std::string &name,
                       bool mute = false);

void ApplyLinearSolver(ILinearAssemble &assemble, Vectors &u, bool mute = false);

void ApplyLinearSolver(ILinearAssemble &assemble, Vectors &u, const std::string &name,
                       bool mute = false);

#endif // of #ifndef _LINEARSOLVER_H_
