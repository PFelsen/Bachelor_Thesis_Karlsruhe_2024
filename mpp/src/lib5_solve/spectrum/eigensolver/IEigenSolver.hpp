#ifndef _EIGENSOLVER_H_
#define _EIGENSOLVER_H_

#include "Algebra.hpp"
#include "Eigenpair.hpp"
#include "LinearSolver.hpp"

#include <functional>

class IEigenSolver {
protected:
  string name;
  double eps;
  int maxstep;
  int verbose;
  int printSteps = 1;
  Date start;
  std::unique_ptr<LinearSolver> solver;
  std::unique_ptr<LinearSolver> solver2;

  double outputShift = 0.0;
  int outputSize = INT_MAX - 100;

  void RitzStep(Eigenfcts &u, Eigenvalues &lambda, Operator &A, Operator &B, Vector &tmp1,
                Vector &tmp2);

  void RitzStep(Eigenfcts &u, Eigenvalues &lambda, Operator &A, Operator &B, SymRMatrix &a,
                SymRMatrix &b, RMatrix &e, Vector &tmp1, Vector &tmp2);

  void print(const Eigenvalues &lambda, int length) const;

  void output(const Eigenvalues &lambda, int length, int step, bool final = false) const;

  void output(const Eigenvalues &lambda, int step = 0, bool final = false) const;

  void outputFinish(const Eigenvalues &lambda, int step = 0) const;

  explicit IEigenSolver(string name, std::unique_ptr<LinearSolver> &&solver,
                        std::unique_ptr<LinearSolver> &&solver2, int verbose, double eps,
                        int maxstep, int printSteps);

  virtual void
  operator()(Eigenfcts &, Eigenvalues &, Matrix &A, Matrix &B,
             std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
                 ritzDefect,
             Vector &tmp, int fev = 0) {
    THROW("Not implemented")
  }

  virtual void
  operator()(Eigenfcts &, Eigenvalues &, Matrix &A, Matrix &B, Operator &P,
             std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
                 ritzDefect,
             Vector &tmp, int fev = 0){THROW("Not implemented")};
public:
  virtual ~IEigenSolver() = default;

  void PrintInfo();

  void setVerbose(int verbose);

  void SetOutputShift(double shift) { outputShift = shift; }

  void SetOutputSize(int size) { outputSize = size; }

  void operator()(Eigenfcts &, Eigenvalues &, Matrix &A, Matrix &B, int evToConverge = 0);

  void operator()(Eigenfcts &, Eigenvalues &, Matrix &A, Matrix &B, Operator &P,
                  int evToConverge = 0);

  void operator()(Eigenfcts &, Eigenvalues &, Matrix &A, Matrix &B, bool initEigenfcts,
                  int evToConverge = 0);

  void operator()(Eigenfcts &, Eigenvalues &, Matrix &A, Matrix &B, Operator &P, bool initEigenfcts,
                  int evToConverge = 0);

  const std::string &Name() const { return name; }

  double MaxStep() const { return maxstep; }

  double Epsilon() const { return eps; }
};

class ILOBPCG : public IEigenSolver {
protected:
  int numAdditional{};
  int evConverged{};
  int firstStep = true;

  std::unique_ptr<SymRMatrix> a = nullptr;
  std::unique_ptr<SymRMatrix> b = nullptr;
  std::unique_ptr<RMatrix> e = nullptr;

  explicit ILOBPCG(string name, int numAdditional, std::unique_ptr<LinearSolver> &&solver,
                   std::unique_ptr<LinearSolver> &&solver2, int verbose, double eps, int maxstep,
                   int printSteps) :
      IEigenSolver(name, std::move(solver), std::move(solver2), verbose, eps, maxstep, printSteps),
      numAdditional(numAdditional) {}

  void
  operator()(Eigenfcts &, Eigenvalues &, Matrix &A, Matrix &B,
             std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
                 ritzDefect,
             Vector &tmp, int fev = 0) override;

  void
  operator()(Eigenfcts &, Eigenvalues &, Matrix &A, Matrix &B, Operator &P,
             std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
                 ritzDefect,
             Vector &tmp, int fev = 0) override;

  virtual void
  operator()(Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B,
             std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
             std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
                 ritzDefect,
             Vector &tmp, int fev);

  void setMatricesSmall(Eigenfcts &u, Vectors &a1, Operator &A, Operator &B, Vector &tmp1,
                        Vector &tmp2);

  void setMatricesLarge(Eigenfcts &u, Vectors &a1, Vectors &a2, Operator &A, Operator &B,
                        Vector &tmp1, Vector &tmp2);

  virtual void ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda,
                              Operator &A, Operator &B, Vector &tmp1, Vector &tmp2) = 0;

  virtual void update(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda,
                      Matrix &A, Matrix &B, Vector &residual,
                      std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
                      double defect, int r){THROW("Not implemented")};
};

class ILOBPCGSelective : public ILOBPCG {
protected:
  int *selection = nullptr;

  explicit ILOBPCGSelective(string name, int numAdditional, std::unique_ptr<LinearSolver> &&solver,
                            std::unique_ptr<LinearSolver> &&solver2, int verbose, double eps,
                            int maxstep, int printSteps) :
      ILOBPCG(name, numAdditional, std::move(solver), std::move(solver2), verbose, eps, maxstep,
              printSteps) {}

  ~ILOBPCGSelective() {
    if (selection) delete[] selection;
  }

  void
  operator()(Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B,
             std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
             std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
                 ritzDefect,
             Vector &tmp, int fev) override;

  void setMatricesSmallSelective(Eigenfcts &u, Vectors &a1, Operator &A, Operator &B, Vector &tmp1,
                                 Vector &tmp2);

  void setMatricesLargeSelective(Eigenfcts &u, Vectors &a1, Vectors &a2, Operator &A, Operator &B,
                                 Vector &tmp1, Vector &tmp2);

  virtual void ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda,
                              Operator &A, Operator &B, Vector &tmp1, Vector &tmp2) = 0;


  virtual bool
  update(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda, Matrix &A, Matrix &B,
         Vector &residual, Vector &tmp,
         std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
         std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
             ritzDefect,
         int fev) = 0;
};

#endif // of #ifndef _EIGENSOLVER_H_
