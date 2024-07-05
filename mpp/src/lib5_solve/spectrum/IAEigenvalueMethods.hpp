#ifndef IAEIGENVALUEMETHODS_H
#define IAEIGENVALUEMETHODS_H

#include "EigenSolverCreator.hpp"
#include "IAEigenvalueAssemble.hpp"
#include "IASpectrum.hpp"

/*!
 * EigenvalueMethod
 *
 * This is the basic class for all eigenvalue methods
 */
template<bool VERIFIED>
class IAEigenvalueMethod {
protected:
  Date start;
  std::string name;
  int verbose;
  int minESSize;
  int addEVs;

  std::unique_ptr<IEigenSolver> ES = nullptr;
  std::unique_ptr<IEigenSolver> ES_coarse = nullptr;
public:
  explicit IAEigenvalueMethod(const std::string &name = "EMethod");

  const std::string &Name() const { return name; }

  void SetVerbose(int v) { verbose = v; };

  constexpr bool IsVerified() const { return VERIFIED; }
protected:
  //! Realizes the computation of approximative eigenpairs
  void computeEigenpairs(IAEigenvalueAssemble<VERIFIED> &assemble, IEigenSolver &esolver,
                         Eigenfcts &U, Eigenvalues &lambda, Matrix &A, Matrix &B, const double &t,
                         bool initBC);

  void computeEigenpairs(IAEigenvalueAssemble<VERIFIED> &assemble, IEigenSolver &esolver,
                         Eigenpairs &EP, Matrix &A, Matrix &B, const double &t, bool initBC);

  void computeEigenpairs(IAEigenvalueAssemble<VERIFIED> &assemble, IEigenSolver &esolver,
                         Eigenpairs &EP, const double &t);

  void computeEigenfcts(IAEigenvalueAssemble<VERIFIED> &assemble, IEigenSolver &esolver,
                        Eigenfcts &U, const double &t);
};

/*!
 * RayleighRitzMethod
 *
 * Implements the Rayleigh Ritz methods to compute upper bounds for the eigenvalues of the
 * eigenvalue problem
 *
 *      \f[ M_t(u,v) = \lambda N_t(u,v) \qquad (\forall v), \f]
 *
 * where t is a parameter which e.g. can be used in a homotopy.
 */
template<bool VERIFIED>
class IARayleighRitzMethod : public IAEigenvalueMethod<VERIFIED> {
public:
  IARayleighRitzMethod() : IAEigenvalueMethod<VERIFIED>("RayleighRitz") {}

  //! Computes (verified) upper bounds depending on template parameter VERIFIED
  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, Eigenfcts &U,
                  IAUpperEigenvalueBounds &Lambda, const double &t = 0.0) {
    (*this)(assemble, U, Lambda, true, t);
  }

  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, Eigenfcts &U,
                  IAUpperEigenvalueBounds &Lambda, bool initEF, const double &t = 0.0);
};

/*!
 * LehmannGoerischMethod
 *
 * Implements the Lehmann Goerisch methods to compute lower bounds for the eigenvalues of the
 * eigenvalue problem
 *
 *      \f[ M_t(u,v) = \lambda N_t(u,v) \qquad (\forall v), \f]
 *
 * where t is a parameter which e.g. can be used in a homotopy.
 */
template<bool VERIFIED>
class IALehmannGoerischMethod : public IAEigenvalueMethod<VERIFIED> {
  std::unique_ptr<LinearSolver> solver = nullptr;
public:
  IALehmannGoerischMethod();

  //! Computes (verified) lower bounds depending on template parameter VERIFIED
  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, int levelGoerisch, Eigenpairs &EP,
                  IALowerEigenvalueBounds &lambda, const double &rho, const double &t = 0.0) {
    (*this)(assemble, levelGoerisch, EP, lambda, rho, true, t);
  }

  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, int levelGoerisch, Eigenpairs &EP,
                  IALowerEigenvalueBounds &lambda, const double &rho, bool initEP,
                  const double &t = 0.0);

  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, Eigenpairs &EP,
                  IALowerEigenvalueBounds &lambda, const double &rho, const double &t = 0.0) {
    (*this)(assemble, -1, EP, lambda, rho, t);
  }

  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, Eigenpairs &EP,
                  IALowerEigenvalueBounds &lambda, const double &rho, bool initEP,
                  const double &t = 0.0) {
    (*this)(assemble, -1, EP, lambda, rho, initEP, t);
  }
private:
  void computeGoerischFcts(IAEigenvalueAssemble<VERIFIED> &assemble, const Eigenpairs &ep,
                           Goerischfcts &w, double t);
};

/*!
 * IAHomotopyMethod
 */
template<bool VERIFIED>
class IAHomotopyMethod : public IAEigenvalueMethod<VERIFIED> {
protected:
  IARayleighRitzMethod<VERIFIED> RM;
  IALehmannGoerischMethod<VERIFIED> GM;

  double stepMin;
  double stepMax;
  double stepSize;
  double step;
  double separationRho;
  double nearRho;
  double separationCluster;
  int numBisectionSteps;
  int numBisectionSteps_coarse;
  int coarseLevel;

  // if non-negative, the homotopy stops if the first appr. eigenvalue is larger than nearZero
  double nearZero = -1.0;
protected:
  void readConfig(bool all);
public:
  explicit IAHomotopyMethod(double stepMin, double stepMax, double stepSize, int coarseLevel,
                            const std::string &name = "");

  explicit IAHomotopyMethod(int coarseLevel, const std::string &name = "");

  explicit IAHomotopyMethod(double stepMin, double stepMax, double stepSize,
                            const std::string &name = "");

  explicit IAHomotopyMethod(const std::string &name = "");

  //! Computes (verified) upper eigenvalue bounds
  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, Eigenfcts &U,
                  IAUpperEigenvalueBounds &Lambda, const double &t = 0.0) {
    RM(assemble, U, Lambda, t);
  }

  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, Eigenfcts &U,
                  IAUpperEigenvalueBounds &Lambda, bool initEF, const double &t = 0.0) {
    RM(assemble, U, Lambda, initEF, t);
  }

  IARayleighRitzMethod<VERIFIED> &GetRayleighRitzMethod() { return RM; }

  //! Computes (verified) lower eigenvalue bounds
  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, int levelGoerisch, Eigenpairs &EP,
                  IALowerEigenvalueBounds &lambda, const double &rho, const double &t = 0.0) {
    GM(assemble, levelGoerisch, EP, lambda, rho, t);
  }

  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, int levelGoerisch, Eigenpairs &EP,
                  IALowerEigenvalueBounds &lambda, const double &rho, bool initEP,
                  const double &t = 0.0) {
    GM(assemble, levelGoerisch, EP, lambda, rho, initEP, t);
  }

  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, Eigenpairs &EP,
                  IALowerEigenvalueBounds &lambda, const double &rho, const double &t = 0.0) {
    GM(assemble, EP, lambda, rho, t);
  }

  void operator()(IAEigenvalueAssemble<VERIFIED> &assemble, Eigenpairs &EP,
                  IALowerEigenvalueBounds &lambda, const double &rho, bool initEP,
                  const double &t = 0.0) {
    GM(assemble, EP, lambda, rho, initEP, t);
  }

  IALehmannGoerischMethod<VERIFIED> &GetLehmannGoerischMethod() { return GM; }

  //! Performes the homotopy steps
  double operator()(IAEigenvalueAssemble<VERIFIED> &assemble, int levelGoerisch, Eigenpairs &EP,
                    double &rho, bool initEP = false);

  double operator()(IAEigenvalueAssemble<VERIFIED> &assemble, Eigenpairs &EP, double &rho,
                    bool initEP = false) {
    return (*this)(assemble, -1, EP, rho, initEP);
  }

  double StepMin() const { return stepMin; }

  void StepMin(double _stepMin) { stepMin = _stepMin; }

  double StepMax() const { return stepMax; }

  void StepMax(double _stepMax) { stepMax = _stepMax; }

  double StepSize() const { return stepSize; }

  void StepSize(double _stepSize) { stepSize = _stepSize; }

  double NearZero() const { return nearZero; };

  void NearZero(double _nearZero) { nearZero = _nearZero; }

  double SeparationRho() const { return separationRho; }

  void SeparationRho(double _separationRho) { separationRho = _separationRho; }

  double NearRho() const { return nearRho; }

  void NearRho(double _nearRho) { nearRho = _nearRho; }

  double SeparationCluster() const { return separationCluster; }

  void SeparationCluster(double _separationCluster) { separationCluster = _separationCluster; }

  void PrintInfo();
private:
  void outputStep(IAEigenvalueType<VERIFIED> rho, int n) const;

  bool computeNextStepValue(IAEigenvalueAssemble<VERIFIED> &assemble, IEigenSolver &esolver,
                            Eigenpairs &EP, double rho, Matrix &A, Matrix &B, int bisectionSteps);

  bool checkRayleighQuotient(IAEigenvalueAssemble<VERIFIED> &assemble, const Eigenfct &u,
                             double rho, double step, bool enableOutput = true);

  bool check(IAEigenvalueType<VERIFIED> rayleighQuotient, double rho);

  void computeNewRho(IAEigenvalueAssemble<VERIFIED> &assemble, int levelGoerisch, Eigenpairs &EP,
                     double &rho);

  void checkHomotopy(IAEigenvalueAssemble<VERIFIED> &assemble, Eigenpairs &EP, const double &rho);
};

#endif // IAEIGENVALUEMETHODS_H
