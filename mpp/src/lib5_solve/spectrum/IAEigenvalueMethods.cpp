#include "IAEigenvalueMethods.hpp"
#include "IAIO.hpp"

template<bool VERIFIED>
IAEigenvalueMethod<VERIFIED>::IAEigenvalueMethod(const std::string &name) :
    name(name), verbose(1), minESSize(5), addEVs(3) {
  Config::Get(name + "Verbose", verbose);
  std::string esolver = "LOBPCGSelectiveVeryFast";
  Config::Get(name + "ESolver", esolver);
  ES = std::unique_ptr<IEigenSolver>(EigenSolverCreator(esolver, name).Create());
  esolver = "LOBPCGSelectiveVeryFast";
  Config::Get(name + "CoarseESolver", esolver);
  ES_coarse = std::unique_ptr<IEigenSolver>(EigenSolverCreator(esolver, name + "Coarse").Create());
  Config::Get(name + "MinimalEVSize", minESSize);
}

template<bool VERIFIED>
void IAEigenvalueMethod<VERIFIED>::computeEigenpairs(IAEigenvalueAssemble<VERIFIED> &assemble,
                                                     IEigenSolver &esolver, Eigenfcts &U,
                                                     Eigenvalues &lambda, Matrix &A, Matrix &B,
                                                     const double &t, bool initBC) {
  vout(0) << "Approximative eigenvalues on level " << U.SpaceLevel() << " (t= " << t << ")" << endl;

  int n = U.size();
  U.resizeData(std::max(minESSize, n + addEVs));

  if (initBC) assemble.BoundaryConditionsEigenSolver(U);
  assemble.MatricesEigenSolver(A, B, t);

  esolver.SetOutputShift(assemble.Shift());
  esolver.SetOutputSize(n);
  esolver(U, lambda, A, B, true);

  U.resizeData(n);
  lambda.resize(n);
  lambda -= assemble.Shift();
}

template<bool VERIFIED>
void IAEigenvalueMethod<VERIFIED>::computeEigenpairs(IAEigenvalueAssemble<VERIFIED> &assemble,
                                                     IEigenSolver &esolver, Eigenpairs &EP,
                                                     Matrix &A, Matrix &B, const double &t,
                                                     bool initBC) {
  computeEigenpairs(assemble, esolver, EP.U, EP.lambda, A, B, t, initBC);
}

template<bool VERIFIED>
void IAEigenvalueMethod<VERIFIED>::computeEigenpairs(IAEigenvalueAssemble<VERIFIED> &assemble,
                                                     IEigenSolver &esolver, Eigenpairs &EP,
                                                     const double &t) {
  Matrix A(EP(0));
  Matrix B(EP(0));
  computeEigenpairs(assemble, esolver, EP.U, EP.lambda, A, B, t, true);
}

template<bool VERIFIED>
void IAEigenvalueMethod<VERIFIED>::computeEigenfcts(IAEigenvalueAssemble<VERIFIED> &assemble,
                                                    IEigenSolver &esolver, Eigenfcts &U,
                                                    const double &t) {
  Matrix A(U[0]);
  Matrix B(U[0]);
  Eigenvalues lambda;
  computeEigenpairs(assemble, esolver, U, lambda, A, B, t, true);
}

template class IAEigenvalueMethod<true>;

template class IAEigenvalueMethod<false>;

template<bool VERIFIED>
void IARayleighRitzMethod<VERIFIED>::operator()(IAEigenvalueAssemble<VERIFIED> &assemble,
                                                Eigenfcts &U, IAUpperEigenvalueBounds &Lambda,
                                                bool initEF, const double &t) {
  if (initEF) this->computeEigenfcts(assemble, *this->ES, U, t);

  mout.StartBlock("Rayleigh Ritz");

  int R = U.size();
  if (Lambda.size() != R) Lambda.resize(R);
  if (this->verbose >= 1) mout << "size= " << R << std::endl;

  SymRMatrixT<IAEigenvalueType<VERIFIED>> a(R);
  SymRMatrixT<IAEigenvalueType<VERIFIED>> b(R);
  assemble.MatricesRayleighRitz(U, a, b, t);

  if constexpr (VERIFIED) verifiedEVreal(a, b, Lambda);
  else EVreal(a, b, Lambda);

  Lambda -= assemble.Shift();
  mout.EndBlock();
}

template class IARayleighRitzMethod<true>;

template class IARayleighRitzMethod<false>;

template<bool VERIFIED>
IALehmannGoerischMethod<VERIFIED>::IALehmannGoerischMethod() :
    IAEigenvalueMethod<VERIFIED>("LehmannGoerisch") {
  std::string linearSolver = "GMRES";
  Config::Get(this->name + "LinearSolver", linearSolver);
  std::string preconditioner = "PS";
  Config::Get(this->name + "Preconditioner", preconditioner);
  solver = GetLinearSolverUnique(linearSolver, GetPC(preconditioner), this->name + "Linear");
}

template<bool VERIFIED>
void IALehmannGoerischMethod<VERIFIED>::operator()(IAEigenvalueAssemble<VERIFIED> &assemble,
                                                   int levelGoerisch, Eigenpairs &EP,
                                                   IALowerEigenvalueBounds &lambda,
                                                   const double &rho, bool initEP,
                                                   const double &t) {
  mout.StartBlock("Lehmann Goerisch");

  if (initEP) this->computeEigenpairs(assemble, *this->ES, EP, t);

  int R = EP.size();
  if (lambda.size() != R) lambda.resize(R);
  if (this->verbose >= 1) mout << "size= " << R << std::endl;

  IAEigenvalueType<VERIFIED> shiftedRho = IAEigenvalueType<VERIFIED>(rho) + assemble.Shift();
  std::shared_ptr<IDiscretization> discGoerisch =
      std::move(assemble.CreateGoerischDisc(EP.GetDisc()));
  Goerischfcts W(R, 0.0, discGoerisch, levelGoerisch);

  std::shared_ptr<IAIDiscretization> iaDiscGoerisch = nullptr;
  if constexpr (VERIFIED) {
    iaDiscGoerisch = std::move(assemble.CreateIAGoerischDisc(EP.GetIADisc()));
    W.SetIADisc(iaDiscGoerisch);
  }
  computeGoerischFcts(assemble, EP, W, t);

  SymRMatrixT<IAEigenvalueType<VERIFIED>> a(R);
  SymRMatrixT<IAEigenvalueType<VERIFIED>> b(R);
  SymRMatrixT<IAEigenvalueType<VERIFIED>> g(R);
  assemble.MatricesRayleighRitz(EP.getEigenfcts(), a, b, t);
  assemble.MatrixGoerisch(EP, W, g, t);

  if (R == 1) {
    if constexpr (VERIFIED) {
      lambda[0] = inf((shiftedRho * b(0, 0) - a(0, 0)) / (shiftedRho * g(0, 0) - b(0, 0))
                      - assemble.Shift());
    } else {
      lambda[0] =
          (shiftedRho * b(0, 0) - a(0, 0)) / (shiftedRho * g(0, 0) - b(0, 0)) - assemble.Shift();
    }

    if (this->verbose >= 4)
      mout << "M_t(u,u)= " << a(0, 0) << endl
           << "N_t(u,u)= " << b(0, 0) << endl
           << "b_t(w,w)= " << g(0, 0) << endl;
    if (this->verbose >= 6)
      mout << "rho_old=  " << IAEigenvalueType<VERIFIED>(rho) << endl
           << "rho_new=  " << IAEigenvalueType<VERIFIED>(lambda[0]) << endl;
  } else if (R > 1) {
    SymRMatrixT<IAEigenvalueType<VERIFIED>> aa(b);
    aa *= (-shiftedRho);
    aa += a;
    SymRMatrixT<IAEigenvalueType<VERIFIED>> bb(g);
    if constexpr (VERIFIED) bb *= sqr(shiftedRho);
    else bb *= shiftedRho * shiftedRho;
    b *= (-2 * shiftedRho);
    bb += b;
    bb += a;

    if constexpr (VERIFIED) {
      IAEigenvalueEnclosures l(R);
      verifiedEVreal(aa, bb, l);

      for (int i = 0; i < R; ++i)
        lambda[i] = inf(-(l[R - i - 1] * rho + assemble.Shift()) / (1 - l[R - i - 1]));
    } else {
      Eigenvalues l(R);
      EVreal(aa, bb, l);

      for (int i = 0; i < R; ++i)
        lambda[i] = -(l[R - i - 1] * rho + assemble.Shift()) / (1 - l[R - i - 1]);
    }

    if (this->verbose >= 6)
      mout << "rho_old=  " << IAEigenvalueType<VERIFIED>(rho) << endl
           << "rho_new=  " << IAEigenvalueType<VERIFIED>(lambda[0]) << endl;
  } else {
    ERROR("Size does not fit!")
  }

  mout.EndBlock();
}

template<bool VERIFIED>
void IALehmannGoerischMethod<VERIFIED>::computeGoerischFcts(
    IAEigenvalueAssemble<VERIFIED> &assemble, const Eigenpairs &ep, Goerischfcts &w, double t) {
  assemble.SetGoerischData(ep, t);
  (*solver)(assemble, w);
  assemble.ClearGoerischData();
}

template class IALehmannGoerischMethod<true>;

template class IALehmannGoerischMethod<false>;

template<bool VERIFIED>
void IAHomotopyMethod<VERIFIED>::readConfig(bool all) {
  if (all) {
    Config::Get(this->name + "StepMin", stepMin);
    Config::Get(this->name + "StepMax", stepMax);
    Config::Get(this->name + "StepSize", stepSize);
  }
  Config::Get(this->name + "SeparationRho", separationRho);
  Config::Get(this->name + "NearRho", nearRho);
  Config::Get(this->name + "SeparationCluster", separationCluster);
  Config::Get(this->name + "NumBisectionSteps", numBisectionSteps);
  Config::Get(this->name + "NumBisectionStepsCoarse", numBisectionSteps_coarse);
  if (separationRho > nearRho) nearRho = 3 * separationRho;
  Config::Get(this->name + "NearZero", nearZero);
}

template<bool VERIFIED>
IAHomotopyMethod<VERIFIED>::IAHomotopyMethod(double stepMin, double stepMax, double stepSize,
                                             int coarseLevel, const std::string &name) :
    IAEigenvalueMethod<VERIFIED>(name + "Homotopy"), stepMin(stepMin), stepMax(stepMax),
    stepSize(stepSize), step(stepMin), separationRho(1e-2), nearRho(3e-2), separationCluster(1e-2),
    numBisectionSteps(5), numBisectionSteps_coarse(10), coarseLevel(coarseLevel) {
  readConfig(false);
}

template<bool VERIFIED>
IAHomotopyMethod<VERIFIED>::IAHomotopyMethod(int coarseLevel, const std::string &name) :
    IAEigenvalueMethod<VERIFIED>(name + "Homotopy"), stepMin(0.0), stepMax(1.0), stepSize(0.05),
    step(0.0), separationRho(1e-2), nearRho(3e-2), separationCluster(1e-2), numBisectionSteps(5),
    numBisectionSteps_coarse(10), coarseLevel(coarseLevel) {
  readConfig(true);
}

template<bool VERIFIED>
IAHomotopyMethod<VERIFIED>::IAHomotopyMethod(double stepMin, double stepMax, double stepSize,
                                             const std::string &name) :
    IAHomotopyMethod<VERIFIED>(stepMin, stepMax, stepSize, -1, name) {}

template<bool VERIFIED>
IAHomotopyMethod<VERIFIED>::IAHomotopyMethod(const std::string &name) :
    IAHomotopyMethod<VERIFIED>(-1, name) {}

std::string int2count(int i) {
  switch (i) {
  case 1:
    return "1st";
  case 2:
    return "2nd";
  case 3:
    return "3rd";
  default:
    return std::to_string(i) + "th";
  }
}

template<bool VERIFIED>
double IAHomotopyMethod<VERIFIED>::operator()(IAEigenvalueAssemble<VERIFIED> &assemble,
                                              int levelGoerisch, Eigenpairs &EP, double &rho,
                                              bool initEP) {
  assemble.PrintInfo();
  PrintInfo();

  mout.StartBlock(this->name);

  step = stepMin;
  outputStep(rho, EP.size());
  int size_start = EP.size();

  std::unique_ptr<Matrix> A = nullptr;
  std::unique_ptr<Matrix> B = nullptr;
  std::unique_ptr<Eigenpairs> EP_coarse = nullptr;

  bool doOnCoarseLevel = (coarseLevel >= 0) && (EP(0).SpaceLevel() != coarseLevel);

  if (doOnCoarseLevel) {
    EP_coarse = std::make_unique<Eigenpairs>(EP.size(), EP.GetSharedDisc(), coarseLevel);
    assemble.BoundaryConditionsEigenSolver(EP_coarse->U);
  }
  assemble.BoundaryConditionsEigenSolver(EP.U);


  while (step < stepMax) {
    double previousStep = step;
    if (doOnCoarseLevel) {
      A = std::make_unique<Matrix>((*EP_coarse)(0));
      B = std::make_unique<Matrix>((*EP_coarse)(0));
      this->computeEigenpairs(assemble, *this->ES_coarse, *EP_coarse, *A, *B, step, initEP);
      if (computeNextStepValue(assemble, *this->ES_coarse, *EP_coarse, rho, *A, *B,
                               numBisectionSteps_coarse))
        break;

      A = nullptr;
      B = nullptr;
    }
    A = std::make_unique<Matrix>(EP(0));
    B = std::make_unique<Matrix>(EP(0));
    this->computeEigenpairs(assemble, *this->ES, EP, *A, *B, step, initEP);

    if (rho - separationRho < EP.LastEigenvalue()) {
      step = std::max(previousStep, step - 3 * stepSize);
      this->computeEigenpairs(assemble, *this->ES, EP, *A, *B, step, initEP);
      if (computeNextStepValue(assemble, *this->ES, EP, rho, *A, *B, numBisectionSteps)) break;
      if (rho - separationRho < EP.LastEigenvalue()) { THROW("Eigenvalue still too large") }
    } else if (EP.LastEigenvalue() < rho - nearRho) {
      if (computeNextStepValue(assemble, *this->ES, EP, rho, *A, *B, numBisectionSteps)) break;
    }
    A = nullptr;
    B = nullptr;

    checkRayleighQuotient(assemble, EP.LastEigenfct(), rho, step);

    computeNewRho(assemble, levelGoerisch, EP, rho);

    outputStep(rho, EP.size());

    if (EP.size() == 0) break;
    if (nearZero >= 0.0 && EP[0] > nearZero) break;

    if (doOnCoarseLevel) EP_coarse->resize(EP.size());
  }
  checkHomotopy(assemble, EP, rho);

  if (this->verbose >= 1)
    mout << "lower bound for the " << int2count(EP.size() + 1) << " eigenvalue" << endl;
  mout << "rho(" << step << ")= " << rho << endl
       << "exactly " << (size_start - EP.size()) << " eigenvalue" << ((EP.size() != 1) ? "s" : "")
       << " lost" << endl;

  mout.EndBlock();
  return step;
}

template<bool VERIFIED>
void IAHomotopyMethod<VERIFIED>::PrintInfo() {
  mout.PrintInfo(this->name, this->verbose, PrintInfoEntry("Verified", VERIFIED ? "yes" : "no"),
                 PrintInfoEntry("EigenSolver", this->ES->Name()),
                 PrintInfoEntry("EigenSolver step max", this->ES->MaxStep()),
                 PrintInfoEntry("EigenSolver epsilon", this->ES->Epsilon()),
                 PrintInfoEntry("Bisections steps", numBisectionSteps),
                 PrintInfoEntry("Coarse EigenSolver", this->ES_coarse->Name()),
                 PrintInfoEntry("Coarse EigenSolver step max", this->ES_coarse->MaxStep()),
                 PrintInfoEntry("Coarse EigenSolver epsilon", this->ES_coarse->Epsilon()),
                 PrintInfoEntry("Coarse Bisections steps", numBisectionSteps_coarse),
                 PrintInfoEntry("EigenSolver minimal size", this->minESSize),
                 PrintInfoEntry("EigenSolver additional eigenvalues", this->addEVs),
                 PrintInfoEntry("Step min", stepMin), PrintInfoEntry("Step max", stepMax),
                 PrintInfoEntry("Step size", stepSize),
                 PrintInfoEntry("Separation rho", separationRho),
                 PrintInfoEntry("Near rho", nearRho),
                 PrintInfoEntry("Separation cluster", separationCluster),
                 PrintInfoEntry("Near zero", nearZero, nearZero >= 0.0 ? 1 : 999),
                 PrintInfoEntry("Coarse level", coarseLevel, coarseLevel >= 0 ? 1 : 999));
}

template<bool VERIFIED>
void IAHomotopyMethod<VERIFIED>::outputStep(IAEigenvalueType<VERIFIED> rho, int n) const {
  if (this->verbose >= 1) mout << "rho(" << step << ")= " << rho << endl;
  if (this->verbose >= 4) mout << "n(" << step << ")=   " << n << endl;
}

template<bool VERIFIED>
bool IAHomotopyMethod<VERIFIED>::computeNextStepValue(IAEigenvalueAssemble<VERIFIED> &assemble,
                                                      IEigenSolver &esolver, Eigenpairs &EP,
                                                      double rho, Matrix &A, Matrix &B,
                                                      int bisectionSteps) {
  if (EP.size() == 0) return true;

  double stepLeft;
  double stepRight = step;
  Eigenpairs EP_left(EP);
  for (int k = 0;; ++k) {
    stepLeft = stepRight;
    stepRight = min(stepLeft + pow(2.0, k) * stepSize, stepMax);
    EP_left = EP;
    this->computeEigenpairs(assemble, esolver, EP, A, B, stepRight, false);
    if (rho - nearRho <= EP.LastEigenvalue()) break;
    if (stepRight == stepMax) {
      step = stepRight;
      return true;
    }
  }
  if (EP.LastEigenvalue() <= rho - separationRho) {
    step = stepRight;
    return false;
  }
  for (int k = 0; k < bisectionSteps; ++k) {
    step = (stepLeft + stepRight) / 2;
    this->computeEigenpairs(assemble, esolver, EP, A, B, step, false);
    if (rho - separationRho < EP.LastEigenvalue()) {
      stepRight = step;
    } else if (EP.LastEigenvalue() < rho - nearRho) {
      stepLeft = step;
      EP_left = EP;
    } else {
      return false;
    }
  }
  step = stepLeft;
  EP = EP_left;
  return false;
}

template<bool VERIFIED>
bool IAHomotopyMethod<VERIFIED>::checkRayleighQuotient(IAEigenvalueAssemble<VERIFIED> &assemble,
                                                       const Eigenfct &u, double rho, double step,
                                                       bool enableOutput) {
  IAEigenvalueType<VERIFIED> a, b;
  assemble.MatricesRayleighRitz(u, a, b, step);
  IAEigenvalueType<VERIFIED> rayleighQuotient = a / b - assemble.Shift();

  if (check(rayleighQuotient, rho)) {
    if (enableOutput) {
      if (this->verbose >= 6)
        mout << "a(" << step << ")= " << a << endl << "b(" << step << ")= " << b << endl;
      if (this->verbose >= 4)
        mout << "RayleighQuotient(" << step << ")= " << rayleighQuotient << endl;
    }
    return true;
  }
  Warning("Rayleigh quotient condition not satisfied! \n quot= " + to_string(rayleighQuotient)
          + ";  rho= " + std::to_string(rho));
  return false;
}

template<>
bool IAHomotopyMethod<true>::check(IAInterval rayleighQuotient, double rho) {
  return sup(rayleighQuotient) < rho;
}

template<>
bool IAHomotopyMethod<false>::check(double rayleighQuotient, double rho) {
  return rayleighQuotient < rho;
}

template<bool VERIFIED>
void IAHomotopyMethod<VERIFIED>::computeNewRho(IAEigenvalueAssemble<VERIFIED> &assemble,
                                               int levelGoerisch, Eigenpairs &EP, double &rho) {
  int n = EP.size();
  int clusterSize = 1;
  while (EP[n - 1] - EP[n - clusterSize - 1] < separationCluster) {
    ++clusterSize;
    if (clusterSize == n) break;
  }

  IALowerEigenvalueBounds lambda;
  while (true) {
    Eigenpairs GEP(EP, n - clusterSize);
    GM(assemble, levelGoerisch, GEP, lambda, rho, false, step);

    if (clusterSize == n) break;
    /// Checks if next eigenvalue is still below new rho, otherwise cluster size is increased and
    /// Goerisch method is compute again with larger cluster of eigenvalues
    if (checkRayleighQuotient(assemble, EP(n - clusterSize - 1), lambda[0], step, false)) break;
    clusterSize++;
  }

  rho = lambda[0];
  n -= clusterSize;
  EP.resize(n);
}

template<bool VERIFIED>
void IAHomotopyMethod<VERIFIED>::checkHomotopy(IAEigenvalueAssemble<VERIFIED> &assemble,
                                               Eigenpairs &EP, const double &rho) {
  if (EP.size() > 0) {
    IAUpperEigenvalueBounds Lambda;
    RM(assemble, EP.U, Lambda, false, step);
    if (Lambda.back() > rho) ERROR("Homotopy failed")
  }
  mout << "exactly " << EP.size() << " eigenvalue" << ((EP.size() != 1) ? "s" : "")
       << " below rho= " << rho << endl;
}

template class IAHomotopyMethod<true>;

template class IAHomotopyMethod<false>;
