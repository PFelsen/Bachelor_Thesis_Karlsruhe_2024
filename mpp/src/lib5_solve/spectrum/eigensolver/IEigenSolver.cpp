#include "IEigenSolver.hpp"

using namespace std::placeholders;

void RandomVector(Vector &u, Operator &B, Operator &IA, Vector &tmp) {
  tmp = u;
  for (int i = 0; i < tmp.size(); ++i)
    tmp[i] = double(rand()) / double(RAND_MAX);
  tmp.ClearDirichletValues();
  Vector b = B * tmp;
  double s = sqrt(b * tmp);
  b /= s;
  u = IA * b;
}

void RandomVectors(Vectors &U, Operator &B, Operator &IA, Vector &tmp) {
  for (int i = 0; i < U.size(); ++i)
    RandomVector(U[i], B, IA, tmp);
}

void RandomVectorP(Vector &u, Operator &B, Operator &IA, Operator &P, Vector &tmp) {
  RandomVector(tmp, B, IA, u);
  u = P * tmp;
}

void RandomVectorsP(Vectors &U, Operator &B, Operator &IA, Operator &P, Vector &tmp) {
  for (int i = 0; i < U.size(); ++i)
    RandomVectorP(U[i], B, IA, P, tmp);
}

double RitzDefect(Vector &u, Vector &residual, double lambda, Operator &A, Operator &B,
                  Vector &tmp) {
  residual = A * u;
  tmp = B * u;
  residual -= lambda * tmp;
  return norm(residual);
}

double RitzDefectB(Vector &u, Vector &residual, double lambda, Operator &A, Operator &B,
                   Operator &IB, Vector &tmp) {
  residual = A * u;
  tmp = B * u;
  residual -= lambda * tmp;
  tmp = IB * residual;
  return sqrt(tmp * residual);
}

void IEigenSolver::RitzStep(Vectors &u, Eigenvalues &lambda, Operator &A, Operator &B, Vector &tmp1,
                            Vector &tmp2) {
  SymRMatrix a(u.size()), b(u.size());
  RMatrix e(u.size());
  RitzStep(u, lambda, A, B, a, b, e, tmp1, tmp2);
}

void IEigenSolver::RitzStep(Vectors &u, Eigenvalues &lambda, Operator &A, Operator &B,
                            SymRMatrix &a, SymRMatrix &b, RMatrix &e, Vector &tmp1, Vector &tmp2) {
  for (int s = 0; s < u.size(); ++s) {
    tmp1 = A * u[s];
    tmp2 = B * u[s];
    for (int r = 0; r <= s; ++r) {
      a(u[r] * tmp1, r, s);
      b(u[r] * tmp2, r, s);
    }
  }

  EVreal(a, b, lambda, e);

  Vectors u_copy(u);
  for (int r = 0; r < u.size(); r++) {
    u[r] = 0;
    for (int s = 0; s < u.size(); ++s)
      u[r] += e[s][r] * u_copy[s];
  }
}

void IEigenSolver::print(const Eigenvalues &lambda, int length) const {
  for (int i = 0; i < std::min(length, outputSize); ++i) {
    double value = lambda[i] - outputShift;
    char buf[128];
    buf[0] = 0;
    if (lambda[i] == 0) mout << "         0";
    else if (abs(value) < 1e-12) mout << "         0";
    else if ((abs(value) > 1000) || (abs(value) < 0.001)) {
      sprintf(buf, "%9.3e", value);
      if (value > 0) mout << " ";
      mout << buf;
    } else {
      sprintf(buf, "%10.5f", value);
      mout << buf;
    }
  }
}

void IEigenSolver::output(const Eigenvalues &lambda, int length, int step, bool final) const {
  if (verbose <= 0) return;
  if (step % printSteps != 0 && !final) return;

  if (verbose == 1) {
    mout << "Step " << step << std::beginl;
  } else if (verbose == 2) {
    if (step < 10) mout << "lambda(" << step << ")=  ";
    else mout << "lambda(" << step << ")= ";
    print(lambda, length);
    mout << std::beginl;
  } else {
    if (step < 10) mout << "lambda(" << step << ")=  ";
    else mout << "lambda(" << step << ")= ";
    print(lambda, length);
    mout << endl;
  }
}

void IEigenSolver::output(const Eigenvalues &lambda, int step, bool final) const {
  output(lambda, lambda.size(), step, final);
}

void IEigenSolver::outputFinish(const Eigenvalues &lambda, int step) const {
  output(lambda, step, true);
  if (verbose == 2) mout << endl;
}

IEigenSolver::IEigenSolver(string name, std::unique_ptr<LinearSolver> &&solver,
                           std::unique_ptr<LinearSolver> &&solver2, int verbose, double eps,
                           int maxstep, int printSteps) :
    name(name), solver(std::move(solver)), solver2(std::move(solver2)), eps(eps), maxstep(maxstep),
    printSteps(printSteps) {
  setVerbose(verbose);
}

void IEigenSolver::operator()(Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B, int fev) {
  (*this)(u, lambda, A, B, true, fev);
}

void IEigenSolver::operator()(Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B, Operator &P,
                              int fev) {
  (*this)(u, lambda, A, B, P, true, fev);
}

void IEigenSolver::operator()(Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B,
                              bool initEigenfcts, int fev) {
  Vector tmp(u[0]);
  if (initEigenfcts) {
    RandomVectors(u, B, (*solver)(A), tmp);
  } else {
    (*solver)(A);
  }
  if (solver2) {
    (*solver2)(B);
    (*this)(u, lambda, A, B, std::bind(RitzDefectB, _1, _2, _3, _4, _5, std::ref(*solver2), _6),
            tmp, fev);
  } else {
    (*this)(u, lambda, A, B,
            std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>(
                RitzDefect),
            tmp, fev);
  }
}

void IEigenSolver::operator()(Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B, Operator &P,
                              bool initEigenfcts, int fev) {
  Vector tmp(u[0]);
  if (initEigenfcts) {
    RandomVectorsP(u, B, (*solver)(A), P, tmp);
  } else {
    (*solver)(A);
  }
  if (solver2) {
    (*solver2)(B);
    (*this)(u, lambda, A, B, P, std::bind(RitzDefectB, _1, _2, _3, _4, _5, std::ref(*solver2), _6),
            tmp, fev);
  } else {
    (*this)(u, lambda, A, B, P,
            std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>(
                RitzDefect),
            tmp, fev);
  }
}

void IEigenSolver::PrintInfo() {
  if (solver2) {
    mout.PrintInfo("EigenSolver " + this->name, this->verbose, PrintInfoEntry("Step max", maxstep),
                   PrintInfoEntry("Epsilon", eps), PrintInfoEntry("LinearSolver", solver->Name()),
                   PrintInfoEntry("LinearSolver2", solver2->Name()));
  } else {
    mout.PrintInfo("EigenSolver " + this->name, this->verbose, PrintInfoEntry("Step max", maxstep),
                   PrintInfoEntry("Epsilon", eps), PrintInfoEntry("LinearSolver", solver->Name()));
  }
}

void IEigenSolver::setVerbose(int verbose) {
  this->verbose = verbose;
  if (verbose <= 3) {
    solver->SetVerbose(-1);
    if (solver2) solver2->SetVerbose(-1);
  }
}

void ILOBPCG::operator()(
    Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B,
    std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
        ritzDefect,
    Vector &tmp, int fev) {
  (*this)(u, lambda, A, B, RandomVector, ritzDefect, tmp, fev);
}

void ILOBPCG::operator()(
    Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B, Operator &P,
    std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
        ritzDefect,
    Vector &tmp, int fev) {
  (*this)(u, lambda, A, B, std::bind(RandomVectorP, _1, _2, _3, P, _4), ritzDefect, tmp, fev);
}

void ILOBPCG::operator()(
    Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B,
    std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
    std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
        ritzDefect,
    Vector &tmp, int fev) {
  mout.StartBlock(Name());
  int R = u.size();
  Vector residual(u[0]); // also used as tmp vector
  std::vector<Vectors> additional{};
  for (int i = 0; i < numAdditional; ++i)
    additional.push_back(Vectors(R, u[0]));
  a = std::make_unique<SymRMatrix>(2 * R);
  b = std::make_unique<SymRMatrix>(2 * R);
  e = std::make_unique<RMatrix>(2 * R);

  firstStep = true;
  int step = 0;
  RitzStep(u, lambda, A, B, residual, tmp);
  while (step < maxstep) {
    output(lambda, R, step);
    ++step;
    evConverged = 0;
    for (int r = 0; r < R; ++r) {
      double defect = ritzDefect(u[r], residual, lambda[r], A, B, tmp);
      update(u, additional, lambda, A, B, residual, randVector, defect, r);
    }
    ritzStepLOBPCG(u, additional, lambda, A, B, residual, tmp);
    if (evConverged == R) break;
  }
  lambda.resize(R);
  outputFinish(lambda, step);
  a = nullptr;
  b = nullptr;
  e = nullptr;
  mout.EndBlock();
}

void ILOBPCG::setMatricesSmall(Eigenfcts &u, Vectors &a1, Operator &A, Operator &B, Vector &tmp1,
                               Vector &tmp2) {
  for (int s = 0, s_shifted = u.size(); s < u.size(); ++s, ++s_shifted) {
    tmp1 = A * u[s];
    tmp2 = B * u[s];
    for (int r = 0; r <= s; ++r) {
      (*a)(u[r] * tmp1, r, s);
      (*b)(u[r] * tmp2, r, s);
    }
    tmp1 = A * a1[s];
    tmp2 = B * a1[s];
    for (int r = 0; r < u.size(); ++r) {
      (*a)(u[r] * tmp1, r, s_shifted);
      (*b)(u[r] * tmp2, r, s_shifted);
    }
    for (int r = 0, r_shifted = u.size(); r <= s; ++r, ++r_shifted) {
      (*a)(a1[r] * tmp1, r_shifted, s_shifted);
      (*b)(a1[r] * tmp2, r_shifted, s_shifted);
    }
  }
}

void ILOBPCG::setMatricesLarge(Eigenfcts &u, Vectors &a1, Vectors &a2, Operator &A, Operator &B,
                               Vector &tmp1, Vector &tmp2) {
  setMatricesSmall(u, a1, A, B, tmp1, tmp2);
  for (int s = 0, s_shifted = 2 * u.size(); s < u.size(); ++s, ++s_shifted) {
    tmp1 = A * a2[s];
    tmp2 = B * a2[s];
    for (int r = 0, r_shifted = u.size(); r < u.size(); ++r, ++r_shifted) {
      (*a)(u[r] * tmp1, s_shifted, r);
      (*b)(u[r] * tmp2, s_shifted, r);
      (*a)(a1[r] * tmp1, s_shifted, r_shifted);
      (*b)(a1[r] * tmp2, s_shifted, r_shifted);
    }
    for (int r = 0, r_shifted = 2 * u.size(); r <= s; ++r, ++r_shifted) {
      (*a)(a2[r] * tmp1, s_shifted, r_shifted);
      (*b)(a2[r] * tmp2, s_shifted, r_shifted);
    }
  }
}

void ILOBPCGSelective::operator()(
    Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B,
    std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
    std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
        ritzDefect,
    Vector &tmp, int fev) {
  mout.StartBlock(Name());
  int R = u.size();
  if (fev == 0) fev = R;
  Vector residual(u[0]); // also used as tmp vector
  std::vector<Vectors> additional{};
  for (int i = 0; i < numAdditional; ++i)
    additional.push_back(Vectors(R, u[0]));
  a = std::make_unique<SymRMatrix>(2 * R);
  b = std::make_unique<SymRMatrix>(2 * R);
  e = std::make_unique<RMatrix>(2 * R);
  selection = new int[R];
  for (int r = 0; r < R; ++r)
    selection[r] = 2;

  firstStep = true;
  int step = 0;
  evConverged = 0;

  RitzStep(u, lambda, A, B, residual, tmp);
  while (step < maxstep) {
    output(lambda, R, step);
    ++step;
    bool finish = update(u, additional, lambda, A, B, residual, tmp, randVector, ritzDefect, fev);
    if (finish) break;
    ritzStepLOBPCG(u, additional, lambda, A, B, residual, tmp);
  }
  lambda.resize(R);
  outputFinish(lambda, step);
  a = nullptr;
  b = nullptr;
  e = nullptr;
  delete[] selection;
  selection = nullptr;
  mout.EndBlock();
}

void ILOBPCGSelective::setMatricesSmallSelective(Eigenfcts &u, Vectors &a1, Operator &A,
                                                 Operator &B, Vector &tmp1, Vector &tmp2) {
  for (int s = 0, s_shifted = u.size(); s < u.size(); ++s, ++s_shifted) {
    tmp1 = A * u[s];
    tmp2 = B * u[s];
    for (int r = 0; r <= s; ++r) {
      if ((selection[r] > 0) || (selection[s] > 0)) {
        (*a)(u[r] * tmp1, r, s);
        (*b)(u[r] * tmp2, r, s);
      }
    }
    tmp1 = A * a1[s];
    tmp2 = B * a1[s];
    for (int r = 0; r < u.size(); ++r) {
      if ((selection[r] > 0) || (selection[s] > 0)) {
        (*a)(u[r] * tmp1, r, s_shifted);
        (*b)(u[r] * tmp2, r, s_shifted);
      }
    }
    for (int r = 0, r_shifted = u.size(); r <= s; ++r, ++r_shifted) {
      if ((selection[r] > 0) || (selection[s] > 0)) {
        (*a)(a1[r] * tmp1, r_shifted, s_shifted);
        (*b)(a1[r] * tmp2, r_shifted, s_shifted);
      }
    }
  }
}

void ILOBPCGSelective::setMatricesLargeSelective(Eigenfcts &u, Vectors &a1, Vectors &a2,
                                                 Operator &A, Operator &B, Vector &tmp1,
                                                 Vector &tmp2) {
  setMatricesSmallSelective(u, a1, A, B, tmp1, tmp2);
  for (int s = 0, s_shifted = 2 * u.size(); s < u.size(); ++s, ++s_shifted) {
    tmp1 = A * a2[s];
    tmp2 = B * a2[s];
    for (int r = 0, r_shifted = u.size(); r < u.size(); ++r, ++r_shifted) {
      if ((selection[r] > 0) || (selection[s] > 0)) {
        (*a)(u[r] * tmp1, s_shifted, r);
        (*b)(u[r] * tmp2, s_shifted, r);
        (*a)(a1[r] * tmp1, s_shifted, r_shifted);
        (*b)(a1[r] * tmp2, s_shifted, r_shifted);
      }
    }
    for (int r = 0, r_shifted = 2 * u.size(); r <= s; ++r, ++r_shifted) {
      if ((selection[r] > 0) || (selection[s] > 0)) {
        (*a)(a2[r] * tmp1, s_shifted, r_shifted);
        (*b)(a2[r] * tmp2, s_shifted, r_shifted);
      }
    }
  }
}