#include "ExponentialIntegrator.hpp"

void ExponentialIntegrator::arnoldiPolynomial(double dt, const Vector &u, const Matrix &massMatrix,
                                              const Matrix &systemMatrix, Vector &du,
                                              int function) {
  std::shared_ptr<LazyVectors> v;
  if (!keepLazy) {
    v = std::make_shared<LazyVectors>(KMax + 1, u);
  } else {
    if (!vInt) vInt = std::make_shared<LazyVectors>(KMax + 1, u);
    v = vInt;
  }
  Vector Av(u);
  RMatrix h(KMax + 1, KMax);
  RVector s(1);
  RVector s_old(1);
  s_old = Scalar(0.0);
  h = Scalar(0.0);
  (*v)[0] = u;
  double beta = sqrt((*v)[0] * (massMatrix * (*v)[0]));
  if (beta < betaEps) {
    if (verbose > 1) mout << "Skipping step because beta is too small " << endl;
    du = 0;
    return;
  }
  (*v)[0] *= 1 / beta;
  s_old[0] = beta;
  RVector s_old_old = s_old;
  int k = 0;
  double eta_rel = infty;
  for (; k < KMax; ++k) {
    polynomialOneStep(systemMatrix, massMatrix, h, Av, *v, k);
    RMatrix H_k(k + 1, k + 1);
    fillMat(h, H_k);
    H_k *= dt;
    evaluateExp(H_k, s, function);
    s *= beta;
    double eta = getEta(s_old, s);
    if (eta == -infty) continue;
    if (k == 1) eta_rel = eta;
    if (eta < Keps || eta / eta_rel < KepsRel) {
      if (verbose > 1) {
        mout << "d(0) eta:" << eta_rel << endl;
        mout << "d(" << H_k.cols() - 1 << ") eta:" << eta << endl;
      }
      break;
    }
  }
  combineSolution(s, *v, du, k, dt);
  if (k == KMax) du[0] = 2e10;
}

void ExponentialIntegrator::polynomialOneStep(const Matrix &systemMatrix, const Matrix &massMatrix,
                                              RMatrix &h, Vector &w, LazyVectors &v, int k) {
  w = systemMatrix * v[k];
  v[k + 1] = (*solvers.at(0)) * w;
  for (int j = 0; j <= k; ++j) {
    h[j][k] = w * v[j];
    v[k + 1] -= h[j][k] * v[j];
  }
  w = massMatrix * v[k + 1];
  h[k + 1][k] = sqrt(v[k + 1] * w);
  v[k + 1] *= (1.0 / h[k + 1][k]);
}

void ExponentialIntegrator::arnoldiRational(double dt, const Vector &u, const Matrix &massMatrix,
                                            const Matrix &systemMatrix, Vector &du, int function) {
  std::shared_ptr<LazyVectors> v;
  std::shared_ptr<LazyVectors> w;
  if (!keepLazy) {
    v = std::make_shared<LazyVectors>(KMax + 1, u);
    w = std::make_shared<LazyVectors>(KMax + 1, u);

  } else {
    if (!vInt) vInt = std::make_shared<LazyVectors>(KMax + 1, u);
    if (!wInt) wInt = std::make_shared<LazyVectors>(KMax + 1, u);
    w = wInt;
    v = vInt;
  }
  RMatrix h(KMax + 1, KMax);
  RVector s(1);
  RVector s_old(1);
  h = Scalar(0.0);
  (*v)[0] = u;
  (*w)[0] = massMatrix * u;
  double beta = sqrt((*v)[0] * (*w)[0]);
  if (beta < betaEps) {
    mout << "Skipping step because beta is too small " << endl;
    du = 0;
    return;
  }
  (*v)[0] *= 1 / beta;
  (*w)[0] *= 1 / beta;
  s_old[0] = beta;
  solvers.at(1)->SetReduction(Keps);
  int k = 0;
  double d = 0;
  double eta_rel = infty;
  for (; k < KMax; ++k) {
    rationalOneStep(systemMatrix, massMatrix, h, *w, *v, k);
    RMatrix S_k(k + 1, k + 1);
    fillMat(h, S_k);
    S_k.Invert();
    S_k *= -1.0;
    for (int j = 0; j <= k; ++j) {
      S_k[j][j] += gamma;
    }
    evaluateExp(S_k, s, function);
    s *= beta;
    double eta = getEta(s_old, s);
    if (eta == -infty) continue;
    if (k == 1) eta_rel = eta;
    solvers.at(1)->SetReduction(Keps / (Keps + eta));
    if (eta < Keps || eta / eta_rel < KepsRel) {
      mout << "d(0) eta:" << eta_rel << endl;
      mout << "d(" << S_k.cols() - 1 << ") eta:" << eta << endl;
      stagesSum += k;
      break;
    }
  }
  combineSolution(s, *v, du, k, dt);
  if (k == KMax) du = 2e10;
}

void ExponentialIntegrator::rationalOneStep(const Matrix &systemMatrix, const Matrix &massMatrix,
                                            RMatrix &h, LazyVectors &w, LazyVectors &v, int k) {
  v[k + 1] = (*solvers.at(1)) * w[k];
  for (int j = 0; j <= k; ++j) {
    h[j][k] = v[k + 1] * w[j];
    v[k + 1] -= h[j][k] * v[j];
  }
  w[k + 1] = massMatrix * v[k + 1];
  h[k + 1][k] = sqrt(v[k + 1] * w[k + 1]);
  v[k + 1] *= (1.0 / h[k + 1][k]);
  w[k + 1] *= (1.0 / h[k + 1][k]);
}

void ExponentialIntegrator::fillMat(const RMatrix &h, RMatrix &H_k) {
  for (int j = 0; j < H_k.cols(); ++j) {
    for (int l = 0; l < H_k.rows(); ++l) {
      H_k[j][l] = h[j][l];
    }
  }
}

double ExponentialIntegrator::getEta(RVector &s_old, RVector &s) {
  double norm_u = s.norm();
  double delta = (s_old - s).norm() / norm_u;
  s_old = s;
  s_old.push_back(0.0);
  double eta = -infty;
  if (s.size() == 1) return eta;
  eta = 1 + norm_u;
  if (delta < 1) eta = min(eta, delta * norm_u / (1 - delta));
  return eta;
}

void ExponentialIntegrator::evaluateExp(RMatrix &H_k, RVector &s, int function) {
  if (function == 0) {
    s = RVector(0.0, H_k.cols());
    s[0] = 1.0;
    H_k.Exp(s, true);
  }
  if (function == 1) {
    RMatrix H_k_k(0.0, H_k.cols() + 1);
    H_k_k.Insert(H_k, 0, 0);
    H_k_k[0][H_k.cols()] = 1.0;
    s = RVector(0.0, H_k_k.cols());
    s[H_k.cols()] = 1.0;
    H_k_k.Exp(s, true);
    s.removeLast();
  }
  if (function == 2) {
    RMatrix H_k_k_k(0.0, H_k.cols() + 2);
    H_k_k_k.Insert(H_k, 0, 0);
    H_k_k_k[0][H_k_k_k.cols() - 2] = 1.0;
    H_k_k_k[0][H_k_k_k.cols() - 1] = 1.0;
    H_k_k_k[H_k_k_k.cols() - 2][H_k_k_k.cols() - 1] = 1.0;
    H_k_k_k.Exp();
    RVector Vec1 = H_k_k_k.col(H_k_k_k.cols() - 2);
    s = H_k_k_k.col(H_k_k_k.cols() - 1) - H_k_k_k.col(H_k_k_k.cols() - 2);
    s.resize(H_k.cols());
  }
  if (function == 3) {
    Exit("optimize this with the formula from Higham 2011 evaluating the exponential like above...")
        H_k.Phi3();
  }
}

void ExponentialIntegrator::combineSolution(const RVector &s, LazyVectors &v, Vector &du, int k,
                                            double dt) {
  du = 0.0;
  for (int j = 0; j <= k; ++j)
    du += dt * s[j] * v[j];
  //    du *= dt;
}

void ExponentialIntegrator::assembleRationalMatrix(double dt, const Matrix &MassMatrix,
                                                   const Matrix &SystemMatrix) {
  solvers.at(0)->operator()(MassMatrix, 1);
  stageMatrix = std::make_unique<Matrix>(MassMatrix);
  *stageMatrix *= gamma;
  *stageMatrix += (-1.0) * dt * SystemMatrix;
  stageMatrix->CreateSparse();
  solvers.at(1)->operator()(*stageMatrix);
  assembled = true;
}

void ExponentialMidpointWithShift::StageFunctions(double dt, const Vector &u,
                                                  const Matrix &massMatrix,
                                                  const Matrix &systemMatrix, Vector &du) {
  if (!assembled) {
    assembleRationalMatrix(dt, massMatrix, systemMatrix);
    dtOld = dt;
  } else if (abs(dtOld - dt) > Eps) Exit("only for uniform timegrids");
  rhs[0] += systemMatrix * u;
  du = (*solvers.at(0)) * rhs[0];
  Vector tmp = du;

  arnoldiRational(dt, tmp, massMatrix, systemMatrix, du, 1);
}

void ExponentialMidpointWithShift::AssembleRHS(ILinearTimeAssemble *assemble, const Vector &u) {
  assemble->RHS(assemble->Time() - 0.5 * assemble->StepSize(), rhs[0]);
}

void ExponentialMidpoint::StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                                         const Matrix &systemMatrix, Vector &du) {
  if (!assembled) {
    solvers.at(0)->operator()(massMatrix, 1);
    assembled = true;
  }
  rhs[0] += systemMatrix * u;
  du = (*solvers.at(0)) * rhs[0];
  Vector tmp = du;
  arnoldiPolynomial(dt, tmp, massMatrix, systemMatrix, du, 1);
}

void ExponentialMidpoint::AssembleRHS(ILinearTimeAssemble *assemble, const Vector &u) {
  assemble->RHS(assemble->Time() - 0.5 * assemble->StepSize(), rhs[0]);
}

void ExponentialEuler::StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                                      const Matrix &systemMatrix, Vector &du) {
  if (!assembled) {
    solvers.at(0)->operator()(massMatrix, 1);
    assembled = true;
  }
  rhs[0] += systemMatrix * u;
  du = (*solvers.at(0)) * rhs[0];
  arnoldiPolynomial(dt, du, massMatrix, systemMatrix, du, 1);
}

void ExponentialEuler::AssembleRHS(ILinearTimeAssemble *assemble, const Vector &u) {
  assemble->RHS(assemble->Time() - assemble->StepSize(), rhs[0]);
}

void ExponentialEulerWithShift::StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                                               const Matrix &systemMatrix, Vector &du) {
  if (!assembled) {
    assembleRationalMatrix(dt, massMatrix, systemMatrix);
    dtOld = dt;
    assembled = true;
  } else if (abs(dtOld - dt) > Eps) Exit("only for uniform timegrids") rhs[0] += systemMatrix * u;
  du = (*solvers.at(0)) * rhs[0];
  arnoldiRational(dt, du, massMatrix, systemMatrix, du, 1);
}
