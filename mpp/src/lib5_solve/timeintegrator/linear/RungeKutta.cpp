#include "RungeKutta.hpp"

int RungeKutta::GetOrder() const {
  CheckButcherTableau();
  double tmp = 0.0;
  int order = 0;
  for (int i = 0; i < b.size(); ++i)
    tmp += b[i];
  if ((tmp - 1.0) > Eps) {
    return order;
  } else {
    order = 1;
    tmp = 0.0;
  }
  for (int i = 0; i < b.size(); ++i)
    tmp += b[i] * c[i];
  if (abs(tmp - 0.5) > Eps) {
    return order;
  } else {
    order = 2;
    tmp = 0.0;
  }
  double tmp2 = 0.0;
  for (int i = 0; i < b.size(); ++i) {
    tmp += b[i] * c[i] * c[i];
    for (int j = 0; j < b.size(); ++j) {
      tmp2 += b[i] * A[i][j] * c[j];
    }
  }
  if (((abs(tmp - 1.0 / 3.0)) > Eps) || (abs(tmp2 - 1.0 / 6.0)) > Eps) {
    return order;
  } else {
    order = 3;
    tmp = 0.0;
    tmp2 = 0.0;
  }
  double tmp3 = 0.0;
  double tmp4 = 0.0;
  for (int i = 0; i < b.size(); ++i) {
    tmp += b[i] * c[i] * c[i] * c[i];
    for (int j = 0; j < b.size(); ++j) {
      tmp2 += b[i] * c[i] * A[i][j] * c[j];
      tmp3 += b[i] * c[j] * A[i][j] * c[j];
      for (int k = 0; k < b.size(); ++k) {
        tmp4 += b[i] * A[i][j] * A[j][k] * c[k];
      }
    }
  }
  if (((abs(tmp - 1.0 / 4.0)) > Eps) || (abs(tmp2 - 1.0 / 8.0)) > Eps
      || (abs(tmp3 - 1.0 / 12.0)) > Eps || (abs(tmp4 - 1.0 / 24.0)) > Eps) {
    return order;
  } else {
    order = 4;
    return order;
  }
}

void RungeKutta::CheckButcherTableau() const {
  if ((c.size() != b.size()) || (A.cols() != A.rows()) || (c.size() != A.cols())) {
    mout << "c.size: " << c.size() << " b.size(): " << b.size() << " A.cols(): " << A.cols()
         << "A.rows(): " << A.rows() << endl;
    THROW("Butcher tableau invalid")
  }
}

void CrankNicolsonHardcoded::StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                                            const Matrix &systemMatrix, Vector &du) {
  if (!assembled) {
    stageMatrix = std::make_unique<Matrix>(massMatrix);
    *stageMatrix += (-0.5) * dt * systemMatrix;
    stageMatrix->CreateSparse();
    solvers.at(0)->operator()(*stageMatrix, 1);
    assembled = true;
  }
  Vector systemMatrixu(systemMatrix * u);
  systemMatrixu += *RHS;
  du *= (1 / dt);
  solvers.at(0)->SolveWithInitialValue(du, systemMatrixu);
  du *= dt;
}

void CrankNicolsonHardcoded::AssembleRHS(ILinearTimeAssemble *assemble, const Vector &u) {
  if (!rhsInitialized) {
    RHS = std::make_unique<Vector>(u);
    rhsInitialized = true;
  }
  Vector tmp(u);
  (*assemble).RHS(assemble->Time(), tmp);
  (*assemble).RHS(assemble->Time() - assemble->StepSize(), *RHS);
  *RHS += tmp;
  *RHS *= 0.5;
}

void ExplicitRungeKutta::StageFunctions(double dt, const Vector &u, const Matrix &massMatrix,
                                        const Matrix &systemMatrix, Vector &du) {
  mout.StartBlock("Stage Function");
  int stages = int(b.size());
  vector<Vector> k(stages, u);
  Vector systemMatrixK(u);
  du = 0.0;
  for (int stage = 0; stage < stages; ++stage) {
    for (int j = 0; j < stage; ++j) {
      k[stage] += dt * A[stage][j] * k[j];
    }
    systemMatrixK = systemMatrix * k[stage] + rhs[stage];
    ;
    k[stage] = (*solvers.at(0)) * systemMatrixK;
    du += b[stage] * dt * k[stage];
  }
  stagesSum += 1;
  mout.EndBlock(verbose <= 1);
}

void DiagonalImplicitRungeKutta::StageFunctions(double dt, const Vector &u,
                                                const Matrix &massMatrix,
                                                const Matrix &systemMatrix, Vector &du) {
  mout.StartBlock("Stage Function");
  int stages = int(b.size());
  vector<Vector> k(stages, u);
  Vector systemMatrixK(u);
  Vector du0(du);
  du = 0.0;
  for (int stage = 0; stage < stages; ++stage) {
    for (int j = 0; j < stage; ++j) {
      k[stage] += dt * A[stage][j] * k[j];
    }
    systemMatrixK = systemMatrix * k[stage] + rhs[stage];
    std::shared_ptr<Matrix> B = getStageMatrix(dt, stage);
    k[stage] = 1.0 / dt * du0;
    solvers.at(stage)->SolveWithInitialValue(k[stage], systemMatrixK);
    if (stage > 3)
      Exit("diagonal implicit runge kutta only for maximal 3 stages!") vout(2)
          << "k stage " << k[stage].norm() << endl;
    du += dt * b[stage] * k[stage];
  }
  vout(2) << "du stage " << du.norm() << endl;
  stagesSum += stages;
  mout.EndBlock(verbose <= 1);
}

void DiagonalImplicitRungeKutta::initializeStageMatrices(double dt) {
  int stages = int(b.size());
  //    stageMatrix =std::vector(stages,std::make_shared<Matrix>(*massMatrix));
  for (int stage = 0; stage < stages; ++stage) {
    stageMatrix.push_back(std::make_shared<Matrix>(*massMatrix));
    double factor = A[stage][stage] * dt;
    *(stageMatrix[stage]) += (-factor) * (*systemMatrix);
    (stageMatrix[stage])->CreateSparse();
    solvers.at(stage)->operator()((*stageMatrix[stage]), 1);
  }
  assembled = true;
}

std::shared_ptr<Matrix> DiagonalImplicitRungeKutta::getStageMatrix(double dt, int stage) {
  if (!assembled) {
    initializeStageMatrices(dt);
  } else {
    if (abs(dt - dtOld) > Eps) {
      MOUT(dt);
      MOUT(dt - dtOld);
      mout << "Change of Stepsize detected! Assembling again " << endl;
      initializeStageMatrices(dt);
    }
  }
  dtOld = dt;
  return stageMatrix[stage];
}
