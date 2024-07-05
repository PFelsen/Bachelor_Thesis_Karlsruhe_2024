#include "Newton.hpp"

void ENewton::operator()(const IAssemble &A, Vector &u) {
  mout.StartBlock("ENewton");
  double E, E_0;
  A.Initialize(u);
  Vector r(u);
  Matrix J(u);
  d = d_0 = A.Residual(u, r);
  E = E_0 = A.Energy(u);
  double eps = defaultEpsilon + defaultReduction * d;
  int LS_cnt = 0;
  int JU_cnt = 0;
  Vector c(r);
  for (iter = 0; iter < max_iter; ++iter) {
    vout(3) << "r(" << iter << ")= " << r << endl;
    if (d < eps) { break; }
    if (d > VeryLarge) { break; }
    if (!(E < VeryLarge)) {
      if (!(E > VeryLarge)) {
        iter = max_iter;
        break;
      }
    }
    vout(1) << "d(" << iter << ")= " << d << "\t E(" << iter << ")= " << E << endl;
    A.Jacobi(u, J);
    JU_cnt = iter;
    c = (*solver)(J, A) * r;
    //    double rc = std::real(r * c);
    double rc = r * c;
    vout(3) << "c(" << iter << ")= " << c << endl;
    u -= c;
    vout(5) << "u-c " << u << endl;
    double d_old = d;
    d = A.Residual(u, r);
    double E_old = E;
    E = A.Energy(u);
    if (d < d_old) continue;
    if (E > E_old + sigma * rc) {
      for (int l = 1; l <= LS_iter; ++l) {
        if (iter == 0)
          if (suppressLS) {
            vout(2) << "line search suppressed" << endl;
            break;
          }
        vout(1) << "line search " << l << ": d(" << iter << ")= " << d << "\t E(" << iter
                << ")= " << E << endl;
        c *= 0.5;
        rc *= 0.5;
        u += c;
        d = A.Residual(u, r);
        E = A.Energy(u);
        if (E < E_old + sigma * rc) break;
      }
    }
    if (E > E_old + sigma * rc) {
      vout(5) << "line search unsuccessful." << endl;
      ++LS_cnt;
      if (LS_cnt == 3) {
        vout(1) << "too many line searches unsuccessful." << endl;
        iter = max_iter;
      }
    }
  }
  vout(0) << "d(" << iter << ")= " << d << "\t E(" << iter << ")= " << E << " rate " << rate()
          << endl;
  mout.EndBlock();
}

bool NewtonMethod(IAssemble &assemble, Vector &u, bool mute) {
  Newton newton;
  if (!mute) newton.PrintInfo();
  if (!mute) assemble.PrintInfo(u);
  newton.operator()(assemble, u);
  return newton.converged();
}

bool NewtonMethod(IAssemble *assemble, Vector &u, bool mute) {
  return NewtonMethod(*assemble, u, mute);
}

std::shared_ptr<Newton> CreateNewton(const std::string &solverName,
                                     const std::string &preconditionerName) {
  auto linearSolverPtr = GetLinearSolver(solverName, GetPC(preconditionerName));
  auto linearSolver = std::unique_ptr<LinearSolver>(linearSolverPtr);
  return std::make_shared<Newton>(std::move(linearSolver));
}
