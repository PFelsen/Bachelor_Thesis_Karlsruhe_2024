#include "BasicNonLinearSolver.hpp"

RMatrix BasicFunction::computeJacobian(const RVector &x) const {
  int n = x.size();
  RVector d(n);
  RVector xh(x);
  RVector xhh(x);
  RMatrix J(n);
  for (int j = 0; j < n; ++j) {
    xh[j] += h;
    xhh[j] -= h;
    d = (*this)(xh) - (*this)(xhh);
    d *= 0.5 / h;
    for (int i = 0; i < n; ++i)
      J[i][j] = d[i];
    xh[j] -= h;
    xhh[j] += h;
  }
  return J;
}

void BasicNewton::operator()(RVector &x, BasicFunction &f) const {
  double d, d1;
  int n = x.size();
  double eps = Eps + Red * x.norm();
  RVector r(n);
  r = f(x);
  d = r.norm();
  double d_previous = d;
  int LS_cnt = 0;
  RMatrix J(n);
  for (int iter = 0; iter < max_iter; ++iter) {
    vout(3) << "  SmallNewton: r(" << iter << ")= " << r << endl;
    if (d < eps) {
      vout(3) << "  SmallNewtonMethod:  d(" << iter << ") = " << d << endl;
      return;
    }
    vout(4) << "  SmallNewtonMethod:  d(" << iter << ") = " << d << endl;

    // compute Jacobian J and determine correction c
    RMatrix J = f.Jacobian(x);
    vout(9) << "  SmallNewtonMethod:  J(" << iter << ") =\n" << J << endl;
    RBasicLinearSolver linearSolver(J, linearSolverName);
    RVector c = linearSolver.Solve(r);
    vout(6) << "  SmallNewtonMethod:  c(" << iter << ") = " << c << endl;
    vout(7) << "  SmallNewtonMethod: _x(" << iter << ") = " << x << endl;
    x -= c;
    vout(5) << "  SmallNewtonMethod: x_(" << iter << ") = " << x << endl;

    // residual
    r = f(x);
    d = r.norm();
    if (d > d_previous) {
      for (int l = 1; l <= LS_iter; ++l) {
        vout(2) << "  SmallNewtonMethod line search " << l << ": d(" << iter << ") = " << d << endl;
        c *= 0.5;
        x += c;
        r = f(x);
        d = r.norm();
      }
    }
    if (d > d_previous) {
      vout(5) << "  SmallNewtonMethod: line search unsuccessful\n" << endl;
      ++LS_cnt;
      if (LS_cnt == 3) {
        vout(1) << "  SmallNewtonMethod: "
                << "Too many line searches unsuccessful." << endl;
        break;
      }
    }
    d_previous = d;
  }
  if (d > eps) {
    vout(0) << "  No convergence in SmallNewtonMethod;"
            << " Defect = " << d << "  epsilon = " << eps << endl;
    vout(1) << " f(x) = " << f(x);
  }
}