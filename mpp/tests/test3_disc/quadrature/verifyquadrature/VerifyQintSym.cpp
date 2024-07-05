#include <fstream>
#include <iomanip>
#include <vector>
#include <l_grad_ari.hpp>
#include <l_nlinsys.hpp>

// To make file include with non-verified quadrature data working
#undef BUILD_IA
#define Exit(s) exit(1)
#include "../../../../src/lib3_disc/quadrature/QuadratureDataInterval.hpp"

static std::string indent = "    ";

void printDataType(std::ostream &stream, const cxsc::l_interval &z) {
  cxsc::interval w(z);
  stream << SaveOpt << "IAInterval(" << RndDown << Inf(w) << ", " << RndUp << Sup(w) << ")"
         << RestoreOpt;
}

l_interval BLOW(double a, double eps) { return a + l_interval(-1, 1) * eps; }

std::ostream &operator<<(std::ostream &stream, const QintSymT_data0<l_interval> &d0) {
  stream << "{";
  printDataType(stream, d0.weight);
  return stream << "}";
}

std::ostream &operator<<(std::ostream &stream, const QintSymT_data1<l_interval> &d1) {
  stream << "{";
  printDataType(stream, d1.weight);
  stream << ", ";
  printDataType(stream, d1.A);
  stream << ", ";
  printDataType(stream, d1.B);
  return stream << "}";
}

int orderGlobal;

void PrintData(std::ostream &stream, int order, const QintSymT_data<l_interval> &qData) {
  stream << SetPrecision(18, 15);
  stream << indent << "case " << (2 * order - 2) << ":" << endl
         << indent << "case " << (2 * order - 1) << ":" << endl;
  if (!(qData.d0.empty() && qData.d1.empty())) {
    stream << indent << "  return QintSymT_data<IAInterval>(" << endl << indent << "      {";
    if (qData.d0.size() > 0) { stream << qData.d0[0]; }
    stream << "}," << endl << indent << "      {";
    for (int i = 0; i < qData.d1.size(); ++i) {
      if (i != 0) { stream << "," << endl << indent << "       "; }
      stream << qData.d1[i];
    }
    stream << "});" << endl;
  }
}

template<typename VEC>
VEC eqn(const VEC &x, int order) {
  VEC r(order);
  for (int i = 1; i <= order; ++i)
    r[i] = 0.0;

  int n = order / 2;

  if (order % 2 == 1) {
    int cnt = 1;
    for (int q = 0; q < order; ++q) {
      if (q == 0) {
        for (int k = 1; k <= n; ++k)
          r[cnt] = r[cnt] + x[k];
        r[cnt] = 2 * r[cnt] + x[order] - 2;
      } else {
        for (int k = 1; k <= n; ++k) {
          if (q == 1) r[cnt] = r[cnt] + x[k] * sqr(x[k + n]);
          else r[cnt] = r[cnt] + x[k] * power(x[k + n], 2 * q);
        }
        r[cnt] = (2 * q + 1) * r[cnt] - 1;
      }
      cnt++;
    }
  } else {
    int cnt = 1;
    for (int q = 0; q < order; ++q) {
      for (int k = 1; k <= n; ++k) {
        if (q == 0) r[cnt] = r[cnt] + x[k];
        else if (q == 1) r[cnt] = r[cnt] + x[k] * sqr(x[k + n]);
        else r[cnt] = r[cnt] + x[k] * power(x[k + n], 2 * q);
      }
      r[cnt] = (2 * q + 1) * r[cnt] - 1;
      cnt++;
    }
  }
  return r;
}

l_GTvector f(const l_GTvector &x) { return eqn(x, orderGlobal); }

l_ivector f(const l_ivector &x) { return eqn(x, orderGlobal); }

l_ivector SetUpVector(int order, const QintSymT_data<double> &qData, bool debugOut,
                      double eps = 1e-8) {
  cout << "Setting up data for verification of quadrature rule..." << endl;
  orderGlobal = order;
  l_ivector x(order);

  if (order % 2 == 1) {
    for (int i = 1; i <= order / 2; ++i) {
      x[i] = 2 * BLOW(qData.d1[i - 1].weight, eps);
      x[order / 2 + i] = 2 * BLOW(qData.d1[i - 1].A, eps) - 1.0;
    }
    x[order] = 2 * BLOW(qData.d0[0].weight, eps);
  } else {
    for (int i = 1; i <= order / 2; ++i) {
      x[i] = 2 * BLOW(qData.d1[i - 1].weight, eps);
      x[order / 2 + i] = 2 * BLOW(qData.d1[i - 1].A, eps) - 1.0;
    }
  }

  return x;
}

using SystemSolPair = std::pair<l_ivector, bool>;

SystemSolPair SolveNonlinearSystem(l_ivector &x, real tol, int numSol = 1) {
  cout << "Solving nonlinear system with " << Ub(x) << " unknowns..." << endl;
  l_imatrix sol;
  intvector unique;
  int N, error;

  try {
    l_AllNLSS(f, x, tol, sol, unique, N, error, numSol);
  } catch (...) {
    cout << "  >> Error: Unknown error!" << endl;
    return SystemSolPair(x, false);
  }

  if (error) {
    if (error != 2) {
      cout << "  >> Error: " << l_AllNLSSErrMsg(error) << endl;
      return SystemSolPair(x, false);
    } else {
      cout << "  >> " << l_AllNLSSErrMsg(error) << endl;
    }
  }

  if (N == 0) {
    cout << "  >> Error: No zero in given interval!" << endl;
    return SystemSolPair(x, false);
  }

  for (int i = 1; i <= N; ++i) {
    if (unique[i]) {
      cout << "  >> Successful: Unique zero!" << endl;
      return SystemSolPair(sol[i], true);
    }
  }
  cout << "  >> Successful: But no unique zero found!" << endl;
  return SystemSolPair(sol[1], true);
}

real CalcRelDiamMax(const QintSymT_data<l_interval> &qData) {
  real d = 0.0;
  if (qData.d0.size() > 0) { d = max(d, RelDiam(interval(qData.d0[0].weight))); }
  for (int i = 0; i < qData.d1.size(); ++i) {
    d = max(d, RelDiam(interval(qData.d1[i].weight)));
    d = max(d, RelDiam(interval(qData.d1[i].A)));
    d = max(d, RelDiam(interval(qData.d1[i].B)));
  }
  return d;
}

QintSymT_data<l_interval> RecalcData(const SystemSolPair &sol) {
  if (!sol.second) {
    cout << "  >> Error: No verified quadrature data computed!" << endl;

    l_ivector defect = f(sol.first);
    cout << "  >> Defect vector:" << endl;
    for (int i = 0; i < orderGlobal; ++i) {
      cout << "     " << i << (i > 9 ? ":  " : ":   ") << defect[i + 1] << endl;
    }
    cout << endl;
    return QintSymT_data<l_interval>{{}, {}};
  }
  cout << "Recalculating verified quadrature data..." << endl;
  QintSymT_data<l_interval> qData{std::vector<QintSymT_data0<l_interval>>(orderGlobal % 2),
                                  std::vector<QintSymT_data1<l_interval>>(orderGlobal / 2)};

  if (orderGlobal % 2 == 1) {
    for (int i = 1; i <= orderGlobal / 2; ++i) {
      qData.d1[i - 1].weight = sol.first[i] / 2.0;
      qData.d1[i - 1].A = (1.0 + sol.first[i + orderGlobal / 2]) / 2.0;
      qData.d1[i - 1].B = (1.0 - sol.first[i + orderGlobal / 2]) / 2.0;
    }
    qData.d0[0].weight = sol.first[orderGlobal] / 2.0;
  } else {
    for (int i = 1; i <= orderGlobal / 2; ++i) {
      qData.d1[i - 1].weight = sol.first[i] / 2.0;
      qData.d1[i - 1].A = (1.0 + sol.first[i + orderGlobal / 2]) / 2.0;
      qData.d1[i - 1].B = (1.0 - sol.first[i + orderGlobal / 2]) / 2.0;
    }
  }
  cout << "  >> Maximal diameter = " << CalcRelDiamMax(qData) << endl;
  return qData;
}

QintSymT_data<l_interval> runVerification(int order, bool debugOut = false) {
  cout << "Starting verification of quadrature rule with order " << order << "..." << endl
       << "=======================================================" << (order > 9 ? "==" : "=")
       << endl;

  QintSymT_data<double> qDataDouble = GetQintSymT_data(2 * order - 1);

  l_ivector x = SetUpVector(order, qDataDouble, debugOut, 1e-8);
  SystemSolPair sol = SolveNonlinearSystem(x, 1e-20); // 1e-18
  QintSymT_data<l_interval> qDataInterval;
  qDataInterval = RecalcData(sol);

  cout << endl;
  return qDataInterval;
}

int main() {
  stagprec = 10;

  std::vector<QintSymT_data<l_interval>> qData(11);
  for (int i = 1; i < qData.size(); ++i) {
    qData[i] = runVerification(i);
  }

  cout << "Verified quadrature data:" << endl;
  for (int i = 1; i < qData.size(); ++i) {
    PrintData(cout, i, qData[i]);
  }

  return 0;
}
