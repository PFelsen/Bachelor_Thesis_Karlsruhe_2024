#include <fstream>
#include <iomanip>
#include <vector>
#include <l_grad_ari.hpp>
#include <l_nlinsys.hpp>

// To make file include with non-verified quadrature data working
#undef BUILD_IA
#define Exit(s) exit(1)
#include "../../../../src/lib3_disc/quadrature/QuadratureDataTriangle.hpp"

/**
 * For analytical background see:
 * [1] D.A. Dunavant
 *     High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle
 *     International Journal for Numerical Methods in Engineering (1985)
 * [2] F.D. Witherden, P.E. Vincent
 *     On the identification of symmetric quadrature rules for finite element methods
 *
 * The coefficients for non-verified quadrature rules are taken from [1] and [2] respectively.
 * They are listed in lib3_disc/quadrature/QuadratureTriangleData.hpp!
 */

static std::string indent = "    ";

void printDataType(std::ostream &stream, const cxsc::l_interval &z) {
  cxsc::interval w(z);
  stream << SaveOpt << "IAInterval(" << RndDown << Inf(w) << ", " << RndUp << Sup(w) << ")"
         << RestoreOpt;
}

l_interval BLOW(double a, double eps) { return a + l_interval(-1, 1) * eps; }

std::ostream &operator<<(std::ostream &stream, const QtriSymT_data0<l_interval> &d0) {
  stream << "{";
  printDataType(stream, d0.weight);
  return stream << "}";
}

std::ostream &operator<<(std::ostream &stream, const QtriSymT_data1<l_interval> &d1) {
  stream << "{";
  printDataType(stream, d1.weight);
  stream << ", ";
  printDataType(stream, d1.A);
  stream << ", ";
  printDataType(stream, d1.B);
  return stream << "}";
}

std::ostream &operator<<(std::ostream &stream, const QtriSymT_data2<l_interval> &d2) {
  stream << "{";
  printDataType(stream, d2.weight);
  stream << ", ";
  printDataType(stream, d2.A);
  stream << ", ";
  printDataType(stream, d2.B);
  stream << ", ";
  printDataType(stream, d2.C);
  return stream << "}";
}

int orderGlobal, n0Global, n1Global, n2Global, numUnknownsGlobal, numEqnsGlobal;

void PrintData(std::ostream &stream, int order, const QtriSymT_data<l_interval> &qData) {
  stream << SetPrecision(18, 15);
  stream << indent << "case " << order << ":" << endl;
  if (!(qData.d0.empty() && qData.d1.empty() && qData.d2.empty())) {
    stream << indent << "  return QtriSymT_data<IAInterval>(" << endl << indent << "      {";
    if (qData.d0.size() > 0) { stream << qData.d0[0]; }
    stream << "}," << endl << indent << "      {";
    for (int i = 0; i < qData.d1.size(); ++i) {
      if (i != 0) { stream << "," << endl << indent << "       "; }
      stream << qData.d1[i];
    }
    stream << "}," << endl << indent << "      {";
    for (int i = 0; i < qData.d2.size(); ++i) {
      if (i != 0) { stream << "," << endl << indent << "       "; }
      stream << qData.d2[i];
    }
    stream << "});" << endl;
  }
}

l_interval w0_fixed;

template<typename VEC>
void eqn(const VEC &x, VEC &r, int n0, int n1, int n2, int m, int j, int /*3*/ k, int v_c, int v_d,
         int idCorrection) {
  int shift = 1 + n0 + idCorrection;
  if (m == 1) {
    if (n0 > 0) {
      if (idCorrection == -1) {
        r[1] = r[1] + w0_fixed;
      } else {
        r[1] = r[1] + x[1];
      }
    }
    for (int i = 0; i < n1; ++i)
      r[1] = r[1] + x[shift + 2 * i];
    for (int i = 0; i < n2; ++i)
      r[1] = r[1] + x[shift + 2 * n1 + 3 * i];
    r[1] = r[1] - 1;
    return;
  }
  if (/*3*/ k == 0) {
    if (j == 2) {
      for (int i = 0; i < n1; ++i)
        r[m] = r[m] + x[shift + 2 * i] * sqr(x[shift + 1 + 2 * i]);
      for (int i = 0; i < n2; ++i)
        r[m] = r[m] + x[shift + 2 * n1 + 3 * i] * sqr(x[shift + 1 + 2 * n1 + 3 * i]);
      r[m] = r[m] - l_interval(v_c) / l_interval(v_d);
    } else {
      for (int i = 0; i < n1; ++i)
        r[m] = r[m] + x[shift + 2 * i] * power(x[shift + 1 + 2 * i], j);
      for (int i = 0; i < n2; ++i)
        r[m] = r[m] + x[shift + 2 * n1 + 3 * i] * power(x[shift + 1 + 2 * n1 + 3 * i], j);
      r[m] = r[m] - l_interval(v_c) / l_interval(v_d);
    }
  } else {
    for (int i = 0; i < n1; ++i)
      r[m] = r[m] + x[shift + 2 * i] * power(x[shift + 1 + 2 * i], j);
    for (int i = 0; i < n2; ++i)
      r[m] = r[m]
             + x[shift + 2 * n1 + 3 * i] * power(x[shift + 1 + 2 * n1 + 3 * i], j)
                   * cos(k * x[shift + 2 + 2 * n1 + 3 * i]);
    r[m] = r[m] - l_interval(v_c) / l_interval(v_d);
  }
}

// Especially for order 16
template<typename VEC>
void eqn16(const VEC &x, VEC &r, int n0, int n1, int n2, int m, int j, int /*3*/ k, int v_c,
           int v_d) {
  int shift = 2;
  if (m == 1) {
    r[1] = r[1] + x[1];
    for (int i = 0; i < n1; ++i)
      r[1] = r[1] + x[shift + 2 * i];
    shift--;
    for (int i = 0; i < n2; ++i)
      r[1] = r[1] + x[shift + 2 * n1 + 3 * i];
    r[1] = r[1] - 1;
    return;
  }
  if (/*3*/ k == 0) {
    if (j == 2) {
      for (int i = 0; i < n1 - 1; ++i)
        r[m] = r[m] + x[shift + 2 * i] * sqr(x[shift + 1 + 2 * i]);
      r[m] = r[m] + x[shift + 2 * (n1 - 1)] * sqr(l_interval(0.5));
      shift--;
      for (int i = 0; i < n2; ++i)
        r[m] = r[m] + x[shift + 2 * n1 + 3 * i] * sqr(x[shift + 1 + 2 * n1 + 3 * i]);
      r[m] = r[m] - l_interval(v_c) / l_interval(v_d);
    } else {
      for (int i = 0; i < n1 - 1; ++i)
        r[m] = r[m] + x[shift + 2 * i] * power(x[shift + 1 + 2 * i], j);
      r[m] = r[m] + x[shift + 2 * (n1 - 1)] * power(l_interval(0.5), j);
      shift--;
      for (int i = 0; i < n2; ++i)
        r[m] = r[m] + x[shift + 2 * n1 + 3 * i] * power(x[shift + 1 + 2 * n1 + 3 * i], j);
      r[m] = r[m] - l_interval(v_c) / l_interval(v_d);
    }
  } else {
    for (int i = 0; i < n1 - 1; ++i)
      r[m] = r[m] + x[shift + 2 * i] * power(x[shift + 1 + 2 * i], j);
    r[m] = r[m] + x[shift + 2 * (n1 - 1)] * power(l_interval(0.5), j);
    shift--;
    for (int i = 0; i < n2; ++i)
      r[m] = r[m]
             + x[shift + 2 * n1 + 3 * i] * power(x[shift + 1 + 2 * n1 + 3 * i], j)
                   * cos(k * x[shift + 2 + 2 * n1 + 3 * i]);
    r[m] = r[m] - l_interval(v_c) / l_interval(v_d);
  }
}

/// see Table I. Polar moments in [1]
static std::vector<std::vector<int>> table = {
    {/*j, 3k, counter of v, denominator of v */}, // m = 0 not used
    {0, 0, 1, 1},                                 // m = 1
    {2, 0, 1, 4},                                 // m = 2
    {3, 3, -1, 10},                               // m = 3
    {4, 0, 1, 10},                                // m = 4
    {5, 3, -2, 35},                               // m = 5
    {6, 0, 29, 560},                              // m = 6
    {6, 6, 1, 28},                                // m = 7
    {7, 3, -1, 28},                               // m = 8
    {8, 0, 11, 350},                              // m = 9
    {8, 6, 1, 40},                                // m = 10
    {9, 3, -37, 1540},                            // m = 11
    {9, 9, -1, 55},                               // m = 12
    {10, 0, 13, 616},                             // m = 13
    {10, 6, 1, 55},                               // m = 14
    {11, 3, -49, 2860},                           // m = 15
    {11, 9, -2, 143},                             // m = 16
    {12, 0, 425, 28028},                          // m = 17
    {12, 6, 137, 10010},                          // m = 18
    {12, 12, 1, 91},                              // m = 19
    {13, 3, -64, 5005},                           // m = 20
    {13, 9, -1, 91},                              // m = 21
    {14, 0, 523, 45760},                          // m = 22
    {14, 6, 85, 8008},                            // m = 23
    {14, 12, 1, 112},                             // m = 24
    {15, 3, -6733, 680680},                       // m = 25
    {15, 9, -109, 12376},                         // m = 26
    {15, 15, -1, 136},                            // m = 27
    {16, 0, 217, 24310},                          // m = 28
    {16, 6, 209, 24752},                          // m = 29
    {16, 12, 1, 136},                             // m = 30
    {17, 3, -2909, 369512},                       // m = 31
    {17, 9, -65, 9044},                           // m = 32
    {17, 15, -2, 323},                            // m = 33
    {18, 0, 66197, 9237800},                      // m = 34
    {18, 6, 8069, 1175720},                       // m = 35
    {18, 12, 317, 51680},                         // m = 36
    {18, 18, 1, 190},                             // m = 37
    {19, 3, -3769, 587860},                       // m = 38
    {19, 9, -77, 12920},                          // m = 39
    {19, 15, -1, 190},                            // m = 40
    {20, 0, 83651, 14226212},                     // m = 41
    {20, 6, 11303, 1989680},                      // m = 42
    {20, 12, 92, 17765},                          // m = 43
    {20, 18, 1, 220}                              // m = 44
};

/// see Table I. Polar moments in [1]
static std::vector<int> tableNumEqns = {
    // 1  2  3  4  5  6  7  8   9   10  11  12  13  14  15  16  14  18  19  10
    -1, 1, 2, 3, 4, 5, 7, 8, 10, 12, 14, 16, 19, 21, 24, 27, 30, 33, 37, 40, 44};

l_GTvector f(const l_GTvector &x) {
  l_GTvector r(numEqnsGlobal);
  for (int m = 1; m <= numEqnsGlobal; ++m) {
    r[m] = 0.0;
    if (orderGlobal == 16) {
      eqn16(x, r, n0Global, n1Global, n2Global, m, table[m][0], table[m][1], table[m][2],
            table[m][3]);
    } else {
      eqn(x, r, n0Global, n1Global, n2Global, m, table[m][0], table[m][1], table[m][2], table[m][3],
          numEqnsGlobal - numUnknownsGlobal);
    }
  }
  return r;
}

l_ivector f(const l_ivector &x) {
  l_ivector r(numEqnsGlobal);
  for (int m = 1; m <= numEqnsGlobal; ++m) {
    r[m] = 0.0;
    if (orderGlobal == 16) {
      eqn16(x, r, n0Global, n1Global, n2Global, m, table[m][0], table[m][1], table[m][2],
            table[m][3]);
    } else {
      eqn(x, r, n0Global, n1Global, n2Global, m, table[m][0], table[m][1], table[m][2], table[m][3],
          numEqnsGlobal - numUnknownsGlobal);
    }
  }
  return r;
}

l_ivector SetUpVector(int order, const QtriSymT_data<double> &qData, bool debugOut,
                      double eps = 1e-8) {
  cout << "Setting up data for verification of quadrature rule..." << endl;
  orderGlobal = order;
  n0Global = qData.d0.size();
  n1Global = qData.d1.size();
  n2Global = qData.d2.size();
  numUnknownsGlobal = n0Global + 2 * n1Global + 3 * n2Global;
  numEqnsGlobal = tableNumEqns[order];
  cout << "  >> (n0, n1, n2)=        (" << n0Global << ", " << n1Global << ", " << n2Global << ")"
       << endl
       << "  >> Number of equations= " << numEqnsGlobal << endl
       << "  >> Number of unknowns=  " << numUnknownsGlobal << endl;

  // Set size of unknown-vector to numEquations since algorithm requires square system. Thus, in
  // case of larger amount of unknowns the number has to be decreased accordingly.
  // Therefore, the weight for point (0.33, 0.33), i.e., n0, is assumed to be good enough for the
  // quadrature rule and can be used as it is. The remaining unknowns can then be computed as
  // described in [1]
  l_ivector x(numEqnsGlobal);

  int idCorrection = numEqnsGlobal - numUnknownsGlobal;
  if (!(idCorrection == 0 || idCorrection == -1)) {
    cout << "  >> Error: Number of unknowns does not fit to number of equations: (n_u= "
         << numUnknownsGlobal << "; n_e= " << numEqnsGlobal << ")!" << endl
         << "            Computation will be aborted!" << endl;
    return x;
  }
  if (idCorrection == -1 && qData.d0.size() == 0) {
    cout << "  >> Error: Tried to adapt number of unknowns but size of d0 is zero!" << endl
         << "            Computation will be aborted!" << endl;
    return x;
  }

  if (idCorrection == -1) {
    cout << "  >> Number of unknowns needs to be adapted. Weight w0 is considered to be enclosed "
            "correctly"
         << endl;
  }

  if (qData.d0.size() > 0) {
    if (idCorrection == -1) {
      w0_fixed = 2 * BLOW(qData.d0[0].weight, 1e-21);
      if (debugOut) {
        cout << "n_0: w= " << w0_fixed << " (not considered in computations)" << endl;
      }
    } else {
      x[1] = 2 * BLOW(qData.d0[0].weight, eps);
      if (debugOut) { cout << "n_0: w= " << x[1] << endl; }
    }
  }
  for (int i = 0; i < qData.d1.size(); ++i) {
    int shift = qData.d0.size() + 2 * i + idCorrection;
    x[1 + shift] = 6 * BLOW(qData.d1[i].weight, eps);
    x[2 + shift] = (1.0 - 3 * BLOW(qData.d1[i].A, eps)) / 2.0;
    if (debugOut) {
      cout << "n_1: w= " << x[1 + shift] << endl << "     r= " << x[2 + shift] << endl;
    }
  }
  for (int i = 0; i < qData.d2.size(); ++i) {
    int shift = qData.d0.size() + 2 * qData.d1.size() + 3 * i + idCorrection;
    l_interval A = BLOW(qData.d2[i].A, eps);
    l_interval B = BLOW(qData.d2[i].B, eps);
    l_interval alpha_i = atan(sqrt(l_interval(3.0)) * (1.0 - A - 2 * B) / (1.0 - 3 * A));
    l_interval r_i = (1.0 - 3 * A) / 2.0 / cos(alpha_i);
    x[1 + shift] = 12 * BLOW(qData.d2[i].weight, eps);
    x[2 + shift] = r_i;
    x[3 + shift] = alpha_i;
    if (debugOut) {
      cout << "n_2: w= " << x[1 + shift] << endl
           << "     r= " << x[2 + shift] << endl
           << "     a= " << x[3 + shift] << endl;
    }
  }
  return x;
}

// Especially for order 16
l_ivector SetUpVector16(const QtriSymT_data<double> &qData, bool debugOut, double eps = 1e-8) {
  cout << "Setting up data for verification of quadrature rule..." << endl;
  orderGlobal = 16;
  n0Global = qData.d0.size();
  n1Global = qData.d1.size();
  n2Global = qData.d2.size();
  numUnknownsGlobal = n0Global + 2 * n1Global + 3 * n2Global;
  numEqnsGlobal = tableNumEqns[16];
  cout << "  >> (n0, n1, n2)=        (" << n0Global << ", " << n1Global << ", " << n2Global << ")"
       << endl
       << "  >> Number of equations= " << numEqnsGlobal << endl
       << "  >> Number of unknowns=  " << numUnknownsGlobal << endl;

  // Set size of unknown-vector to numEquations since algorithm requires square system. Thus, in
  // case of larger amount of unknowns the number has to be decreased accordingly.
  // This is also required for order 16, however, in this special case, we have to use a different
  // approach since for n1 the value A is exactly zero for one pair. This can be exploited here
  // to reduce the number of unknowns. The remaining unknowns can then be computed as
  // described in [1]
  l_ivector x(numEqnsGlobal);

  cout << "  >> Number of unknowns needs to be adapted. For n1 one r is exactly zero." << endl;

  x[1] = 2 * BLOW(qData.d0[0].weight, eps);
  if (debugOut) { cout << "n_0: w= " << x[1] << endl; }
  for (int i = 0; i < qData.d1.size() - 1; ++i) {
    int shift = 1 + 2 * i;
    x[1 + shift] = 6 * BLOW(qData.d1[i].weight, eps);
    x[2 + shift] = (1.0 - 3 * BLOW(qData.d1[i].A, eps)) / 2.0;
    if (debugOut) {
      cout << "n_1: w= " << x[1 + shift] << endl << "     r= " << x[2 + shift] << endl;
    }
  }
  x[2 * qData.d1.size()] = 6 * BLOW(qData.d1[qData.d1.size() - 1].weight, eps);
  if (debugOut) {
    cout << "n_1: w= " << x[2 * qData.d1.size()] << endl
         << "     r= " << l_interval(0.5) << " (not considered in computations)" << endl;
  }
  for (int i = 0; i < qData.d2.size(); ++i) {
    int shift = 2 * qData.d1.size() + 3 * i;
    l_interval A = BLOW(qData.d2[i].A, eps);
    l_interval B = BLOW(qData.d2[i].B, eps);
    l_interval alpha_i = atan(sqrt(l_interval(3.0)) * (1.0 - A - 2 * B) / (1.0 - 3 * A));
    l_interval r_i = (1.0 - 3 * A) / 2.0 / cos(alpha_i);
    x[1 + shift] = 12 * BLOW(qData.d2[i].weight, eps);
    x[2 + shift] = r_i;
    x[3 + shift] = alpha_i;
    if (debugOut) {
      cout << "n_2: w= " << x[1 + shift] << endl
           << "     r= " << x[2 + shift] << endl
           << "     a= " << x[3 + shift] << endl;
    }
  }
  return x;
}

l_ivector ImproveInitialVector(const l_ivector &x, double eps = 1e-9) {
  cout << "Perform Newton steps to improve solution:" << endl;
  l_rvector Y(numEqnsGlobal);
  for (int i = 1; i <= numEqnsGlobal; ++i)
    Y[i] = mid(x[i]);
  for (int i = 0; i < 10; ++i) {
    cout << "  >> Step " << (i + 1) << " ..." << endl;
    NewtonStep(f, Y);
  }
  cout << "  >> Newton steps finished!" << endl;
  l_ivector X(numEqnsGlobal);
  for (int i = 1; i <= numEqnsGlobal; ++i)
    X[i] = Y[i] + l_interval(-1, 1) * eps;
  return X;
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

real CalcRelDiamMax(const QtriSymT_data<l_interval> &qData) {
  real d = 0.0;
  if (qData.d0.size() > 0) { d = max(d, RelDiam(interval(qData.d0[0].weight))); }
  for (int i = 0; i < qData.d1.size(); ++i) {
    d = max(d, RelDiam(interval(qData.d1[i].weight)));
    d = max(d, RelDiam(interval(qData.d1[i].A)));
    d = max(d, RelDiam(interval(qData.d1[i].B)));
  }
  for (int i = 0; i < qData.d2.size(); ++i) {
    d = max(d, RelDiam(interval(qData.d2[i].weight)));
    d = max(d, RelDiam(interval(qData.d2[i].A)));
    d = max(d, RelDiam(interval(qData.d2[i].B)));
    d = max(d, RelDiam(interval(qData.d2[i].C)));
  }
  return d;
}

QtriSymT_data<l_interval> RecalcData(const SystemSolPair &sol) {
  if (!sol.second) {
    cout << "  >> Error: No verified quadrature data computed!" << endl;

    l_ivector defect = f(sol.first);
    cout << "  >> Defect vector:" << endl;
    for (int i = 0; i < numEqnsGlobal; ++i) {
      cout << "     " << i << (i > 9 ? ":  " : ":   ") << defect[i + 1] << endl;
    }
    cout << endl;
    return QtriSymT_data<l_interval>{{}, {}, {}};
  }
  cout << "Recalculating verified quadrature data..." << endl;
  QtriSymT_data<l_interval> qData{std::vector<QtriSymT_data0<l_interval>>(n0Global),
                                  std::vector<QtriSymT_data1<l_interval>>(n1Global),
                                  std::vector<QtriSymT_data2<l_interval>>(n2Global)};

  int idCorrection = numEqnsGlobal - numUnknownsGlobal;
  if (n0Global > 0) {
    if (idCorrection == -1) {
      qData.d0[0].weight = w0_fixed / 2.0;
    } else {
      qData.d0[0].weight = sol.first[1] / 2.0;
    }
  }
  int shift = 1 + n0Global + idCorrection;
  for (int i = 0; i < n1Global; ++i) {
    qData.d1[i].weight = sol.first[shift + 2 * i] / 6.0;
    qData.d1[i].A = (1.0 - 2.0 * sol.first[shift + 1 + 2 * i]) / l_interval(3);
    qData.d1[i].B = (1.0 - qData.d1[i].A) / l_interval(2);
  }
  for (int i = 0; i < n2Global; ++i) {
    l_interval r_i = sol.first[shift + 1 + 2 * n1Global + 3 * i];
    l_interval alpha_i = sol.first[shift + 2 + 2 * n1Global + 3 * i];
    qData.d2[i].weight = sol.first[shift + 2 * n1Global + 3 * i] / 12.0;
    qData.d2[i].A = (1 - 2 * r_i * cos(alpha_i)) / l_interval(3);
    qData.d2[i].B =
        (1 + r_i * cos(alpha_i) - sqrt(l_interval(3.0)) * r_i * sin(alpha_i)) / l_interval(3);
    qData.d2[i].C = 1 - qData.d2[i].A - qData.d2[i].B;
  }
  cout << "  >> Maximal diameter = " << CalcRelDiamMax(qData) << endl;
  return qData;
}

// Especially for order 16
QtriSymT_data<l_interval> RecalcData16(const SystemSolPair &sol) {
  if (!sol.second) {
    cout << "  >> Error: No verified quadrature data computed!" << endl;

    l_ivector defect = f(sol.first);
    cout << "  >> Defect vector:" << endl;
    for (int i = 0; i < numEqnsGlobal; ++i) {
      cout << "     " << i << (i > 9 ? ":  " : ":   ") << defect[i + 1] << endl;
    }
    cout << endl;
    return QtriSymT_data<l_interval>{{}, {}, {}};
  }
  cout << "Recalculating verified quadrature data..." << endl;
  QtriSymT_data<l_interval> qData{std::vector<QtriSymT_data0<l_interval>>(n0Global),
                                  std::vector<QtriSymT_data1<l_interval>>(n1Global),
                                  std::vector<QtriSymT_data2<l_interval>>(n2Global)};

  qData.d0[0].weight = sol.first[1] / 2.0;

  int shift = 2;
  for (int i = 0; i < n1Global - 1; ++i) {
    qData.d1[i].weight = sol.first[shift + 2 * i] / 6.0;
    qData.d1[i].A = (1.0 - 2.0 * sol.first[shift + 1 + 2 * i]) / l_interval(3);
    qData.d1[i].B = (1.0 - qData.d1[i].A) / l_interval(2);
  }
  qData.d1[n1Global - 1].weight = sol.first[2 * n1Global] / 6.0;
  qData.d1[n1Global - 1].A = l_interval(0.0);
  qData.d1[n1Global - 1].B = l_interval(0.5);

  shift--;
  for (int i = 0; i < n2Global; ++i) {
    l_interval r_i = sol.first[shift + 1 + 2 * n1Global + 3 * i];
    l_interval alpha_i = sol.first[shift + 2 + 2 * n1Global + 3 * i];
    qData.d2[i].weight = sol.first[shift + 2 * n1Global + 3 * i] / 12.0;
    qData.d2[i].A = (1 - 2 * r_i * cos(alpha_i)) / l_interval(3);
    qData.d2[i].B =
        (1 + r_i * cos(alpha_i) - sqrt(l_interval(3.0)) * r_i * sin(alpha_i)) / l_interval(3);
    qData.d2[i].C = 1 - qData.d2[i].A - qData.d2[i].B;
  }
  cout << "  >> Maximal diameter = " << CalcRelDiamMax(qData) << endl;
  return qData;
}

QtriSymT_data<l_interval> runVerification(int order, bool debugOut = false) {
  cout << "Starting verification of quadrature rule with order " << order << "..." << endl
       << "=======================================================" << (order > 9 ? "==" : "=")
       << endl;

  QtriSymT_data<double> qDataDouble = GetQtriSymT_data(order);

  l_ivector x;
  if (order == 16) {
    x = SetUpVector16(qDataDouble, debugOut, 1e-8);
  } else {
    x = SetUpVector(order, qDataDouble, debugOut, 1e-8);
  }
  l_ivector x_improved = ImproveInitialVector(x, order > 18 ? 1e-14 : 1e-9);
  SystemSolPair sol = SolveNonlinearSystem(x_improved, order > 18 ? 1e-11 : 1e-20); // 1e-18
  QtriSymT_data<l_interval> qDataInterval;
  if (order == 16) {
    qDataInterval = RecalcData16(sol);
  } else {
    qDataInterval = RecalcData(sol);
  }

  cout << endl;
  return qDataInterval;
}

int main() {
  stagprec = 10;

  std::vector<QtriSymT_data<l_interval>> qData(21);
  for (int i = 1; i < qData.size(); ++i) {
    qData[i] = runVerification(i);
  }

  cout << "Verified quadrature data:" << endl;
  for (int i = 1; i < qData.size(); ++i) {
    PrintData(cout, i, qData[i]);
  }

  return 0;
}
