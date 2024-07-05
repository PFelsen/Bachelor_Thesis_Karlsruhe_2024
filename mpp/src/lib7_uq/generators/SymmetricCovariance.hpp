#ifndef M_COVARIANCEFUNCTION_H
#define M_COVARIANCEFUNCTION_H

#include "Assertion.hpp"
#include "RVector.hpp"
#include "CVector.hpp"
#include "RMatrix.hpp"
#include "CMatrix.hpp"
#include "Config.hpp"
#include "RTensor.hpp"


/*
 * Todo padding should be done here
 */

RVector MOCE(RVector &toepRow, RVector &toepCol);

RMatrix BMOCE(RMatrix &toepRows, RMatrix &toepCols);

RTensor BBMOCE(RTensor &toepRows, RTensor &toepCols);

template<typename T>
class SymmetricCovariance {
public:
  virtual double Function(double *tau) = 0;

  virtual void ToeplitzMatrix(T &toepRow, T &toepCol) = 0;

  virtual T EmbeddedToeplitzMatrix(T &toepRow, T &toepCol) = 0;

  virtual ~SymmetricCovariance() {};
};

class CovarianceFunction1D : public SymmetricCovariance<RVector> {
private:
  using T = RVector;

  double sigma = 1.0;

  std::vector<double> lambda = {0.3, 0.3};

  double smoothing = 1.0;

  double norm(const double *x) { return std::abs(x[0]); }

public:
  CovarianceFunction1D() {
    Config::Get("sigma", sigma);
    Config::Get("lambda", lambda);
    Config::Get("smoothing", smoothing);
  }

  CovarianceFunction1D(double sigma, double lambda, double smoothing)
    : sigma(sigma), lambda(1, lambda), smoothing(smoothing) {}

  double Function(double *tau) override {
    tau[0] /= lambda[0];
    return std::pow(sigma, 2) * std::exp(-std::pow(norm(tau), smoothing));
  }

  void ToeplitzMatrix(T &toepRow, T &toepCol) override;

  T EmbeddedToeplitzMatrix(T &toepRow, T &toepCol) override {
    ToeplitzMatrix(toepRow, toepCol);
    return MOCE(toepRow, toepCol);
  }
};

class CovarianceFunction2D : public SymmetricCovariance<RMatrix> {
private:
  using T = RMatrix;

  int norm_p = 2;

  double sigma = 1.0;

  std::vector<double> lambda = {0.15, 0.15};

  double smoothing = 1.5;

  double norm(const double *x) {
    if (norm_p == 1) return std::abs(x[0]) + std::abs(x[1]);
    if (norm_p == 2) return sqrt(std::pow(x[0], 2) + std::pow(x[1], 2));
    else Warning("p not valid, 2-norm is taken")
    return 2;
  }

public:
  CovarianceFunction2D() {
    Config::Get("sigma", sigma);
    Config::Get("norm_p", norm_p);
    Config::Get("lambda", lambda);
    Config::Get("smoothing", smoothing);
  }

  CovarianceFunction2D(double sigma, std::vector<double> lambda,
                       double smoothing, int norm_p)
    : sigma(sigma), lambda(lambda), smoothing(smoothing), norm_p(norm_p) {}

  double Function(double *tau) override {
    tau[0] /= lambda[0];
    tau[1] /= lambda[1];
//      np.power(linalg.norm(tau, ord=self.p) / self.lam, self.alpha))
    return std::pow(sigma, 2) * std::exp(-std::pow(norm(tau), smoothing));
  }

  void ToeplitzMatrix(T &toepRows, T &toepCols) override;

  T EmbeddedToeplitzMatrix(T &toepRows, T &toepCols) override {
    ToeplitzMatrix(toepRows, toepCols);
    return BMOCE(toepRows, toepCols);
  }

  std::string ToString() const {
    return "norm_p=" + std::to_string(norm_p) + "  sigma=" + std::to_string(sigma)
      + "  lambda=" + std::to_string(lambda[0]) + "  lambda=" + std::to_string(lambda[1])
      + "  smoothing=" + std::to_string(smoothing);
  }
};


class CovarianceFunction3D : public SymmetricCovariance<RTensor> {
private:
  using T = RTensor;

  int norm_p = 2;

  double sigma = 1.0;

  std::vector<double> lambda = {0.15, 0.15, 0.15};

  double smoothing = 1.5;

  double norm(const double *x) {
    if (norm_p == 1) return std::abs(x[0]) + std::abs(x[1]) + std::abs(x[2]);
    if (norm_p == 2) return sqrt(std::pow(x[0], 2) + std::pow(x[1], 2) + std::pow(x[2], 2));
    else Warning("p not valid, 2-norm is taken")
    return 2;
  }

public:
  CovarianceFunction3D() {
    Config::Get("sigma", sigma);
    Config::Get("norm_p", norm_p);
    Config::Get("lambda", lambda);
    Config::Get("smoothing", smoothing);
  }

  CovarianceFunction3D(double sigma, std::vector<double> lambda,
                       double smoothing, int norm_p)
    : sigma(sigma), lambda(lambda), smoothing(smoothing), norm_p(norm_p) {}

  double Function(double *tau) override {
    tau[0] /= lambda[0];
    tau[1] /= lambda[1];
    tau[2] /= lambda[2];
    return std::pow(sigma, 2) * std::exp(-std::pow(norm(tau), smoothing));
  }

  void ToeplitzMatrix(T &toepRows, T &toepCols) override;
//need to look if input needs to be changed
  T EmbeddedToeplitzMatrix(T &toepRows, T &toepCols) override {
    ToeplitzMatrix(toepRows, toepCols);
    return BBMOCE(toepRows, toepCols);
  }

  std::string ToString() const {
    return "norm_p=" + std::to_string(norm_p) + "  sigma=" + std::to_string(sigma)
      + "  lambda=" + std::to_string(lambda[0]) + "  lambda=" + std::to_string(lambda[1])
      + "  lambda=" + std::to_string(lambda[2]) + "  smoothing=" + std::to_string(smoothing);
  }
};


std::vector<double> linspace(const double &start, const double &end, int num);

//void CirculantEmbedding1D::padding(RVector &ar,
//                                   RVector &ac) {
//    double delta = (domain_end - domain_start) / (numberOfCells - 1);
//
//    // Add 50% more points
//    int N_padded = numberOfCells + int(numberOfCells / 2);
//
//    double extended_domain_start = domain_start;
//    double extended_domain_end = delta * (N_padded - 1);
//
//    ar.resize(N_padded);
//    ac.resize(N_padded);
//
//    const vector<double> &coordinate1 = linspace(extended_domain_start,
//                                                 extended_domain_end, N_padded);
//    double y_rows = coordinate1.front();
//    for (int i = numberOfCells; i < N_padded; i++) {
//        double x_rows = coordinate1[i];
//        double tau_rows = x_rows - y_rows;
//
//        double x_columns = coordinate1.front();
//        double y_columns = coordinate1[i];
//        double tau_columns = x_columns - y_columns;
//
//        ar[i] = Function(&tau_rows);
//        ac[i] = Function(&tau_columns);
//    }
//}

//void CirculantEmbedding2D::padding(BlockToeplitzRow &ar,
//                                   BlockToeplitzColumn &ac) {
//    double delta1 = (domain1_end - domain1_start) / (numberOfCells[0] - 1);
//    double delta2 = (domain2_end - domain2_start) / (numberOfCells[1] - 1);
//
//    // Add 50% more points
//    int N_padded1 = numberOfCells[0] + int(numberOfCells[0] / 2);
//    int N_padded2 = numberOfCells[1] + int(numberOfCells[1] / 2);
//
//    double extended_domain1_start = domain1_start;
//    double extended_domain1_end = delta1 * (N_padded1 - 1);
//    double extended_domain2_start = domain2_start;
//    double extended_domain2_end = delta2 * (N_padded2 - 1);
//
//    ar.resize(N_padded2);
//    ac.resize(N_padded2);
//
//    const vector<double> &coordinate1 = linspace(extended_domain1_start,
//                                                 extended_domain1_end, N_padded1);
//    const vector<double> &coordinate2 = linspace(extended_domain2_start,
//                                                 extended_domain2_end, N_padded2);
//    double y_rows[2] = {coordinate1.front(), coordinate2.front()};
//    for (int index_c2 = 0; index_c2 < N_padded2; index_c2++) {
//        ar[index_c2].resize(N_padded1);
//        ac[index_c2].resize(N_padded1);
//        for (int index_c1 = numberOfCells[0]; index_c1 < N_padded1; index_c1++) {
//            double c1 = coordinate1[index_c1];
//            double c2 = coordinate2[index_c2];
//            double x_rows[2] = {c1, c2};
//            double tau_rows[2] = {x_rows[0] - y_rows[0],
//                                  x_rows[1] - y_rows[1]};
//
//            double x_columns[2] = {coordinate1.front(), c2};
//            double y_columns[2] = {c1, coordinate2.front()};
//            double tau_columns[2] = {x_columns[0] - y_columns[0],
//                                     x_columns[1] - y_columns[1]};
//
//            ac[index_c2][index_c1] = Function(tau_columns);
//            ar[index_c2][index_c1] = Function(tau_rows);
//        }
//    }
//    for (int index_c2 = numberOfCells[1]; index_c2 < N_padded2; index_c2++) {
//        for (int index_c1 = 0; index_c1 < numberOfCells[0]; index_c1++) {
//            double c1 = coordinate1[index_c1];
//            double c2 = coordinate2[index_c2];
//            double x_rows[2] = {c1, c2};
//            double tau_rows[2] = {x_rows[0] - y_rows[0],
//                                  x_rows[1] - y_rows[1]};
//
//            double x_columns[2] = {coordinate1.front(), c2};
//            double y_columns[2] = {c1, coordinate2.front()};
//            double tau_columns[2] = {x_columns[0] - y_columns[0],
//                                     x_columns[1] - y_columns[1]};
//
//            ac[index_c2][index_c1] = Function(tau_columns);
//            ar[index_c2][index_c1] = Function(tau_rows);
//        }
//    }
//}


#endif //M_COVARIANCEFUNCTION_H
