#ifndef VISCOELASTICPROBLEMS_HPP
#define VISCOELASTICPROBLEMS_HPP

#include "ProblemBase.hpp"

class TProblemViscoElastic : public ProblemBase {
public:
  double rho;
  double lambda;
  double mu;
  double tau;
  int dim;
  double tau_p;
  double tau_s;
  int numL;

  TProblemViscoElastic(int d = 3, int nL = 1) : ProblemBase(),
                                                dim(d), rho(1.0), lambda(1.0), mu(1.0), tau(1.0),
                                                numL(nL),
                                                tau_p(0.1), tau_s(0.1) {
    Config::Get("rho", rho);
    Config::Get("lambda", lambda);
    Config::Get("mu", mu);
    Config::Get("tau", tau);
    Config::Get("numL", numL);
  }

  std::string Name() const override {
    return "TProblemViscoElastic";
  }

  virtual bool HasExactSolution() const { return false; }

  virtual double Rho(const Point &x) const { return rho; }

  virtual double Kappa(const Point &x) const { return 2 / 3 * mu + lambda; }

  virtual double Mu(const Point &x, int var = -1) const {
    if (var == -1) return mu * (1 + numL * tau_s);
    if (var == 1) return mu;
    else return mu * tau_s;
  }

  virtual double Lambda(const Point &x, int var = -1) const {
    double kappa = Kappa(x);
    if (var == -1)
      return kappa * (1 + numL * tau_p) - 2 / 3 * mu * (1 + numL * tau_s);
    if (var == 1) return kappa - 2 / 3 * mu;
    else return kappa * tau_p - 2 / 3 * mu * tau_s;

    return lambda;
  }

  virtual double Tau(const Point &x, int var = -1) const { return tau; }

  virtual double C_s(const Point &x) const { return sqrt(Mu(x) / Rho(x)); }

  virtual double C_p(const Point &x) const {
    return sqrt((2.0 * Mu(x) + Lambda(x)) / Rho(x));
  }

  virtual bool HasRHS() const { return false; }

  virtual double F(const Point &x, COMPONENT comp) const { return 0.0; }

  virtual int BndID(const Point &z) const { return 1; }

  virtual double ut(const Point &x, COMPONENT comp) const {
    Exit("Not implemented")
  }

  virtual double ut(const Point &x, const Cell &c, COMPONENT comp) const {
    Exit("Not implemented")
  }

  VectorField ut_Velocity(const Point &x) const {
    return VectorField(ut(x, COMPONENT::V_X), ut(x, COMPONENT::V_Y));
  }

  Tensor ut_Stress(const Point &x) const {
    Exit("fix! dimension");
    return Tensor(ut(x, COMPONENT::S0_x), ut(x, COMPONENT::S0_xy),
                  ut(x, COMPONENT::S0_xy), ut(x, COMPONENT::S0_y));
  }

  virtual double ut_dual(const Point &, COMPONENT comp, const Point &, const Point &) const {
    return 0;
  }

  int Dim() const {
    int td = std::max(1, dim + (dim - 1) + (dim - 2));
    return dim + td + td * numL;
  }

  int nL() const {
    return numL;
  }

  std::vector<COMPONENT> GetComponents() const override {
    if (dim == 2) {
      return std::vector{COMPONENT::V_X, COMPONENT::V_Y,
                         COMPONENT::S0_x, COMPONENT::S0_y, COMPONENT::S0_xy,
                         COMPONENT::S1_x, COMPONENT::S1_y, COMPONENT::S1_xy,
                         COMPONENT::S2_x, COMPONENT::S2_y, COMPONENT::S2_xy,
                         COMPONENT::S3_x, COMPONENT::S3_y, COMPONENT::S3_xy,
                         COMPONENT::S4_x, COMPONENT::S4_y, COMPONENT::S4_xy};
    }
    Exit("Check");
  }

/*
  std::string dimensionName(int m) const override {
    if (m < dim) {
      if (m == 0) { return "V_x"; }
      else if (m == 1) { return "V_y"; }
      else if (m == 2) { return "V_z"; }
    } else if (dim == 2) {
      if (m == dim) { return "S0_x"; }
      else if (m == dim + 1) { return "S0_y"; }
      else if (m == dim + 2) { return "S0_xy"; }

      else if (m == dim + 3) { return "S1_x"; }
      else if (m == dim + 4) { return "S1_y"; }
      else if (m == dim + 5) { return "S1_xy"; }

      else if (m == dim + 6) { return "S2_x"; }
      else if (m == dim + 7) { return "S2_y"; }
      else if (m == dim + 8) { return "S2_xy"; }

      else if (m == dim + 9) { return "S3_x"; }
      else if (m == dim + 10) { return "S3_y"; }
      else if (m == dim + 11) { return "S3_xy"; }

      else if (m == dim + 12) { return "S4_x"; }
      else if (m == dim + 13) { return "S4_y"; }
      else if (m == dim + 14) { return "S4_xy"; }

      else return "tooMany";
    } else Exit("dim = 3 ToDo")

    return "NotImplementedDimensionName";
  }*/
};

class TProblemElasticity : public ProblemBase {
  int verbose;
  double rho;
  double lambda;
  double mu;
  int dim;
public:
  TProblemElasticity(int d = 3) : ProblemBase(),
                                  dim(d), verbose(0), mu(1), rho(1), lambda(1) {
    Config::Get("lambda", lambda);
    Config::Get("mu", mu);
    Config::Get("rho", rho);
  }

  std::string Name() const override {
    return "TProblemElasticity";
  }

  virtual bool HasExactSolution() const { return false; }

  virtual double Rho(const Point &x) const { return rho; }

  virtual double Lambda(const Point &x) const { return lambda; }

  virtual double Mu(const Point &x) const { return mu; }

  virtual double Kappa(const Point &x) const {Exit("Only Acoustic Wave"); }

  double C_s(const Point &x) const { return sqrt(Mu(x) / Rho(x)); }

  double C_p(const Point &x) const { return sqrt((2.0 * Mu(x) + Lambda(x)) / Rho(x)); }

  virtual bool HasRHS() const { return false; }

  virtual double F(const Point &x, COMPONENT comp) const { return 0.0; }

  VectorField F_Velocity(const Point &x) const {
    return VectorField(F(x, COMPONENT::V_X), F(x, COMPONENT::V_Y));
  }

  Tensor F_Stress(const Point &x) const {
    return Tensor(F(x, COMPONENT::S0_x), F(x, COMPONENT::S0_xy),
                  F(x, COMPONENT::S0_xy), F(x, COMPONENT::S0_y));
  }

  virtual double ut(const Point &x, COMPONENT comp) const {Exit("Not implemented"); }

  virtual double ut(const Point &x,
                    const Cell &c,
                    COMPONENT comp) const {Exit("Not implemented"); }

  VectorField ut_Velocity(const Point &x) const {
    return VectorField(ut(x, COMPONENT::V_X), ut(x, COMPONENT::V_Y));
  }

  Tensor ut_Stress(const Point &x) const {
    return Tensor(ut(x, COMPONENT::S0_x), ut(x, COMPONENT::S0_xy),
                  ut(x, COMPONENT::S0_xy), ut(x, COMPONENT::S0_y));
  }

  virtual double ut_dual(const Point &,
                         COMPONENT comp,
                         const Point &,
                         const Point &) const { return 0; }

  virtual VectorField ut_dual_Velocity(const Point &x,
                                       const Point &a,
                                       const Point &b) const {
    return VectorField(ut_dual(x, COMPONENT::V_X, a, b),
                       ut_dual(x, COMPONENT::V_Y, a, b));
  }

  virtual Tensor ut_dual_Stress(const Point &x,
                                const Point &a,
                                const Point &b) const {
    return Tensor(ut_dual(x, COMPONENT::S0_x, a, b),
                  ut_dual(x, COMPONENT::S0_xy, a, b),
                  ut_dual(x, COMPONENT::S0_xy, a, b),
                  ut_dual(x, COMPONENT::S0_y, a, b));
  }

  virtual bool Dirichlet(const Point &z) const {Exit("Not implemented"); }

  virtual void SetDirichlet(Vector &u) const {Exit("Not implemented"); }

  int Dim() const {
    if (dim == 2) return 5;
    Exit("Not implemented");
  }

  std::vector<COMPONENT> GetComponents() const override {
    if (dim != 2) Exit("Not implemented.")
    return std::vector{COMPONENT::S0_x, COMPONENT::S0_y, COMPONENT::S0_xy,
                       COMPONENT::V_X, COMPONENT::V_Y};
  }

/*
  std::string dimensionName(int m) const override {
    if (dim == 2) {
      if (m == 0) { return "S_x"; }
      else if (m == 1) { return "S_y"; }
      else if (m == 2) { return "S_xy"; }
      else if (m == 3) { return "V_x"; }
      else if (m == 4) { return "V_y"; }
    }
    return "NotImplementedDimensionName";
  }*/
};

TProblemElasticity *CreateElasticProblem(const string &name);

std::unique_ptr<TProblemElasticity> CreateElasticProblemUnique(const string &name);

TProblemViscoElastic *CreateViscoElasticProblem(const string &name);

std::unique_ptr<TProblemViscoElastic> CreateViscoElasticProblemUnique(const string &name);

#endif //VISCOELASTICPROBLEMS_HPP
