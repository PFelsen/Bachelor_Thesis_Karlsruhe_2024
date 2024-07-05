#ifndef MPP_ACOUSTIC_PROBLEMS_H
#define MPP_ACOUSTIC_PROBLEMS_H

#include "IProblem.hpp"
#include "RVector.hpp"
#include "SeismogramData.hpp"

double Ricker(double tau, double duration);

double CosHat(double dist, double width);

double GaussHat(double dist, double width);

// TODO maybe use RVector instead of DampingVector
class DampingVector : public std::valarray<double> {
public:
  DampingVector(const size_t size)
      : std::valarray<double>(0.0, std::max((size_t)1, size)) {}
};

inline double absNorm(const DampingVector &dampingVector) {
  return std::accumulate(std::begin(dampingVector), std::end(dampingVector),
                         0.0, [](double in1, double in2) {
                           return in1 + std::abs(in2);
                         });
}

inline double maxNorm(const DampingVector &dampingVector) {
  return *std::max_element(std::begin(dampingVector), std::end(dampingVector),
                           [](double in1, double in2) {
                             return std::abs(in1) < std::abs(in2);
                           });
}

class AcousticProblem : virtual public IProblem {
protected:
  int dim;

  int numL = 0;

  double cfl = 0.5;

  double t0 = 0.0;

  double T = 1.0;

  double rho = 1.0;
  double kappa = 1.0;
  double tau = 0.1 / M_PI;
  double tau_p = 0.1;
  std::vector<COMPONENT> components;

  static std::vector<COMPONENT> createComponents(int dim, int nL) {
    std::vector<COMPONENT> components{COMPONENT::V_X};
    if (dim > 1)
      components.push_back(COMPONENT::V_Y);
    if (dim > 2)
      components.push_back(COMPONENT::V_Z);
    components.push_back(COMPONENT::P0);
    if (nL > 0)
      components.push_back(COMPONENT::P1);
    if (nL > 1)
      components.push_back(COMPONENT::P2);
    if (nL > 2)
      components.push_back(COMPONENT::P3);
    if (nL > 3)
      components.push_back(COMPONENT::P4);
    if (nL > 4)
      components.push_back(COMPONENT::P5);
    if (nL > 5) {
      Exit("implement more damping variables.")
    }
    return components;
  }

public:
  AcousticProblem(int dim, int numL)
      : dim(dim), numL(numL), components(createComponents(dim, numL)) {
    Config::Get("CFL", cfl);
    Config::Get("t0", t0);
    Config::Get("T", T);
  }

  virtual double GetStepSize(int level) const { return cfl * pow(2, -level); }

  virtual double GetStartTime() const { return t0; }

  virtual double GetEndTime() const { return T; }

  double GetCFL() const { return cfl; }

  virtual int BndID(const Point &z) const { return 2; }

  int SpaceDim() const { return dim; }

  int NumDamping() const { return numL; }

  bool Damping() const { return numL > 0; }

  virtual double ut(double t, const Point &x, const Cell &c, int i) const {
    Exit("Not implemented: ut(double t, const Point &x, const Cell &c, int i)");
  }

  virtual double ut(double t, const Point &x, const Cell &c,
                    COMPONENT comp) const {
    return 0.0;
  }

  virtual VectorField ut_Velocity(double t, const Point &x,
                                  const Cell &c) const {
    if (dim == 2) {
      return VectorField{ut(t, x, c, COMPONENT::V_X),
                         ut(t, x, c, COMPONENT::V_Y)};
    } else if (dim == 1) {
      return VectorField{ut(t, x, c, COMPONENT::V_X)};
    } else if (dim == 3) {
      return VectorField{ut(t, x, c, COMPONENT::V_X),
                         ut(t, x, c, COMPONENT::V_Y),
                         ut(t, x, c, COMPONENT::V_Z)};
    }
    Exit("ut_Velocity not implemented for dim=" + std::to_string(dim))
  }

  virtual double ut_Pressure(double t, const Point &x, const Cell &c) const {
    return ut(t, x, c, COMPONENT::P0);  // ut_Pressure = sum_{i=0}^numL P_i?
  }

  virtual DampingVector ut_DampingPressure(double t, const Point &x,
                                           const Cell &c) const {
    DampingVector DP(numL);
    for (int i = 0; i < numL; i++) {
      COMPONENT comp = GetDampingComponent(i);
      DP[i] = ut(t, x, c, comp);
    }
    return DP;
  }

  virtual double ut_dual(const Point &x, COMPONENT comp, const Point &a,
                         const Point &b) const {
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1])
      return 0.0;
    Point ba = b - a;
    if (comp != COMPONENT::P0)
      return 0;
    else
      return 1.0 / (ba[0] * ba[1]);
  }

  double Dirichlet(double t, const Cell &c, const Point &z, int i) const {
    return ut(t, z, c, dim + i);
  }

  double Neumann(double t, const Cell &c, const Point &z,
                 const VectorField &N) const {
    VectorField V = zero;
    for (int i = 0; i < dim; ++i) {
      V[i] = ut(t, z, c, i);
    }
    return V * N;
  }

  virtual bool HasRHS() const { return false; }

  virtual bool hasPointSource(Point &x) const { return false; }

  virtual double F_PointSource(double t, const Point &x) const {
    Exit("Not implemented: F_PointSource(const Point &x)")
  }
  virtual double ForceP_i(double t, const Cell &c, const Point &x,
                          int i) const {
    return 0.0;
  }

  virtual VectorField ForceV(double t, const Cell &c, const Point &x) const {
    return zero;
  }

  virtual void Force(VectorField &V_rhs, RVector &Force_, double t,
                     const Cell &c, const Point &z) const {
    V_rhs = ForceV(t, c, z);
    for (int n = 0; n < Force_.size(); ++n) {
      Force_[n] = ForceP_i(t, c, z, n);
    }
  }

  virtual double Kappa(const Cell &c, const Point &x) const final {
      double kappa = 0.0;
      for (int i = 0; i < 1 + numL; ++i) {
        kappa += Kappa_i(c, x, i);
      }
      return kappa;
  }

  virtual double Kappa_i(const Cell &c, const Point &x, int i) const {
    if (i == 0) return kappa;
    Exit("not implemented here!")
  }

  virtual double Tau_i(const Cell &c, const Point &x, int i) const {
    return tau;
  }

  virtual double Tau_p(const Cell &, const Point &) const {
    Exit("Tau_p not implemented here!");
  }

  virtual double alpha(const Cell &, const Point &) const {
    Exit("alpha not implemented here!");
  }
  virtual double V_p(const Cell &, const Point &) const {
    Exit("V_p not implemented here!");
  }

  virtual double Rho(const Cell &c, const Point &x) const { return rho; }

  virtual bool InRegionOfInterest(const Point &x) const { return true; }

  virtual double u0(const Cell &c, const Point &x, int i) const {
    return ut(0, x, c, i);
  }

  virtual bool hasInitialData() const { return false; }

  virtual double Mdtut(double t, const Point &x, const Cell &c,
                       COMPONENT comp) const {
    Exit("not implemented for " + Name()) return 0.0;
  }

  virtual double Mdtut_Pressure(double t, const Point &x, const Cell &c) const {
    return Mdtut(t, x, c, COMPONENT::P0);
  }

  virtual VectorField Mdtut_Velocity(double t, const Point &x,
                                     const Cell &c) const {
    if (dim == 2) {
      return {Mdtut(t, x, c, COMPONENT::V_X), Mdtut(t, x, c, COMPONENT::V_Y)};
    } else if (dim == 1) {
      return {Mdtut(t, x, c, COMPONENT::V_X)};
    } else if (dim == 3) {
      return {Mdtut(t, x, c, COMPONENT::V_X), Mdtut(t, x, c, COMPONENT::V_Y),
              Mdtut(t, x, c, COMPONENT::V_Z)};
    }
    Exit("Mdtut_Velocity not implemented for dim=" + std::to_string(dim))
  }

  virtual double Dut(double t, const Point &x, const Cell &c,
                     COMPONENT comp) const {
    return 0.0;
  }

  virtual double Aut(double t, const Point &x, const Cell &c,
                     COMPONENT comp) const {
    Exit("not implemented for " + Name()) return 0.0;
  }

  virtual double Aut_Pressure(double t, const Point &x, const Cell &c) const {
    return Aut(t, x, c, COMPONENT::P0);
  }

  virtual VectorField Aut_Velocity(double t, const Point &x,
                                   const Cell &c) const {
    if (dim == 2) {
      return {Aut(t, x, c, COMPONENT::V_X), Aut(t, x, c, COMPONENT::V_Y)};
    } else if (dim == 1) {
      return {Aut(t, x, c, COMPONENT::V_X)};
    } else if (dim == 3) {
      return {Aut(t, x, c, COMPONENT::V_X), Aut(t, x, c, COMPONENT::V_Y),
              Aut(t, x, c, COMPONENT::V_Z)};
    }
    Exit("Aut_Velocity not implemented for dim=" + std::to_string(dim))
  }

  virtual double F(double t, const Cell &c, const Point &x, COMPONENT comp) const {
    if (HasExactSolution()) {
      return Mdtut(t, x, c, comp) + Aut(t, x, c, comp) + Dut(t, x, c, comp);
    }
    return 0.0;
  }

  virtual DampingVector F_DampingPressure(double t, const Cell &c,
                                          const Point &x) const {
    DampingVector DV(numL);
    for (int i = 0; i < numL; i++) {
      DV[i] = F(t, c, x, GetDampingComponent(i));
    }
    return DV;
  }

  virtual double F_Pressure(double t, const Cell &c, const Point &x) const {
    return F(t, c, x, COMPONENT::P0);
  }

  virtual VectorField F_Velocity(double t, const Cell &c,
                                 const Point &x) const {
    if (dim == 2) {
      return {F(t, c, x, COMPONENT::V_X), F(t, c, x, COMPONENT::V_Y)};
    } else if (dim == 1) {
      return {F(t, c, x, COMPONENT::V_X)};
    } else if (dim == 3) {
      return {F(t, c, x, COMPONENT::V_X), F(t, c, x, COMPONENT::V_Y),
              F(t, c, x, COMPONENT::V_Z)};
    }
    Exit("F_Velocity not implemented for dim=" + std::to_string(dim))
  }

  virtual double F_Pkt(const double &t, const Cell &c, const int &i) const {
    return 0.0;
  }

  int Dim() const { return dim + 1 + numL; }

  int nL() const { return numL; }

  std::vector<COMPONENT> GetComponents() const { return components; };

  virtual std::shared_ptr<SeismogramData> GetSourceData() const {
    return nullptr;
  }
};

AcousticProblem *CreateAcousticProblem(const string &name);

std::unique_ptr<AcousticProblem> CreateAcousticProblemUnique(const string &name);

std::shared_ptr<AcousticProblem> CreateAcousticProblemShared(const string &name);

std::shared_ptr<AcousticProblem> CreateAcousticProblemShared(const string &name, std::shared_ptr<Meshes> meshes);

struct FlexibleAcousticProblemData {
public:
  std::shared_ptr<Meshes> meshes;
  std::shared_ptr<AcousticProblem> problem;
  std::function<double(const Point &)> kappa = nullptr;
  std::function<double(const Point &)> rho = nullptr;
  std::shared_ptr<SeismogramData> sourceData;
  std::function<double(double t, const Cell &c, const Point &)> FPressure =
      nullptr;
  std::string name;
};

class FlexibleAcousticProblem : public AcousticProblem {
public:
  const FlexibleAcousticProblemData &builder;
  FlexibleAcousticProblem(FlexibleAcousticProblemData &builderParam)
      : IProblem(builderParam.meshes),
        AcousticProblem(builderParam.problem->Dim(),
                        builderParam.problem->nL()),
        builder(builderParam) {}

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    double kappa = builder.kappa ? builder.kappa(x) : 1.0;
    if (i == 0) {
      return kappa;
    }
    return Tau_p(c, x) * kappa;
  }

  double Rho(const Cell &c, const Point &x) const override {
    return builder.rho ? builder.rho(c()) : 1.0;
  }

  double F(double t, const Cell &c, const Point &x, COMPONENT comp) const override {
    if (comp == COMPONENT::P0 && builder.FPressure) {
      return builder.FPressure(t, c, x);
    }
    return 0.0;
  }

  bool HasRHS() const override { return builder.FPressure != nullptr; }

  std::shared_ptr<SeismogramData> GetSourceData() const override {
    return builder.sourceData ? builder.sourceData : nullptr;
  }

  std::string Name() const override {
    return builder.name.empty() ? "FlexibleAcousticProblem without name"
                                : builder.name;
  }
};

#endif  // MPP_ACOUSTIC_PROBLEMS_H
