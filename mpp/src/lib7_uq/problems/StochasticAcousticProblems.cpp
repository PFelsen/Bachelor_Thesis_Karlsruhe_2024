#include "StochasticAcousticProblems.hpp"


class GaussHatAndRicker2D : public IStochasticAcousticProblem {
protected:
  double startTime = 0.0;

  double factor = 10.0;

  Point shotLocation;

public:
  explicit GaussHatAndRicker2D(Point shotLocation = {0.5, 0.75}) :
      IProblem("Square"), IStochasticAcousticProblem(2, 0), shotLocation(shotLocation) {}

  double ForceP_i(double t, const Cell &c, const Point &z, int n) const override {
    if (n == 0) {
      return factor * Ricker(t - startTime, 0.1) * GaussHat(dist(z, shotLocation), 0.1);
    }
    return 0.0;
  }

  bool InRegionOfInterest(const Point &x) const override {
    if (0.25 <= x[0] && x[0] <= 0.75 && x[1] <= 0.25) return true;
    else return false;
  }

  double ut(double t, const Point &x, const Cell &c, int i) const override { return 0.0; }

  bool HasRHS() const override { return true; }
};

class StochasticGaussHatAndRicker2D : public GaussHatAndRicker2D {
  CirculantEmbedding2D circEmb;
public:
  StochasticGaussHatAndRicker2D() :
      IProblem("Square"),
      GaussHatAndRicker2D(Point{0.5, 0.75}), circEmb(CirculantEmbedding2D()) {}

  void drawSample(const SampleID &id) override { }

  double Rho(const Cell &c, const Point &z) const override {
    return 0.0;
//    return this->circEmb.EvalSample(c.Center());
  }
};

class StochasticRHSGaussHatAndRicker2D : public GaussHatAndRicker2D {
  KLExpansionGenerator KLE;
public:
  StochasticRHSGaussHatAndRicker2D() :
      IProblem("Square"),

      GaussHatAndRicker2D(Point{0.5, 0.75}), KLE(KLExpansionGenerator()) {}

  void drawSample(const SampleID &id) override { KLE.DrawSample(id, meshes); }

  double ForceP_i(double t, const Cell &c, const Point &z, int n) const override {
    if (n == 0) {
      return factor * KLE.EvalSample()[t] * GaussHat(dist(z, shotLocation), 0.1);
    }
    return 0.0;
  }
};

class RiemannWave2D : public IStochasticAcousticProblem {
  VectorField normal = VectorField(1.1, 0.7);
  double kappa_L = 1;
  double kappa_R = 1;
  double rho_L = 1;
  double rho_R = 1;
  VectorField v_L = VectorField(0.0, 0.0);
  VectorField v_R = VectorField(0.0, 0.0);
  double p_L = 1;
  double p_R = -1;
  double c_L;
  double c_R;
  double Z_L;
  double Z_R;
  double beta_L;
  double beta_R;
public:
  RiemannWave2D() :
      IProblem("Todo"),
      IStochasticAcousticProblem(2, 0) {
    Config::Get("kappa_L", kappa_L);
    Config::Get("kappa_R", kappa_R);
    Config::Get("rho_L", rho_L);
    Config::Get("rho_R", rho_R);
    Config::Get("normal_x", normal[0]);
    Config::Get("normal_y", normal[1]);
    normal /= norm(normal);
    c_L = sqrt(kappa_L / rho_L);
    c_R = sqrt(kappa_R / rho_R);
    Z_L = rho_L * c_L;
    Z_R = rho_R * c_R;
    beta_L = (p_R - p_L + Z_R * normal * (v_R - v_L)) / (Z_L + Z_R);
    beta_R = (p_R - p_L - Z_L * normal * (v_R - v_L)) / (Z_L + Z_R);
  }

  bool HasRHS() const override { return false; }

  double Kappa_i(const Cell &c, const Point &z, int i) const override {
    if (i != 0) { THROW("Kappa_i called with i != 0")}
    if (normal * z < 0) return kappa_L;
    return kappa_R;
  }

  double Rho(const Cell &c, const Point &z) const override {
    if (normal * z < 0) return rho_L;
    return rho_R;
  }

  bool HasExactSolution() const override { return true; }

  double ut(double t, const Point &z, const Cell &c, int i) const override {
    double x = z * normal;
    if (x <= -c_L * t) {
      switch (i) {
        case 0:
          return v_L[0];
        case 1:
          return v_L[1];
        case 2:
          return p_L;
        default: Exit("Component not implemented");
      }
    } else if (x < 0) {
      switch (i) {
        case 0:
          return v_L[0] + beta_L * normal[0];
        case 1:
          return v_L[1] + beta_L * normal[1];
        case 2:
          return p_L + beta_L * Z_L;
        default: Exit("Component not implemented");
      }
    } else if (x < c_R * t) {
      switch (i) {
        case 0:
          return v_R[0] + beta_R * normal[0];
        case 1:
          return v_R[1] + beta_R * normal[1];
        case 2:
          return p_R - beta_R * Z_R;
        default: Exit("Component not implemented");
      }
    } else {
      switch (i) {
        case 0:
          return v_R[0];
        case 1:
          return v_R[1];
        case 2:
          return p_R;
        default: Exit("Component not implemented");
      }
    }
  }
};

class Linear : public IStochasticAcousticProblem {
public:
  Linear() : IProblem("Square"),
             IStochasticAcousticProblem(2, 0) {}

  bool HasExactSolution() const override { return false; }

  double ut(double t, const Point &x, const Cell &c, int i) const override {
    switch (i) {
      case 0:
        return -t + x[0];
      case 1:
        return x[1];
      case 2:
        return t - x[0];
      default: Exit("Component not implemented");
    }
  }
};

class Quadratic : public IStochasticAcousticProblem {
public:
  Quadratic() :
      IProblem("Todo"),
      IStochasticAcousticProblem(2, 0) {}

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }

  double ut(double t, const Point &x, const Cell &c, int i) const override {
    switch (i) {
      case 0:
        return t * x[0];
      case 1:
        return x[1] * x[1];
      case 2:
        return x[0] * x[1];
      default: Exit("Component not implemented");
    }
  }

  VectorField ForceV(double t, const Cell &c, const Point &x) const override {
    return {x[0] - x[1], -x[0]};
  }

  double ForceP_i(double t, const Cell &c, const Point &x,int i) const override {
    return -t - 2 * x[1];
  }
};

class CRC : public IStochasticAcousticProblem {
protected:

  double aa = 0.5;

  double bb = 1.0;

  void drawSample(const SampleID &id) override {}

public:
  CRC() :
      IProblem("Todo"),

      IStochasticAcousticProblem(2, 0) {
    T = 8.0;
  }

  double GetStepSize(int level) const override { return 0.25 * pow(2, -level); }

  double ut(double t, const Point &x, const Cell &c, int i) const override {
    if (x[0] - bb > aa || x[0] < bb) return 0.0;
    double w = sqrt(aa * aa + bb * bb);
    switch (i) {
      case 0:
        return sin(w * Pi * t) * sin(2.0 * Pi * ((x[0] - bb) / aa) + Pi);
      case 1:
        return 0;
      case 2:
        return cos(w * Pi * t) * cos(2.0 * Pi * ((x[0] - bb) / aa) + Pi) + 1.0;
      default :
        return 0.0;
    }
  }

  int BndID(const Point &z) const override { return 2; }

  bool InRegionOfInterest(const Point &x) const override {
    if (4.5 <= x[0] && x[0] <= 5.5) return true;
    else return false;
  }
};

#if 0
class SGInitialValueCRC : public IStochasticAcousticProblem {
protected:
  SparseGridGenerator sparseGridGen;

  double startTime = 0.0;

  double factor = 10.0;

  Point center;

  Point shotLocation;

  void drawSample(const SampleID &id) override {
    RVector sample = 0.5 * sparseGridGen.DrawAndEvalSample(id);
    shotLocation = {center[0] + sample[0], center[1] + sample[1]};
  }

public:
  explicit SGInitialValueCRC(Point center = {1.0, 0.0}) :
    center(center), IProblem("deformed"),

    IStochasticAcousticProblem(2, 0),
    sparseGridGen(GridDomain({{-1, -1}, {1, 1}})) {
    T = 8.0;
  }

  void InitGenerator(int init) override {
    sparseGridGen.InitGenerator(meshes, init);
  }

  double SampleWeight() const override {
    return sparseGridGen.SampleWeight(CurrentID());
  }

  int NumOfSamples() const override {
    return sparseGridGen.GetNumPoints();
  }

  double GetStepSize(int level) const override { return 0.25 * pow(2, -level); }

  double ForceP_i(double t, const Cell &c, const Point &z, int n) const override {
    if (n == 0) {
      return factor * Ricker(t - startTime) * GaussHat(dist(z, shotLocation));
    }
    return 0.0;
  }

  bool InRegionOfInterest(const Point &x) const override {
    if (4.5 <= x[0] && x[0] <= 5.5) return true;
    else return false;
  }

  double ut(double t, const Point &x, const Cell &c, int i) const override { return 0.0; }

  bool HasRHS() const override { return true; }
};
#endif

IStochasticAcousticProblem *CreateStochasticAcousticProblem(const string &problemName) {
  IStochasticAcousticProblem *problem = nullptr;
#if 0
  if (problemName == "SGInitialValueCRC")
    problem = new SGInitialValueCRC();
#endif
  if (problemName == "StochasticGaussHatAndRicker2D")
    problem = new StochasticGaussHatAndRicker2D();
  if (problemName == "StochasticRHSGaussHatAndRicker2D")
    problem = new StochasticRHSGaussHatAndRicker2D();

  return problem;
}

std::unique_ptr<IStochasticAcousticProblem>
CreateStochasticAcousticProblemUnique(const string &problemName) {
  return std::unique_ptr<IStochasticAcousticProblem>(
      CreateStochasticAcousticProblem(problemName));
}

std::shared_ptr<IStochasticAcousticProblem>
CreateStochasticAcousticProblemShared(const std::string &problemName) {
  return std::shared_ptr<IStochasticAcousticProblem>(
      CreateStochasticAcousticProblem(problemName)
  );
}
