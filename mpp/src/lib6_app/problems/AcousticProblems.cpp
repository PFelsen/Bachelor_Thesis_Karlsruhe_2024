#include "AcousticProblems.hpp"
#include "ParameterImage.hpp"
#include "CImg.hpp"

double Ricker(double tau, double duration) {
  double f2 = 1.0 / (duration * duration);
  double exponent = M_PI * M_PI * f2 * tau * tau;
  return (1.0 - 2.0 * exponent) * exp(-exponent);
}

double CosHat(double dist, double width) {
  if (dist >= width) return 0.0;
  return pow(cos(M_PI / 2.0 * dist / width), 6);
}

double GaussHat(double dist, double width) {
  if (dist >= width) return 0;
  dist /= width;
  return exp(-1 / (1 - dist * dist)) / 0.466512050417577 * (1.0 / (width * width));
}

class AcousticBenchmarkC2 : public AcousticProblem {
  //Timestepping
  double v_bg = 3500.0;
  double rho_bg = 2000.0;
  double tau_p_bg = 0.02;
  double impedance_bg = 1;
  double size_perturb = 200.0;
  double x_perturb = 900.0;
  double y_perturb_0 = 1200.0;
  double y_perturb_1 = 600.0;
  double y_perturb_2 = -200.0;
  double perturb_0 = 9;
  double perturb_1 = 20;
  double perturb_2 = 10;
  // Force function parameters (Ricker source)
  double startTime = 0.2; //3*M_PI/w_0
  double factor = 0.01;
  double x_PointSource = 375.0;
  double y_PointSource = 375.0;
  double freq_source = 25; // w_0/2*M_PI
  double w_0 = freq_source * 2 * M_PI; // w_0/2*M_PI
  int Qfactor = 20;

  std::vector<double> tau_i;
  double alphaValue = 0.0;

public :
  static int GetDampingFromConfig() {
    int numL = 0;
    Config::Get("numL", numL);
    return numL;
  }

  explicit AcousticBenchmarkC2(int numL) : IProblem("SeismogramGrid"), AcousticProblem(2, numL) {
    Config::Get("x_perturb", x_perturb);
    Config::Get("y_perturb_0", y_perturb_0);
    Config::Get("y_perturb_1", y_perturb_1);
    Config::Get("y_perturb_2", y_perturb_2);
    Config::Get("y_perturb_0_value", perturb_0);
    Config::Get("y_perturb_1_value", perturb_1);
    Config::Get("y_perturb_2_value", perturb_2);
    Config::Get("size_perturb", size_perturb);

    // Parameters for background domain
    Config::Get("rho_bg", rho_bg);
    Config::Get("v_bg", v_bg);
    Config::Get("tau_p_bg", tau_p_bg);
    Config::Get("impedance_bg", impedance_bg);
    Config::Get("Qfactor", Qfactor);

    std::vector<double> relax_freq;
    if (numL == 1) {
      relax_freq = {18.3294};
      tau_p_bg = 0.1278;
    } else if (numL == 2) {
      relax_freq = {2.8288, 37.5512};
      tau_p_bg = 0.1016;
    } else if (numL == 3) {
      relax_freq = {1.4618, 13.3875, 84.4705};
      tau_p_bg = 0.0779;
    } else if (numL == 4) {
      relax_freq = {0.8837, 5.1206, 32.3769, 146.9207};
      tau_p_bg = 0.0661;
    } else if (numL == 5) {
      relax_freq = {0.5110, 2.0534, 10.0022, 45.1536, 126.9230};
      tau_p_bg = 0.0552;
    }
    for (double f: relax_freq) {
      tau_i.push_back(1.0 / (f * 2.0 * M_PI));
    }
    calculateAlphaValue();
  }

  void calculateAlphaValue() {
    for (int l = 0; l < numL; ++l) {
      alphaValue += (w_0 * w_0 * tau_i[l] * tau_i[l]) / (1 + (w_0 * w_0 * tau_i[l] * tau_i[l]));
    }
  }

  std::string Name() const override { return "AcousticBenchmarkC2L" + std::to_string(numL); }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return false; }

  bool hasInitialData() const override { return false; }

  bool hasPointSource(Point &x) const override {
    x = Point(x_PointSource, y_PointSource);
    return false;
  }

  // Checking the perturbed domain
  bool isInPerturb(const Point &x, double y_perturb) const {
    if (x[0] >= x_perturb && x[0] <= x_perturb + size_perturb) {
      if (x[1] >= y_perturb && x[1] <= y_perturb + size_perturb) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  // Allocation of Rho to the perturbed domain
  double Rho(const Cell &c, const Point &x) const override {
    if (isInPerturb(x,
                    y_perturb_0)) {          // y_perturb_0 mit flags oder block number für die Identifizierung
      return rho_bg * perturb_0;
    } else if (isInPerturb(x, y_perturb_1)) {
      return rho_bg * perturb_1;
    } else if (isInPerturb(x, y_perturb_2)) {
      return rho_bg * perturb_2;
    }
    return rho_bg * impedance_bg;
  }

  double Tau_p(const Cell &c, const Point &x) const override {
    if (isInPerturb(x, y_perturb_0)) {
      return tau_p_bg * perturb_0;
    } else if (isInPerturb(x, y_perturb_1)) {
      return tau_p_bg * perturb_1;
    } else if (isInPerturb(x, y_perturb_2)) {
      return tau_p_bg * perturb_2;
    }
    return tau_p_bg;
  }

  double V_p(const Cell &c, const Point &x) const override {
    if (isInPerturb(x, y_perturb_0)) {
      return v_bg * perturb_0;
    } else if (isInPerturb(x, y_perturb_1)) {
      return v_bg * perturb_1;
    } else if (isInPerturb(x, y_perturb_2)) {
      return v_bg * perturb_2;
    }
    return v_bg;
  }

  double Tau_i(const Cell &c, const Point &x, int i) const override {
    return tau_i[i - 1];
  }

  double alpha(const Cell &c, const Point &x) const override {
    return alphaValue;
  }

  double Kappa_0(const Cell &c, const Point &x) const {
    double v = V_p(c, x);
    double kappa_0 = ((numL == 0) ? Rho(c, x) * v * v : Rho(c, x) * v * v /
                                                        (1 + alpha(c, x) * Tau_p(c, x)));
    return kappa_0;
  }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    return (i == 0) ? Kappa_0(c, x) : Tau_p(c, x) * Kappa_0(c, x);
  }

  double ForceP_i(double t, const Cell &c, const Point &x, int i) const override {
    switch (i) {
      case 0:
        if (t > 0.1 && t <= 0.34 && 100 <= x[0] && x[0] <= 300 && 400 <= x[1] && x[1] <= 600)
          return 0.0001;
        else return 0.0;
      case 1:
        return 0.0;
      case 2:
        return 0;
      case 3:
        return 0;
      case 4:
        return 0;
      case 5:
        return 0;
      default: Exit("Component not implemented");
    }
  }

  double F(double t, const Cell &c, const Point &x, COMPONENT comp) const override {
    if (c.IsSpaceTime() && t != x.t()) THROW("Check Time Component for F")

    switch (comp) {
      case COMPONENT::V_X:
        return 0.0;
      case COMPONENT::V_Y:
        return 0.0;
      case COMPONENT::P0:
        if (t > 0.1 && t <= 0.3 && 100 <= x[0] && x[0] <= 300 && 400 <= x[1] && x[1] <= 600)
          return 0.0001;
        else return 0.0;
      case COMPONENT::P1:
        //        if (500 <= x[0] && x[0] <= 700 && 800 <= x[1] && x[1] <= 1000) return 1.0;
        //        else return 0.0; //Ungleiche Lösung bei unterschiedlichen P_i
        return 0.0;
      case COMPONENT::P2:
        return 0;
      case COMPONENT::P3:
        return 0;
      case COMPONENT::P4:
        return 0;
      case COMPONENT::P5:
        return 0;
      default: Exit("Component not implemented");
    }
  }

  double ut(double t, const Point &x, const Cell &c, int i) const override {
    return 0.0;
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override { return 0.0; }
};

typedef cimg_library::CImg<unsigned char> MyImg;

class AcousticWaveDamping : public AcousticProblem {
  double rho = 10.0;
  double kappa_0 = 100.0;
  double tau_p = 1.05;
  std::vector<double> tau_i = {0.6, 0.3, 0.1, 0.05, 0.02};
public:
  static int GetDampingFromConfig() {
    int numL = 0;
    Config::Get("numL", numL);
    return numL;
  }

  explicit AcousticWaveDamping(int numL) : IProblem("UnitSquare"), AcousticProblem(2, numL) {
  }

  std::string Name() const override { return "AcousticWaveDamping"; }

  bool HasExactSolution() const override { return true; }

  bool HasRHS() const override { return true; }

  bool hasInitialData() const override { return true; }

  double Rho(const Cell &c, const Point &x) const override {
    return rho;
  }

  double Tau_i(const Cell &c, const Point &x, int i) const override {
    return tau_i[i - 1];
  }

  double Tau_p(const Cell &c, const Point &x) const override {
    return tau_p;
  }

  double Kappa_0(const Cell &c, const Point &x) const {
    return kappa_0;
  }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    return (i == 0) ? Kappa_0(c, x) : Tau_p(c, x) * Kappa_0(c, x);
  }

  VectorField ForceV(double t, const Cell &c, const Point &x) const override {
    if (numL == 0) {
      return {Rho(c, x), -Rho(c, x)};
    } else if (numL == 1) {
      return {Rho(c, x), -Rho(c, x) - 12};
    } else {
      return {Rho(c, x) + 1, -Rho(c, x) - 13};
    }
  }

  double ForceP_i(double t, const Cell &c, const Point &x, int i) const override {
    switch (i) {
    case 0:
      return 50 * (1 / Kappa_i(c, x, i)) - 19;
    case 1:
      return 14 / Kappa_i(c, x, i) + (14 * t + 12 * x[1]) / (Tau_i(c, x, i) * Kappa_i(c, x, i)) -
             19;
    case 2:
      return (-x[0] + x[1]) / (Tau_i(c, x, i) * Kappa_i(c, x, i)) - 19;
    case 3:
      return 1 / (Tau_i(c, x, i) * Kappa_i(c, x, i)) - 19;
    case 4:
      return 0.05 / (Tau_i(c, x, i) * Kappa_i(c, x, i)) - 19;
    case 5:
      return 0.01 / (Tau_i(c, x, i) * Kappa_i(c, x, i)) - 19;
    default: Exit("Component not implemented");
    }
  }

  double F(double t, const Cell &c, const Point &x,
           COMPONENT comp) const override {
    if (c.IsSpaceTime() && t != x.t()) THROW("Check Time Component for F")

    switch (comp) {
    case COMPONENT::V_X:
      if (numL == 0 || numL == 1) {
        return Rho(c, x);
      } else {
        return Rho(c, x) + 1;
      }
    case COMPONENT::V_Y:
      if (numL == 0) {
        return -Rho(c, x);
      } else if (numL == 1) {
        return -Rho(c, x) - 12;
      } else {
        return -Rho(c, x) - 13;
      }
    case COMPONENT::P0:
      return 50 * (1 / Kappa_i(c, x, 0)) - 19;
    case COMPONENT::P1:
      return 14 / Kappa_i(c, x, 1) + (14 * t + 12 * x[1]) / (Tau_i(c, x, 1) * Kappa_i(c, x, 1)) -
             19;
    case COMPONENT::P2:
      return (-x[0] + x[1]) / (Tau_i(c, x, 2) * Kappa_i(c, x, 2)) - 19;
    case COMPONENT::P3:
      return 1 / (Tau_i(c, x, 3) * Kappa_i(c, x, 3)) - 19;
    case COMPONENT::P4:
      return 0.05 / (Tau_i(c, x, 4) * Kappa_i(c, x, 4)) - 19;
    case COMPONENT::P5:
      return 0.01 / (Tau_i(c, x, 5) * Kappa_i(c, x, 5)) - 19;
    default: Exit("Component not implemented");
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

  double ut(double t, const Point &x, const Cell &c, int i) const override {
    return ut(t, x, c, components[i]);
  }

  double ut(double t, const Point &xx, const Cell &c, COMPONENT comp) const override {

    if (c.IsSpaceTime() && t != xx.t()) THROW("Check Time Component for ut")

    double x = xx[0];
    double y = xx[1];

    double vx = 17 * x + y + t;
    double vy = 2 * y + x - t;

    double p0 = 50 * t;
    double p1 = 14 * t + 12 * y;
    double p2 = -x + y;
    double p3 = 1;
    double p4 = 0.05;
    double p5 = 0.01;

    switch (comp) {
    case COMPONENT::V_X:
      return vx;
    case COMPONENT::V_Y:
      return vy;
    case COMPONENT::P0:
      return p0;
    case COMPONENT::P1:
      return p1;
    case COMPONENT::P2:
      return p2;
    case COMPONENT::P3:
      return p3;
    case COMPONENT::P4:
      return p4;
    case COMPONENT::P5:
      return p5;

    default: Exit("No.");
    }
  }
};
class AcousticWaveProcTest : public AcousticProblem {
  double rho = 1;
  double kappa_0 = 3000*3000;
  double tau_p = 1.05;
  std::vector<double> tau_i = {0.6, 0.3, 0.1, 0.05, 0.02};
public:
  static int GetDampingFromConfig() {
    int numL = 0;
    Config::Get("numL", numL);
    return numL;
  }
  explicit AcousticWaveProcTest(int numL) : IProblem("ST_large_rectangle"), AcousticProblem(2, numL) {
  }

  std::string Name() const override { return "AcousticWaveProcTest"; }

  bool HasExactSolution() const override { return false; }

  bool HasRHS() const override { return true; }

  bool hasInitialData() const override { return true; }

  double Rho(const Cell &c, const Point &x) const override {
    return rho;
  }

  double Tau_i(const Cell &c, const Point &x, int i) const override {
    return tau_i[i - 1];
  }

  double Tau_p(const Cell &c, const Point &x) const override {
    return tau_p;
  }

  double Kappa_0(const Cell &c, const Point &x) const {
    return kappa_0;
  }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    return (i == 0) ? Kappa_0(c, x) : Tau_p(c, x) * Kappa_0(c, x);
  }
  VectorField ForceV(double t, const Cell &c, const Point &x) const override {
    return {0, 0};
  }

  double F(double t, const Cell &c, const Point &x, COMPONENT comp) const override {
    if (c.IsSpaceTime() && t != x.t()) THROW("Check Time Component for F")

    switch (comp) {
    case COMPONENT::V_X:
      return 0.0;
    case COMPONENT::V_Y:
      return 0.0;
    case COMPONENT::P0:
      if (t > 0.1 && t <= 0.3 && -100 <= x[0] && x[0] <= 1000 && -650 <= x[1] && x[1] <= 500)
        return 100;
      else return 0.0;
    case COMPONENT::P1:
      return 0;
    case COMPONENT::P2:
      return 0;
    case COMPONENT::P3:
      return 0;
    case COMPONENT::P4:
      return 0;
    case COMPONENT::P5:
      return 0;
    default: Exit("Component not implemented");
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }
};
class TestProblem3DST : public AcousticProblem {
public:
  explicit TestProblem3DST(bool dirichlet = false) : AcousticProblem(3, 0) {}

  std::string Name() const override {
    return "TestProblem3DST";
  }

  int BndID(const Point &z) const override { return 1; }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return 1;
      case COMPONENT::V_X:
        return 2;
      case COMPONENT::V_Y:
        return 3;
      case COMPONENT::V_Z:
        return 4;
      default:
        return 0.0;
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }
};


class GaussHatAndRicker2D : public AcousticProblem {
protected:
  double startTime = 0.0;

  double factor = 10.0;

  Point shotLocation;

public:
  explicit GaussHatAndRicker2D(Point shotLocation = {0.5, 0.75}) :
      IProblem("Square"), AcousticProblem(2, 0), shotLocation(shotLocation) {}

  std::string Name() const { return "GaussHatAndRicker2D"; }

  double ForceP_i(double t, const Cell &c, const Point &z, int n) const override {
    if (n == 0) {
      return factor * Ricker(t - startTime, 0.1) * GaussHat(dist(z, shotLocation),0.1);
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

class RiemannWave2D : public AcousticProblem {
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
  RiemannWave2D() : IProblem("QD"), AcousticProblem(2, 0) {
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

  std::string Name() const override { return "RiemannWave2D"; }

  bool HasRHS() const override { return false; }

  double Kappa_i(const Cell &c, const Point &z, int i) const override {
    if (i == 0) {
      if (normal * z < 0) return kappa_L;
      return kappa_R;
    }
    THROW("Do not get here!")
  }

  double Rho(const Cell &c, const Point &z) const override {
    if (normal * z < 0) return rho_L;
    return rho_R;
  }

  bool HasExactSolution() const override { return true; }

  double ut(double t, const Point &z, const Cell &c, int i) const override {
    double x = z * normal;
    if (x <= -c_L * t) {
      switch(i) {
        case 0:
          return v_L[0];
        case 1:
          return v_L[1];
        case 2:
          return p_L;
        default: Exit("Component not implemented");
      }
    } else if (x < 0) {
      switch(i) {
        case 0:
          return v_L[0] + beta_L * normal[0];
        case 1:
          return v_L[1] + beta_L * normal[1];
        case 2:
          return p_L + beta_L * Z_L;
        default: Exit("Component not implemented");
      }
    } else if (x < c_R * t) {
      switch(i) {
        case 0:
          return v_R[0] + beta_R * normal[0];
        case 1:
          return v_R[1] + beta_R * normal[1];
        case 2:
          return p_R - beta_R * Z_R;
        default: Exit("Component not implemented");
      }
    } else {
      switch(i) {
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

class Linear : public AcousticProblem {
public:
  Linear() : IProblem("Square"), AcousticProblem(2, 0) {}

  std::string Name() const override { return "Linear"; }

  bool HasExactSolution() const override { return false; }

  double ut(double t, const Point &x, const Cell &c, int i) const override {
    switch(i) {
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

class Quadratic : public AcousticProblem {
public:
  Quadratic() : IProblem("Square"), AcousticProblem(2, 0) {}

  std::string Name() const { return "Quadratic"; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }

  double ut(double t, const Point &x, const Cell &c, int i) const override {
    switch(i) {
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

  double ForceP_i(double t, const Cell &c, const Point &x, int) const override {
    return -t - 2 * x[1];
  }
};

class CRC : public AcousticProblem {
protected:

  double aa = 0.5;

  double bb = 1.0;

public:
  CRC() : IProblem("deformed"), AcousticProblem(2, 0) {
    T = 8.0;
  }

  std::string Name() const override { return "CRC"; }

  double GetStepSize(int level) const override { return 0.25 * pow(2, -level); }

  double ut(double t, const Point &x, const Cell &c, int i) const override {
    if (x[0] - bb > aa || x[0] < bb) return 0.0;
    double w = sqrt(aa * aa + bb * bb);
    switch(i) {
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



class SphericalProblem3DST : public AcousticProblem {
public:
  const bool dirichlet;

  explicit SphericalProblem3DST(bool dirichlet = false) : AcousticProblem(3, 0),
                                                          dirichlet(dirichlet) {}

  std::string Name() const override {
    return "SphericalProblem3D";
  }

  int BndID(const Point &z) const override { return dirichlet ? 1 : 2; }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return false; }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (x.t() != 0 || comp != COMPONENT::P0) {
      return 0.0;
    }
    double d = dist(x, {0.5, 0.5, 0.5});
    return CosHat(d, 0.1);
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }
};


class PolynomialProblem1D : public AcousticProblem {
public:
  const int degree;
  const bool dirichlet;

  explicit PolynomialProblem1D(int degree, bool dirichlet = false) : AcousticProblem(1, 0),
                                                                     degree(degree),
                                                                     dirichlet(dirichlet) {}

  std::string Name() const override {
    return "Polynomial with degree " + std::to_string(degree);
  }

  int BndID(const Point &z) const override { return dirichlet ? 1 : 2; }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    return 1.0;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return pow(x[0], degree) + pow(x.t(), degree);
      case COMPONENT::V_X:
        return pow(x[0], degree) + pow(x.t(), degree);
      default: Exit("Component " + to_string(comp) + " not defined for " + Name());
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (degree == 0) {
      return 0.0;
    }
    return degree * pow(x.t(), degree - 1);
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (degree == 0) {
      return 0.0;
    }
    switch (comp) {
      case COMPONENT::P0:
        return -degree * pow(x[0], degree - 1);
      case COMPONENT::V_X:
        return -degree * pow(x[0], degree - 1);
      default: Exit("do not get here!");
        return 0.0;
    }
  }
};

class Riemann1DST : public AcousticProblem {
public:
  Riemann1DST() : AcousticProblem(1, 0) {}

  std::string Name() const override {
    return "Riemann1DST";
  }

  int BndID(const Point &z) const override { return 1; }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return false; }

  bool HasExactSolution() const override { return false; }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (comp == COMPONENT::P0)
      if (x[0] < 1 / 3.0) return 1;
    return 0;
  }
};

class PolynomialProblem3D : public AcousticProblem {
public:
  const int degree;
  const bool dirichlet;

  explicit PolynomialProblem3D(int degree, bool dirichlet = false) : AcousticProblem(3, 0),
                                                                     degree(degree),
                                                                     dirichlet(dirichlet) {}

  std::string Name() const override {
    return "Polynomial with degree " + std::to_string(degree);
  }

  int BndID(const Point &z) const override { return dirichlet ? 1 : 2; }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }

  double Rho(const Cell &c, const Point &x) const override {
    return 1.0;
  }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    return 1.0;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x[2], degree) + pow(x.t(), degree);
      case COMPONENT::V_X:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x[2], degree) + pow(x.t(), degree);
      case COMPONENT::V_Y:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x[2], degree) + pow(x.t(), degree);
      case COMPONENT::V_Z:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x[2], degree) + pow(x.t(), degree);
      default:
        return 0.0;
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (degree == 0) {
      return 0.0;
    }
    return degree * pow(x.t(), degree - 1);
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (degree == 0) {
      return 0.0;
    }
    switch (comp) {
      case COMPONENT::P0:
        return -degree * (pow(x[0], degree - 1) + pow(x[1], degree - 1) + pow(x[2], degree - 1));
      case COMPONENT::V_X:
        return -degree * pow(x[0], degree - 1);
      case COMPONENT::V_Y:
        return -degree * pow(x[1], degree - 1);
      case COMPONENT::V_Z:
        return -degree * pow(x[2], degree - 1);
      default:
        return 0.0;
    }
  }
};

class PolynomialProblem2D : public AcousticProblem {
public:
  const int degree;
  const bool dirichlet;

  explicit PolynomialProblem2D(int degree, bool dirichlet = false) : AcousticProblem(2, 0),
                                                                     degree(degree),
                                                                     dirichlet(dirichlet) {}

  std::string Name() const override {
    return "Polynomial with degree " + std::to_string(degree);
  }

  int BndID(const Point &z) const override { return dirichlet ? 1 : 2; }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    return 1.0;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x.t(), degree);
      case COMPONENT::V_X:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x.t(), degree);
      case COMPONENT::V_Y:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x.t(), degree);
      default:
        return 0.0;
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (degree == 0) {
      return 0.0;
    }
    return degree * pow(x.t(), degree - 1);
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (degree == 0) {
      return 0.0;
    }
    switch (comp) {
      case COMPONENT::P0:
        return -degree * (pow(x[0], degree - 1) + pow(x[1], degree - 1));
      case COMPONENT::V_X:
        return -degree * pow(x[0], degree - 1);
      case COMPONENT::V_Y:
        return -degree * pow(x[1], degree - 1);
      default:
        return 0.0;
    }
  }
};

// Beispiel mit numL
class PolynomialProblemNumL : public AcousticProblem {
public:
  const int degree;
  const bool dirichlet;

  explicit PolynomialProblemNumL(int degree, int numL, bool dirichlet = false)
      : AcousticProblem(2, numL),
        degree(degree),
        dirichlet(dirichlet) {
    tau_p = 1.0;
    tau = 1.0;
    kappa = 1.0;
    rho = 1.0;
  }

  std::string Name() const override {
    return "Polynomial with degree " + std::to_string(degree);
  }

  int BndID(const Point &z) const override { return dirichlet ? 1 : 2; }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    return 1.0;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x.t(), degree);
      case COMPONENT::P1:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x.t(), degree);
      case COMPONENT::P2:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x.t(), degree);
      case COMPONENT::P3:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x.t(), degree);
      case COMPONENT::V_X:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x.t(), degree);
      case COMPONENT::V_Y:
        return pow(x[0], degree) + pow(x[1], degree) + pow(x.t(), degree);
      default:
        return 0.0;
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (degree == 0) {
      return 0.0;
    }
    return degree * pow(x.t(), degree - 1);
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (degree == 0) {
      return 0.0;
    }
    switch (comp) {
      case COMPONENT::P0:
      case COMPONENT::P1:
      case COMPONENT::P2:
      case COMPONENT::P3:
        return -degree * (pow(x[0], degree - 1) + pow(x[1], degree - 1));
      case COMPONENT::V_X:
        return -(numL + 1) * degree * pow(x[0], degree - 1);
      case COMPONENT::V_Y:
        return -(numL + 1) * degree * pow(x[1], degree - 1);
      default:
        return 0.0;
    }
  }

  double Dut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P1:
      case COMPONENT::P2:
      case COMPONENT::P3:
        return ut(t, x, c, comp);
      default:
        return 0.0;
    }
  }

};

class ZeroProblem : public AcousticProblem {
public:
  const bool dirichlet;

  explicit ZeroProblem(bool dirichlet = false) : AcousticProblem(2, 0),
                                                 dirichlet(dirichlet) {}

  std::string Name() const override { return "ZeroProblem"; }

  int BndID(const Point &z) const override { return dirichlet ? 1 : 2; }

  bool hasInitialData() const override { return false; }

  bool HasRHS() const override { return false; }

  bool HasExactSolution() const override { return true; }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }
};

class MaterialJump : public AcousticProblem {
  VectorField normal = VectorField(1.0, 0.0, 0.0);
  double threshold = 2.0 / 3.0;
  bool constantMaterial = false;
public:
  MaterialJump() : AcousticProblem(2, 0) {
    Config::Get("normal_x", normal[0]);
    Config::Get("normal_y", normal[1]);
    Config::Get("threshold", threshold);
    Config::Get("constantMaterial", constantMaterial);
    normal /= norm(normal);
  }

  std::string Name() const override {
    std::stringstream ss;
    ss << "MaterialJump with threshold=" << threshold
       << " and n=[" << normal[0] << ", " << normal[1] << "]";
    return ss.str();
  }

  bool hasInitialData() const override { return false; }

  bool HasExactSolution() const override { return true; }

  int BndID(const Point &z) const override { return 2; }

  double Rho(const Cell &c, const Point &x) const override {
    double xw = x * normal;
    if (xw <= threshold || constantMaterial) {
      return 1.0;
    } else {
      return 2.0;
    }
  }

  double F(double t, const Cell &c, const Point &x, COMPONENT comp) const override {
    return 0.0;
  }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    if (i != 0) { THROW("Kappa_i called with i != 0")}
    return 1.0 / Rho(c, x);
  }

  double a0(double xw) const {
    if (0 <= xw && xw <= 1.0 / 3.0) {
      double s = sin(3.0 * Pi * xw);
      return pow(s, 2);
    }
    return 0;
  }

  double d_a0(double xw) const {
    if (0 <= xw && xw <= 1.0 / 3.0) {
      double s = sin(3.0 * Pi * xw);
      double d_s = 3.0 * Pi * cos(3.0 * Pi * xw);
      return d_s * pow(s, 1);
    }
    return 0;
  }

  double u0(const Point &x, COMPONENT comp) const {
    double xw = x * normal;
    double a = a0(xw);
    switch (comp) {
      case COMPONENT::P0:
        return a;
      case COMPONENT::V_X:
        return -a * normal[0];
      case COMPONENT::V_Y:
        return -a * normal[1];
    }
    Exit("u0 not defined for component " + to_string(comp));
    return 0.0;
  }

  double ut(const Point &xt, COMPONENT comp) const {
    double t = xt.t();
    if (t == 0) return u0(xt, comp);
    Point x = Point(xt[0], xt[1]);
    double xw = x * normal;
    double a = a0(xw);
    VectorField arg;
    if (xw <= threshold || constantMaterial) {
      arg = x - t * normal;
    } else {
      arg = 2 * x - (t + threshold) * normal;
    }
    return u0(Point(arg[0], arg[1]), comp);
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    double xw = x * normal;
    if (xw <= threshold || constantMaterial) {
      double d_a = d_a0(xw - t);
      switch (comp) {
        case COMPONENT::P0:
          return -d_a;
        case COMPONENT::V_X:
          return d_a * normal[0];
        case COMPONENT::V_Y:
          return d_a * normal[1];
      }
    } else {
      double d_a = d_a0(2 * xw - t - threshold);
      switch (comp) {
        case COMPONENT::P0:
          return -2 * d_a;
        case COMPONENT::V_X:
          return 2 * d_a * normal[0];
        case COMPONENT::V_Y:
          return 2 * d_a * normal[1];
      }
    }
    return 0;
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    double xw = x * normal;
    if (xw <= threshold || constantMaterial) {
      double d_a = d_a0(xw - t);
      switch (comp) {
        case COMPONENT::P0:
          return d_a;
        case COMPONENT::V_X:
          return -d_a * normal[0];
        case COMPONENT::V_Y:
          return -d_a * normal[1];
      }
    } else {
      double d_a = d_a0(2 * xw - t - threshold);
      switch (comp) {
        case COMPONENT::P0:
          return 2 * d_a;
        case COMPONENT::V_X:
          return -2 * d_a * normal[0];
        case COMPONENT::V_Y:
          return -2 * d_a * normal[1];
      }
    }
    return 0;
  }
};


class CircularInclusion : public AcousticProblem {
  Point initialPoint;
  double initialWidth;
  Point inclusionPoint;
  double inclusionWidth;
  double inclusionKappa;
  double inclusionRho;

public:
  CircularInclusion() : AcousticProblem(2, 0),
                        initialPoint(0.25, 0.25, 0.0),
                        initialWidth(0.1),
                        inclusionPoint(0.75, 0.75, 0.0),
                        inclusionWidth(0.3),
                        inclusionKappa(1.0),
                        inclusionRho(3.0) {
    Config::Get("initialPoint", initialPoint);
    Config::Get("inclusionPoint", inclusionPoint);
  }

  double Rho(const Cell &c, const Point &x) const override {
    return dist(x, inclusionPoint) < inclusionWidth ? inclusionRho : 1;
  }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    if (i != 0) { THROW("Kappa_i called with i != 0")}
    return dist(x, inclusionPoint) < inclusionWidth ? inclusionKappa : 1;
  }

  double A0(double s) const {
    if (s < initialWidth) return pow(cos(s * 0.5 * M_PI / initialWidth), 6.0);
    return 0.0;
  }

  double A(const Point &x) const {
    double r = abs(initialPoint[0] - x[0]);
    return A0(r);
  }

  double u0(const Point &x, COMPONENT comp) const {
    const double tmp = A(x);
    switch (comp) {
      case COMPONENT::P0:
        return tmp;
      case COMPONENT::V_X:
        return -tmp;
      default:
        return 0.0;
    }
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (x.t() > 0) return 0.0;
    return u0(x, comp);
  }

//  double ut(double t, const Point &x, const Cell &tc, COMPONENT comp) const override {
//    return ut(t, x, c, comp);
//  }

  std::string Name() const override {
    return "ViscoAcoustic_CircularInclusion";
  }

  bool HasExactSolution() const override { return false; }

  bool hasInitialData() const override { return true; }
};

class Maxwell : public AcousticProblem {
private:
  Point Psource = {0.5, 0.5, 0.05};
  double width = 0.1;
  double pulse_duration = 0.2;
  double source_factor = 1.0;
  double pml = 25;
public:
  Maxwell() : AcousticProblem(2, 1) {
    tau_p = 1.0;
  }

  double PML(const Point &z) const {
    double pml_size = 0.25;
    // Distance of quad with sides xl, xr, yt, yb
    // (max(x-xr,0)^2+max(y-yt,0)^2+max(xl-x,0)^2+max(yb-y,0)^2)^0.5
    double xl = 0 + pml_size;
    double xr = 3 - pml_size;
    double yb = -1 + pml_size;
    double yt = 2 - pml_size;
    double x = z[0];
    double y = z[1];
    double xx = pow(max(x - xr, 0.0), 2) + pow(max(xl - x, 0.0), 2);
    double yy = pow(max(y - yt, 0.0), 2) + pow(max(yb - y, 0.0), 2);
    double s = pow(xx + yy, 0.5);
    if (s == 0) {
      return 1.0;
    }
    return exp(pml * s);
  }

  double Rho(const Cell &c, const Point &x) const override { return rho / PML(x); }

  double Kappa_i(const Cell &c, const Point &x, int var) const override {
    double kappa = this->kappa * PML(x);
    if (var == -1) return kappa * (1 + numL * tau_p);
    if (var == 1) return kappa;
    else return kappa * tau_p;
  }

  double Tau_i(const Cell &c, const Point &x, int var) const override {
    Point mid{1, 0.5};
    if (norm(x - mid) < 0.25) {
      return 1.0 / 100000.0;
    }

    if (numL != 1) Exit("MaxwellVA only for numL=1");
    return 1;
  }

  virtual double TimeFunction(double t_diff) const {
    double pi2 = M_PI * M_PI;
    double f2 = 1.0 / (pulse_duration * pulse_duration);
    double a = pi2 * f2 * t_diff * t_diff;
    return (1 - 2 * a) * exp(-a);
  }

  virtual double SpaceFunction(const double &dist) const {
    if (dist >= width) return 0.0;
    return pow(cos(M_PI / 2 * dist / width), 6);
  }

  bool HasRHS() const override { return true; }

  double F(double t, const Cell &c, const Point &x,
           COMPONENT comp) const override {
    if (comp == COMPONENT::P0)
      return source_factor * TimeFunction(x.t() - Psource.t())
             * SpaceFunction(dist(Psource, x));
    return 0.0;
  }

  virtual double ut(double t, const Point &x, int) const { return 0.0; }

  std::string Name() const override { return "Maxwell"; };
};


class HighDensityProblem : public AcousticProblem {
  Point shotLocation{0.25, 0.25, 0.0};
  double factor = 10.0;
public:
  HighDensityProblem() : AcousticProblem(2, 0) {}

  std::string Name() const override { return "HighDensityProblem"; }

  double Rho(const Cell &c, const Point &x) const override {

    return 1 + 50 * pow(sin(10 * x[0]), 2) + 30 * pow(cos(10 * x[1]), 2);

    /*if (x[0] < 0.05 || x[0] > 0.95 || x[1] > 0.95 || x[1] < 0.05){
       return 1e20;
     }
     double r = dist(x, Point(0.5, 0.5));
     if (r > 0.25) {
       return 1;
     }
     double r2 = r * r;
     return 1 + (1.0 / 16 - r2) * 5000 * std::pow(3, -r2);*/
  }

  bool HasRHS() const override { return true; }

  double F(double t, const Cell &c, const Point &x,
           COMPONENT comp) const override {
    if (comp != COMPONENT::P0) return 0.0;
    return factor * Ricker(x.t() - shotLocation.t(), 0.1) * CosHat(dist(shotLocation, x), 0.1);
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override { return 0.0; }
};

class QuadraticST : public AcousticProblem {
public:
  QuadraticST() : AcousticProblem(2, 0) {}

  std::string Name() const override { return "QuadraticST"; }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return x[0] * x[1];
      case COMPONENT::V_X:
        return x[0] * x.t();
      case COMPONENT::V_Y:
        return x[1] * x[1];
      default:
        return 0.0;
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return 0 * 2 * x.t();
      case COMPONENT::V_X:
        return x[0];
      case COMPONENT::V_Y:
        return 0.0;
      default:
        return 0.0;
    }
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return -(x.t() + 2 * x[1]);
      case COMPONENT::V_X:
        return -x[1];
      case COMPONENT::V_Y:
        return -x[0];
      default:
        return 0.0;
    }
  }
};

class QuadraticDirichlet : public QuadraticST {
public:
  int BndID(const Point &z) const override { return 1; }
};

class QuadraticPressure : public AcousticProblem {
public:
  QuadraticPressure() : AcousticProblem(2, 0) {}

  std::string Name() const override {
    return "AcousticQuadraticPressure";
  }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return x[0] * x[0] + x[1] * x[1] + x.t() * x.t();
      case COMPONENT::V_X:
        return 0.0;
      case COMPONENT::V_Y:
        return 0.0;
      default:
        return 0.0;
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return 2.0 * x.t();
      case COMPONENT::V_X:
        return 0.0;
      case COMPONENT::V_Y:
        return 0.0;
      default:
        return 0.0;
    }
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return 0;
      case COMPONENT::V_X:
        return -2.0 * x[0];
      case COMPONENT::V_Y:
        return -2.0 * x[1];
      default:
        return 0.0;
    }
  }
};

class QuadraticV_X : public AcousticProblem {
public:
  QuadraticV_X() : AcousticProblem(2, 0) {}

  std::string Name() const override {
    return "AcousticQuadraticV_X";
  }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return 0;
      case COMPONENT::V_X:
        return x[0] * x[0];
      case COMPONENT::V_Y:
        return 0;
      default:
        return 0.0;
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return 0.0;
      case COMPONENT::V_X:
        return 0.0;
      case COMPONENT::V_Y:
        return 0.0;
      default:
        return 0.0;
    }
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return -2.0 * (x[0]);
      case COMPONENT::V_X:
        return 0.0;
      case COMPONENT::V_Y:
        return 0.0;
      default:
        return 0.0;
    }
  }
};

class QuadraticV_Y : public AcousticProblem {
public:
  QuadraticV_Y() : AcousticProblem(2, 0) {}

  std::string Name() const override {
    return "AcousticQuadraticV_Y";
  }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return 0;
      case COMPONENT::V_X:
        return x[1] * x[1];
      case COMPONENT::V_Y:
        return 0;
      default:
        return 0.0;
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return 0.0;
      case COMPONENT::V_X:
        return 0.0;
      case COMPONENT::V_Y:
        return 0.0;
      default:
        return 0.0;
    }
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return 0.0;
      case COMPONENT::V_X:
        return 0.0;
      case COMPONENT::V_Y:
        return 0.0;
      default:
        return 0.0;
    }
  }
};

class RiemannJump : public AcousticProblem {
  double left_value{};
  double right_value{};
  double jumpX;
  double kappa;
public:
  RiemannJump()
      : AcousticProblem(2, 0),
        left_value{1.0},
        right_value{0.0},
        jumpX{0.25},
        kappa{4} {
    Config::Get("jumpX", jumpX);
    Config::Get("kappa", kappa);
  }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    if (i != 0) { THROW("Kappa_i called with i != 0")}
    return x[0] > jumpX ? kappa : 1;
  }

  bool HasExactSolution() const override {
    return true;
  }

  int BndID(const Point &z) const override { return 2; }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    double xx = x[0];
    double zr = sqrt(Kappa(c, x) * Rho(c, x));
    if (xx > jumpX) {
      double dx = xx - jumpX;
      double dt = t - 1 * jumpX;
      if (2 * t > -dx + 1) {
        switch (comp) {
          case COMPONENT::V_X:
            return -1.0 / 3;
            break;
          case COMPONENT::V_Y:
            return 0;
            break;
          case COMPONENT::P0:
            return 4.0 / 3;
            break;
        }
      }
      if (2 * dt > dx) {
        switch (comp) {
          case COMPONENT::V_X:
            return -1.0 / 3;
            break;
          case COMPONENT::V_Y:
            return 0;
            break;
          case COMPONENT::P0:
            return 2.0 / 3;
            break;
        }
      }
    } else {
      if (t > -xx + 2 * jumpX) {
        switch (comp) {
          case COMPONENT::V_X:
            return -1.0 / 3;
            break;
          case COMPONENT::V_Y:
            return 0;
            break;
          case COMPONENT::P0:
            return 2.0 / 3;
            break;
        }
      }

      if (t * 1 > abs(xx)) {
        switch (comp) {
          case COMPONENT::V_X:
            return -1.0 / 2;
            break;
          case COMPONENT::V_Y:
            return 0;
            break;
          case COMPONENT::P0:
            return 1.0 / 2;
            break;
        }
      }
    }
    //if (x.t() > 0 || i != COMPONENT::P0) return 0.0;
    switch (comp) {
      case COMPONENT::V_X:
        return 0;
        break;
      case COMPONENT::V_Y:
        return 0;
        break;
      case COMPONENT::P0:
        if (c()[0] < 0) {
          return left_value;
        } else {
          return right_value;
        }
        break;
      default: Exit("Problem component not defined");
    }
  }

  bool hasInitialData() const override { return true; }

  std::string Name() const override {
    return "RiemannJumpProblem";
  }
};

class DoubleRiemann : public AcousticProblem {
  VectorField normal = VectorField(1.0, 0.0);
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
  DoubleRiemann()
      : AcousticProblem(2, 0) {
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
    //mout << DOUT(beta_R) << DOUT(beta_L) << endl;
  }

  bool HasExactSolution() const override {
    return true;
  }

  int BndID(const Point &z) const override { return 1; }

  double ut(double t, const Point &z, const Cell &tc, COMPONENT comp) const override {
    double x = z * normal;
    if (x <= -c_L * t) {
      switch (comp) {
        case COMPONENT::V_X:
          //mout << "vx <-ct" << v_L[0] << endl;
          return v_L[0];
        case COMPONENT::V_Y:
          //mout << "vy <-ct" << v_L[1]<< endl;
          return v_L[1];
        case COMPONENT::P0:
          return p_L;
        default: Exit("i > 2");
      }
    } else if (x < 0) {
      switch (comp) {
        case COMPONENT::V_X:
          //mout << "vx <0" << v_L[0] + beta_L * normal[0]<< endl;
          return v_L[0] + beta_L * normal[0];
        case COMPONENT::V_Y:
          //mout << "vy <0" <<  v_L[1] + beta_L * normal[1]<< endl;
          return v_L[1] + beta_L * normal[1];
        case COMPONENT::P0:
          return p_L + beta_L * Z_L;
        default: Exit("i > 2");
      }
    } else if (x < c_R * t) {
      switch (comp) {
        case COMPONENT::V_X:
          //mout << "vx <ct" << v_R[0] + beta_R * normal[0]<< endl;
          return v_R[0] + beta_R * normal[0];
        case COMPONENT::V_Y:
          //mout << "vy <ct" << v_R[1] + beta_R * normal[1]<< endl;
          return v_R[1] + beta_R * normal[1];
        case COMPONENT::P0:
          return p_R - beta_R * Z_R;
        default: Exit("i > 2");
      }
    } else {
      switch (comp) {
        case COMPONENT::V_X:
          //mout << "vx <else" << v_R[0]<< endl;
          return v_R[0];
        case COMPONENT::V_Y:
          //mout << "vy <else" <<  v_R[1]<< endl;
          return v_R[1];
        case COMPONENT::P0:
          return p_R;
        default: Exit("i > 2");
      }
    }
  }

  double Mdtut(double t, const Point &z, const Cell &c, COMPONENT comp) const override {
    return 0.0;
    double x = z * normal;
    if (x == -c_L * t) return infty;
    if (x == c_R * t) return infty;
    return 0.0;
  }

  double Aut(double t, const Point &z, const Cell &c, COMPONENT comp) const override {
    return 0.0;
    double x = z * normal;
    if (x == -c_L * t) return infty;
    if (x == c_R * t) return infty;
    return 0.0;
  }

  double F(double t, const Cell &c, const Point &x,
           COMPONENT comp) const override { return 0.0; }

//  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
//    Exit("call other ut with cell")
//  }

  bool hasInitialData() const override { return true; }

  std::string Name() const override {
    std::stringstream ss;
    ss << "DoubleRiemannProblem with n=[" << normal[0] << ", " << normal[1] << "]";
    return ss.str();
  }
};

class ConstantST : public AcousticProblem {
public:
  ConstantST() : AcousticProblem(2, 0) {}

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::V_X:
        return 3.0;
      case COMPONENT::V_Y:
        return 2.0;
      case COMPONENT::P0:
        return 1.0;
      default:
        return 0.0;
    }
  }

  bool HasExactSolution() const override { return true; }

  bool hasInitialData() const override { return true; }

  std::string Name() const override {
    return "ConstantProblem, Solution: (3,0,0,0) ";
  }
};

class LinearST : public AcousticProblem {
public:
  LinearST() : AcousticProblem(2, 0) {}

  std::string Name() const override { return "LinearST"; }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return x[1] - x[0];
      case COMPONENT::V_X:
        return x.t();
      case COMPONENT::V_Y:
        return x[1];
      default:
        return 0.0;
    }
  }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::V_X:
        return 1.0;
      default:
        return 0.0;
    }
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    switch (comp) {
      case COMPONENT::P0:
        return 0;
      case COMPONENT::V_X:
        return 0;
      case COMPONENT::V_Y:
        return 0;
      default:
        return 0.0;
    }
  }

  bool hasInitialData() const override { return true; }

  bool HasRHS() const override { return true; }

  bool HasExactSolution() const override { return true; }
};

class SinCosST : public AcousticProblem {
public:

  SinCosST() : AcousticProblem(2, 0) {}

  bool HasExactSolution() const override { return true; }

  bool HasRHS() const override { return false; }

  bool hasInitialData() const override { return true; }

  std::string Name() const override { return "SinCos"; }

  double F(double t, const Cell &c, const Point &x,
           COMPONENT comp) const override {
    return 0.0;
  }

  double Mdtut(double t, const Point &xx, const Cell &c, COMPONENT comp) const override {
    double x = xx[0];
    double y = xx[1];
    double s13 = sqrt(13);
    double s13i = 1.0 / s13;
    const double pi = M_PI;
    //double p = -           sin(M_PI * x) * sin(1.5 * M_PI * y) * sin(s13 * M_PI * t / 2);
    //double vx = 2 * s13i * cos(M_PI * x) * sin(1.5 * M_PI * y) * cos(s13 * M_PI * t / 2);
    //double vy = 3 * s13i * sin(M_PI * x) * cos(1.5 * M_PI * y) * cos(s13 * M_PI * t / 2);

    //double p_dt =  -s13 * M_PI/2 * sin(M_PI * x) * sin(1.5 * M_PI * y) * cos(s13 * M_PI * t / 2);
    //double vx_dt = -M_PI *         cos(M_PI * x) * sin(1.5 * M_PI * y) * sin(s13 * M_PI * t / 2);
    //double vy_dt = -1.5 * M_PI *   sin(M_PI * x) * cos(1.5 * M_PI * y) * sin(s13 * M_PI * t / 2);
    double p_dt =
        -sqrt(13) * pi * sin(pi * x) * sin(3 * pi * y / 2) * cos(sqrt(13) * pi * t / 2) / 2;
    double vx_dt = -pi * cos(pi * x) * sin(3 * pi * y / 2) * sin(sqrt(13) * pi * t / 2);
    double vy_dt = -3 * pi * sin(pi * x) * cos(3 * pi * y / 2) * sin(sqrt(13) * pi * t / 2) / 2;

    switch (comp) {
      case COMPONENT::P0:
        return p_dt;
      case COMPONENT::V_X:
        return vx_dt;
      case COMPONENT::V_Y:
        return vy_dt;
      default: Exit("No.");
    }
  }

  double Aut(double t, const Point &xx, const Cell &c, COMPONENT comp) const override {
    double x = xx[0];
    double y = xx[1];
    double s13 = sqrt(13);
    double s13i = 1.0 / s13;
    const double pi = M_PI;
    //double p = -                    sin(M_PI * x)    * sin(1.5 * M_PI * y)    * sin(s13 * M_PI * t / 2);
    double p_dx = -cos(pi * x) * sin(3 * pi * y / 2) * sin(sqrt(13) * pi * t / 2);
    double p_dy = -1.5 * sin(pi * x) * cos(3 * pi * y / 2) * sin(sqrt(13) * pi * t / 2);
    //double vx = 2 * s13i *          cos(M_PI * x)    * sin(1.5 * M_PI * y)    * cos(s13 * M_PI * t / 2);
    double vx_dx = -2 * pi * s13i * sin(pi * x) * sin(3 * pi * y / 2) * cos(sqrt(13) * pi * t / 2);
    //double vy = 3 * s13i *          sin(M_PI * x)    * cos(1.5 * M_PI * y)    * cos(s13 * M_PI * t / 2);
    double vy_dy =
        -4.5 * M_PI * s13i * sin(pi * x) * sin(3 * pi * y / 2) * cos(sqrt(13) * pi * t / 2);

    switch (comp) {
      case COMPONENT::P0:
        return -(vx_dx + vy_dy);
      case COMPONENT::V_X:
        return -p_dx;
      case COMPONENT::V_Y:
        return -p_dy;
      default: Exit("No.");
    }
  }

  double ut(double t, const Point &xx, const Cell &c, COMPONENT comp) const override {
    double x = xx[0];
    double y = xx[1];
    double s13 = sqrt(13);
    double s13i = 1.0 / s13;

    double p = -sin(M_PI * x) * sin(1.5 * M_PI * y) * sin(s13 * M_PI * t / 2);
    double vx = 2 * s13i * cos(M_PI * x) * sin(1.5 * M_PI * y) * cos(s13 * M_PI * t / 2);
    double vy = 3 * s13i * sin(M_PI * x) * cos(1.5 * M_PI * y) * cos(s13 * M_PI * t / 2);
    switch (comp) {
      case COMPONENT::P0:
        return p;
      case COMPONENT::V_X:
        return vx;
      case COMPONENT::V_Y:
        return vy;
      default: Exit("No.");
    }
  }
};

class Benchmark1DST : public AcousticProblem {
  double Rho1 = 0.5;
  double Rho2 = 2.0;
  double startRho1 = 0.25;
  double startRho2 = 0.5;
  double shift = 0.1;
  double width = 0.1;
  int linear_weight = 1;
public:

  Benchmark1DST() : AcousticProblem(1, 0) {
    Config::Get("width", width);
    Config::Get("shift", shift);
    Config::Get("linear_weight", linear_weight);
  }

  std::string Name() const override {
    return "Benchmark1DST";
  }

  double Rho(const Cell &c, const Point &x) const override {
    if (x[0] < startRho1) return 1.0;
    else if (x[0] < startRho2) return Rho1;
    else return Rho2;
  }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    if (i != 0) { THROW("Kappa_i called with i != 0")}
    return 1 / Rho(c, x);
  }

  double A0(double s) const {
    if (abs(s) < width) return pow(cos(s * 0.5 * M_PI / width), 6.0);
    return 0;
  }

  bool hasInitialData() const override {
    return true;
  }

  double A(double s) const { return A0(s - shift); }

  double DA0(double s) const {
    if (abs(s) < width)
      return 3.0 * M_PI / width * pow(cos(s * 0.5 * M_PI / width), 5.0)
             * sin(s * 0.5 * M_PI / width);
    else return 0.0;
  }

  double DA(double s) const { return DA0(s - shift); }

  double u0(const Point &x, COMPONENT comp) const {
    const double tmp = A(x[0]);
    switch (comp) {
      case COMPONENT::P0:
        return tmp;
      case COMPONENT::V_X:
        return -tmp;
      default:
        return 0.0;
    }
  }

  double ut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    if (t == 0)
      return u0(x, comp);
    if (x[0] < 0.25) {
      x[0] -= t;
      return u0(x, comp);
    } else if (x[0] < 0.5) {
      x[0] = Rho1 * (x[0] + startRho1) - t;
      return u0(x, comp);
    } else {
      x[0] = Rho1 + Rho2 * (x[0] + startRho2) - t;
      return u0(x, comp);
    }
  }

  bool HasExactSolution() const override {
    return true;
  }

  double Mdtut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    x[0] -= x.t();
    const double tmp = DA(x[0]);
    switch (comp) {
      case COMPONENT::P0:
        return -tmp;
      case COMPONENT::V_X:
        return tmp;
      default:
        return 0.0;
    }
  }

  double Aut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    x[0] -= x.t();
    const double tmp = DA(x[0]);
    switch (comp) {
      case COMPONENT::P0:
        return -tmp;
      case COMPONENT::V_X:
        return tmp;
      default:
        return 0.0;
    }
  }

  double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const override {
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch (linear_weight) {
      case 1: {
        Point ba = b - a;
        if (isVelocity(comp))
          return 0;
        else return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }
};

class Benchmark2DST : public AcousticProblem {
  double Rho1 = 0.5;
  double Rho2 = 2.0;
  double shift = -1;
  double width = 1.0;
  int linear_weight = 1;
public:
  Benchmark2DST() : AcousticProblem(2, 0) {
    Config::Get("width", width);
    Config::Get("shift", shift);
    Config::Get("linear_weight", linear_weight);
  }

  std::string Name() const override { return "Benchmark2DST"; }

  double Rho(const Cell &c, const Point &x) const override {
    if (x[0] < 0) return 1.0;
    else if (x[0] < 1) return Rho1;
    else return Rho2;
  }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    if (i != 0) { THROW("Kappa_i called with i != 0")}
    return 1 / Rho(c, x);
  }

  double A0(double s) const {
    if (abs(s) < width) return pow(cos(s * 0.5 * M_PI / width), 6.0);
    return 0;
  }

  bool hasInitialData() const override { return true; }

  double A(double s) const { return A0(s - shift); }

  double DA0(double s) const {
    if (abs(s) < width)
      return 3.0 * M_PI / width * pow(cos(s * 0.5 * M_PI / width), 5.0)
             * sin(s * 0.5 * M_PI / width);
    else return 0.0;
  }

  double DA(double s) const { return DA0(s - shift); }

  double u0(const Point &x, COMPONENT comp) const {
    const double tmp = A(x[0]);
    switch (comp) {
      case COMPONENT::P0:
        return tmp;
      case COMPONENT::V_X:
        return -tmp;
      default:
        return 0.0;
    }
  }

  double ut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    if (x.t() == 0) return u0(x, comp);

    if (x[0] < 0) {
      x[0] -= t;
      return u0(x, comp);
    } else if (x[0] < 1) {
      x[0] = Rho1 * x[0] - t;
      return u0(x, comp);
    } else {
      x[0] = Rho1 + Rho2 * (x[0] - 1) - t;
      return u0(x, comp);
    }
  }

  bool HasExactSolution() const override { return true; }

  double Mdtut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    x[0] -= x.t();
    const double tmp = DA(x[0]);
    switch (comp) {
      case COMPONENT::P0:
        return -tmp;
      case COMPONENT::V_X:
        return tmp;
      default:
        return 0.0;
    }
  }

  double Aut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    x[0] -= x.t();
    const double tmp = DA(x[0]);
    switch (comp) {
      case COMPONENT::P0:
        return -tmp;
      case COMPONENT::V_X:
        return tmp;
      default:
        return 0.0;
    }
  }

  double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const override {
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch (linear_weight) {
      case 1: {
        Point ba = b - a;
        if (isVelocity(comp)) return 0;
        else return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }
};

class SphericalWaveST : public AcousticProblem {
  Point Pmid;
  int linear_weight = 1;
  double width;
public:
  SphericalWaveST()
      : AcousticProblem(2, 0), Pmid(0, 0, 0), width(0.25) {
    Config::Get("ProblemMid", Pmid);
    Config::Get("linear_weight", linear_weight);
  }

  virtual std::string Name() const override {
    return "SphericalWaveST";
  }

  double A0(double s) const {
    if (s < width) return pow(cos(s * 0.5 * M_PI / width), 6.0);
    else return 0.0;
  }

  double A(const Point &x) const {
    double r = dist(Pmid, x);
    return A0(r * r);
  }

  double u0(const Point &x, COMPONENT comp) const {
    const double tmp = A(x);
    if (comp == COMPONENT::V_X) return (x[0] - Pmid[0]) * tmp;
    if (comp == COMPONENT::V_Y) return (x[1] - Pmid[1]) * tmp;
    return 0.0;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (x.t() > 0) return 0.0;
    else return u0(x, comp);
  }

  double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const override {
    if (comp != COMPONENT::P0) return 0.0;
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch (linear_weight) {
      case 1: {
        Point ba = b - a;
        return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }

  virtual bool hasInitialData() const override { return true; }
};

class WaveProblemFromImages2d_RHS_VA : public AcousticProblem {
  double left, right, top, bottom;
  double min_rho, max_rho, min_kappa, max_kappa;

  int problemLevel;
  double mlLength; //macro-level length of cell
  short *RhoPixel;
  short *KappaPixel;
  int RhoWidth;
  int RhoHeight;
  int KappaWidth;
  int KappaHeight;

  double pow2L1;
  double pow2L;

  double pml;

public:
  std::string Name() const override { return "WaveProblemFromImages2d_RHS"; }

  WaveProblemFromImages2d_RHS_VA(int numL, const Point &minP, const Point &maxP)
      : AcousticProblem(2, numL),
        left(minP[0]), right(maxP[0]), top(maxP[1]), bottom(minP[1]),
        min_rho(1.0), max_rho(2.8),
        min_kappa(1.0), max_kappa(5.0),
        problemLevel(-1), mlLength(1.0), pml(0.0) {
    vout(1) << "WaveProblemFromImages2d_RHS: mesh bounding box:"
            << "(" << left << ", " << top << ")" << " - "
            << "(" << right << ", " << bottom << ")" << endl;

    string rho_image = "marmousi_rho";
    string kappa_image = "marmousi_kappa";
    Config::Get("ModelImageRho", rho_image);
    Config::Get("ModelImageKappa", kappa_image);
    if (PPM->Master(0)) {
      MyImg ImgRho, ImgKappa;
      std::stringstream rho_pp, mu_pp, kappa_pp;
      rho_pp << string(ProjectMppDir) + "/conf/models/" << rho_image << ".png";
      kappa_pp << string(ProjectMppDir) + "/conf/models/" << kappa_image << ".png";
      ImgRho.load_png(rho_pp.str().c_str());
      ImgKappa.load_png(kappa_pp.str().c_str());
      RhoWidth = ImgRho.width();
      RhoHeight = ImgRho.height();
      RhoPixel = new short[RhoWidth * RhoHeight];
      for (int i = 0; i < RhoWidth; ++i)
        for (int j = 0; j < RhoHeight; ++j)
          RhoPixel[i * RhoHeight + j] = ImgRho(i, j);
      KappaWidth = ImgKappa.width();
      KappaHeight = ImgKappa.height();
      KappaPixel = new short[KappaWidth * KappaHeight];
      for (int i = 0; i < KappaWidth; ++i)
        for (int j = 0; j < KappaHeight; ++j)
          KappaPixel[i * KappaHeight + j] = ImgKappa(i, j);
      PPM->BroadcastInt(RhoWidth);
      PPM->BroadcastInt(RhoHeight);
      PPM->Broadcast(RhoPixel, (RhoWidth * RhoHeight) * sizeof(short), 0);
      PPM->BroadcastInt(KappaWidth);
      PPM->BroadcastInt(KappaHeight);
      PPM->Broadcast(KappaPixel, (KappaWidth * KappaHeight) * sizeof(short), 0);
    } else {
      RhoWidth = PPM->BroadcastInt();
      RhoHeight = PPM->BroadcastInt();
      RhoPixel = new short[RhoWidth * RhoHeight];
      PPM->Broadcast(RhoPixel, (RhoWidth * RhoHeight) * sizeof(short), 0);
      KappaWidth = PPM->BroadcastInt();
      KappaHeight = PPM->BroadcastInt();
      KappaPixel = new short[KappaWidth * KappaHeight];
      PPM->Broadcast(KappaPixel, (KappaWidth * KappaHeight) * sizeof(short), 0);
    }

    vout(1) << "WaveProblemFromImages2d_RHS:: rho model read   '"
            << rho_image << "' (size: " << RhoWidth << "x"
            << RhoHeight << ")" << endl;
    vout(1) << "WaveProblemFromImages2d_RHS:: kappa model read '"
            << kappa_image << "' (size: " << KappaWidth << "x"
            << KappaHeight << ")" << endl;

    Config::Get("MinRho", min_rho);
    Config::Get("MaxRho", max_rho);
    Config::Get("MinKappa", min_kappa);
    Config::Get("MaxKappa", max_kappa);
    Config::Get("ProblemLevel", problemLevel);
    if (problemLevel != -1) {
      Config::Get("mlLength", mlLength);
      pow2L1 = pow(0.5, problemLevel + 1);
      pow2L = pow(2, problemLevel);
    }
    Config::Get("pml", pml);
  }

  double PML(const Point &z) const {
    double s = z[0] - 3.0;
    if (s < 0) return exp(pml * s);
    s = 7.0 - z[0];
    if (s >= 0) return 1;
    return exp(pml * s);
  }

  int toImgCoordHorizontal(double c, int height, int width) const {
    if (problemLevel < 0)
      return int(round((c - left) / (right - left) * width));

    c /= mlLength;
    double ganzzahl = int(c);
    double dezimal = c - ganzzahl;
    for (int k = 0; k <= pow2L; k++) {
      double fac = (1 + 2 * k) * pow2L1;
      if (abs(dezimal - fac) < pow2L1) {
        //mout << DOUT(ganzzahl)<< DOUT(dezimal)<< DOUT(pow2L)<< DOUT(pow2L1)
        //     << DOUT(k)<< DOUT((1 + 2 * k) * pow2L1)<< DOUT((right - left))
        //     << DOUT(width) << endl;
        return int(round(((ganzzahl + fac) * mlLength - left) / (right - left) * width));
      }
    }
    Exit("out of bounds");
  }

  int toImgCoordVertical(double c, int height, int width) const {
    if (problemLevel < 0)
      return int(round((c - bottom) / (top - bottom) * height));

    c /= mlLength;
    double ganzzahl = int(c);
    double dezimal = c - ganzzahl;
    for (int k = 0; k <= pow2L; k++) {
      double fac = (1 + 2 * k) * pow2L1;
      if (abs(dezimal - fac) < pow2L1) {
        return int(round(((ganzzahl + fac) * mlLength - bottom) / (top - bottom) * height));
      }
    }
    Exit("out of bounds");
  }

  double Rho(const Cell &c, const Point &P) const override {
    double x = P[0], y = P[1];
    int ix = toImgCoordHorizontal(x, RhoHeight, RhoWidth);
    int iy = toImgCoordVertical(y, RhoHeight, RhoWidth);
    return (min_rho + RhoPixel[ix * RhoHeight + iy] * (max_rho - min_rho) / 255) / PML(P);
  }

  double Kappa_i(const Cell &c, const Point &P, int i) const override {
    double x = P[0], y = P[1];
    int ix = toImgCoordHorizontal(x, KappaHeight, KappaWidth);
    int iy = toImgCoordVertical(y, KappaHeight, KappaWidth);
    double kappaTemp = KappaPixel[ix * KappaHeight + iy] * (max_kappa - min_kappa) / 255;
    const double Vp = PML(P) * (min_kappa + kappaTemp);
    const double kappa = Vp * Vp * Rho(c, P) / (1 + numL * tau_p);
    if (i == 0) return kappa;
    else return kappa * tau_p;
  }
};

class WaveProblemFromImages2d_RHS_VA2 : public AcousticProblem {
  ParameterImage rhoImage;
  ParameterImage vpImage;
  double pml = 0;
public:
  std::string Name() const override { return "WaveProblemFromImages2d_RHS"; }

  WaveProblemFromImages2d_RHS_VA2(int numL, const Point &minP, const Point &maxP)
      : AcousticProblem(2, numL),
        rhoImage("Rho",
                 string(ProjectMppDir) + "/conf/models/marmousi2-density.png",
                 std::make_pair(minP[0], maxP[0]),
                 std::make_pair(minP[1], maxP[1])),
        vpImage("Kappa",
                string(ProjectMppDir) + "/conf/models/marmousi2-vp.png",
                std::make_pair(minP[0], maxP[0]),
                std::make_pair(minP[1], maxP[1])) {
    vout(1) << "WaveProblemFromImages2d_RHS: mesh bounding box:"
            << "(" << minP[0] << ", " << maxP[0] << ")" << " - "
            << "(" << minP[1] << ", " << maxP[1] << ")" << endl;
    Config::Get("pml", pml);
  }

  double PML(const Point &z) const {
    double s = z[0] - 3.0;
    if (s < 0) return exp(pml * s);
    s = 7.0 - z[0];
    if (s >= 0) return 1;
    return exp(pml * s);
  }


  double Rho(const Cell &c, const Point &P) const override {
    return rhoImage.ValueAtPoint(P) / PML(P);
  }

  double Kappa_i(const Cell &c, const Point &P, int i) const override {
    double vpValue = this->vpImage.ValueAtPoint(P);
    const double Vp = PML(P) * vpValue;
    const double kappa = Vp * Vp * Rho(c, P) / (1 + numL * tau_p);
    if (i == 0) return kappa;
    else return kappa * tau_p;
  }
};


class WaveProblemFromImages2d_RHS_Ricker_VA : public WaveProblemFromImages2d_RHS_VA {
public:
  Point Pmid;
  double width;
  double pulse_duration;
  double source_factor;
  int linear_weight;

  explicit WaveProblemFromImages2d_RHS_Ricker_VA(int numL)
      : WaveProblemFromImages2d_RHS_VA(numL, Point(0, 0), Point(9, 3)),
        Pmid(-1.0, 1.0, 0.0, 0.4), width(0.2), pulse_duration(0.2),
        source_factor(100.0),
        linear_weight(1) {
    Config::Get("ProblemMid", Pmid);
    Config::Get("width", width);
    Config::Get("source_factor", source_factor);
    Config::Get("source_duration", pulse_duration);
    Config::Get("linear_weight", linear_weight);
  }

  std::string Name() const override {
    return "WaveProblemFromImages2d_RHS_Ricker_VA";
  }

  virtual double TimeFunction(double t_diff) const {
    double pi2 = M_PI * M_PI;
    double f2 = 1.0 / (pulse_duration * pulse_duration);
    double a = pi2 * f2 * t_diff * t_diff;
    return (1 - 2 * a) * exp(-a);
  }

  virtual double SpaceFunction(const double &dist) const {
    if (dist >= width) return 0.0;
    return pow(cos(M_PI / 2 * dist / width), 6);
  }

  bool HasRHS() const override { return true; }

  double F(double t, const Cell &c, const Point &x,
           COMPONENT comp) const override {
    if (comp == COMPONENT::P0) {
      double tf = TimeFunction(x.t() - Pmid.t());
      double sf = SpaceFunction(dist(Pmid, x));
      return source_factor * tf * sf;
    }
    return 0.0;
  }

  bool hasPointSource(Point &P) const override {
    P = Pmid;
    return false;
  }

  virtual double F_Pkt(const int &i, const double &t) const {
    if (i == 2) return source_factor * TimeFunction(t - Pmid.t());
    return 0.0;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override { return 0.0; }

  double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const override {
    if (comp != COMPONENT::P0) return 0;
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch (linear_weight) {
      case 1: {
        Point ba = b - a;
        return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }
};

class Marmousi2_VA : public WaveProblemFromImages2d_RHS_VA {
public:
  Point Pmid{};
  double width{};
  double pulse_duration{};
  double source_factor{};
  int linear_weight{};

  Marmousi2_VA()
      : WaveProblemFromImages2d_RHS_VA(3, Point(-4000.0, 0), Point(13000.0, 3500.0)),
        Pmid(0.0, 0.0, 0.0, 0.0), pulse_duration(0.1), source_factor(100.0),
        linear_weight(1), width(0.1) {
    Config::Get("ProblemMid", Pmid);
    Config::Get("source_factor", source_factor);
    Config::Get("source_duration", pulse_duration);
    Config::Get("linear_weight", linear_weight);

    tau = pulse_duration / (2 * M_PI);
  }

  std::string Name() const override {
    return "Marmousi2_VA";
  }

  virtual double TimeFunction(double t_diff) const {
    double pi2 = M_PI * M_PI;
    double f2 = 1.0 / (pulse_duration * pulse_duration);
    double a = pi2 * f2 * t_diff * t_diff;
    return (1 - 2 * a) * exp(-a);
  }

  bool hasPointSource(Point &P) const override {
    P = Pmid;
    return true;
  }

  virtual double F_Pkt(const int &i, const double &t) const {
    //        Exit("why? check!");
    if (i == 2) return source_factor * TimeFunction(t - Pmid.t());
    return 0.0;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override { return 0.0; }

  double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const override {
    if (comp != COMPONENT::P0) return 0.0;
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch (linear_weight) {
      case 1: {
        Point ba = b - a;
        return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }

  double Tau_i(const Cell &c, const Point &x, int i) const override {
    if (numL <= 1) return pulse_duration / (2 * M_PI);
    else if (numL == 3) {
      if (i == 1) return 1 / (0.151 * 2 * M_PI);
      if (i == 2) return 1 / (1.93 * 2 * M_PI);
      if (i == 3) return 1 / (18.9 * 2 * M_PI);
      Exit("ToDo");
    }
    Exit("ToDo");
  }
};

class Marmousi2_VA_rhs : public WaveProblemFromImages2d_RHS_VA {
public:
  Point Pmid;
  double width;
  double pulse_duration;
  double source_factor;
  int linear_weight;

  explicit Marmousi2_VA_rhs(int numL)
      : WaveProblemFromImages2d_RHS_VA(numL, Point(-4.0, 0), Point(13.0, 3.5)),
        Pmid(0.0, 0.0, 0.0, 0.0), pulse_duration(0.1), source_factor(100), width(0.1),
        linear_weight(1) {
    Config::Get("ProblemMid", Pmid);
    Config::Get("width", width);
    Config::Get("source_factor", source_factor);
    Config::Get("source_duration", pulse_duration);
    Config::Get("linear_weight", linear_weight);

    tau = pulse_duration / (2 * M_PI);
  }

  std::string Name() const override {
    return "Marmousi2_VA_rhs";
  }

  bool HasExactSolution() const override {
    return false;
  }

  double TimeFunction(double t_diff) const {
    double pi2 = M_PI * M_PI;
    double f2 = 1.0 / (pulse_duration * pulse_duration);
    double a = pi2 * f2 * t_diff * t_diff;
    return (1 - 2 * a) * exp(-a);
  }

  double SpaceFunction(const double &dist) const {
    if (dist >= width) return 0.0;
    return pow(cos(M_PI / 2 * dist / width), 6);
  }

  bool HasRHS() const override { return true; }

  double F(double t, const Cell &c, const Point &x, COMPONENT comp) const override {
    if (comp == COMPONENT::P0) {
      double tf = TimeFunction(x.t() - Pmid.t());
      double sf = SpaceFunction(dist(Pmid, x));
      return source_factor * tf * sf;
    }
    return 0.0;
  }

  bool hasPointSource(Point &P) const override {
    P = Pmid;
    return false;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override { return 0.0; }

  double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const override {
    if (comp != COMPONENT::P0) return 0.0;
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch (linear_weight) {
      case 1: {
        Point ba = b - a;
        return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }

  double Tau_i(const Cell &c, const Point &x, int i) const override {
    if (numL <= 1) return pulse_duration / (2 * M_PI);
    else if (numL == 3) {
      if (i == 1) return 1 / (0.151 * 2 * M_PI);
      if (i == 2) return 1 / (1.93 * 2 * M_PI);
      if (i == 3) return 1 / (18.9 * 2 * M_PI);
      Exit("ToDo");
    }
    Exit("ToDo");
  }

  int BndID(const Point &z) const override {
    if (z[1] == 0) return 2;
    return 3;
  }
};


class Marmousi2_VA_rhsIMAGE : public WaveProblemFromImages2d_RHS_VA2 {
public:
  Point Pmid;
  double width;
  double pulse_duration;
  double source_factor;
  int linear_weight;

  explicit Marmousi2_VA_rhsIMAGE(int numL)
      : WaveProblemFromImages2d_RHS_VA2(numL, Point(-4.0, 0), Point(13.0, 3.5)),
        Pmid(0.0, 0.0, 0.0, 0.0), pulse_duration(0.1), source_factor(100), width(0.1),
        linear_weight(1) {
    Config::Get("ProblemMid", Pmid);
    Config::Get("width", width);
    Config::Get("source_factor", source_factor);
    Config::Get("source_duration", pulse_duration);
    Config::Get("linear_weight", linear_weight);

    tau = pulse_duration / (2 * M_PI);
  }

  std::string Name() const override {
    return "Marmousi2_VA_rhs";
  }

  bool HasExactSolution() const override {
    return false;
  }

  double TimeFunction(double t_diff) const {
    double pi2 = M_PI * M_PI;
    double f2 = 1.0 / (pulse_duration * pulse_duration);
    double a = pi2 * f2 * t_diff * t_diff;
    return (1 - 2 * a) * exp(-a);
  }

  double SpaceFunction(const double &dist) const {
    if (dist >= width) return 0.0;
    return pow(cos(M_PI / 2 * dist / width), 6);
  }

  bool HasRHS() const override { return true; }

  double F(double t, const Cell &c, const Point &x,
           COMPONENT comp) const override {
    if (comp == COMPONENT::P0) {
      double tf = TimeFunction(x.t() - Pmid.t());
      double sf = SpaceFunction(dist(Pmid, x));
      return source_factor * tf * sf;
    }
    return 0.0;
  }

  bool hasPointSource(Point &P) const override {
    P = Pmid;
    return false;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override { return 0.0; }

  double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const override {
    if (comp != COMPONENT::P0) return 0.0;
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch (linear_weight) {
      case 1: {
        Point ba = b - a;
        return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }

  double Tau_i(const Cell &c, const Point &x, int i) const override {
    if (numL <= 1) return pulse_duration / (2 * M_PI);
    else if (numL == 3) {
      if (i == 1) return 1 / (0.151 * 2 * M_PI);
      if (i == 2) return 1 / (1.93 * 2 * M_PI);
      if (i == 3) return 1 / (18.9 * 2 * M_PI);
      Exit("ToDo");
    }
    Exit("ToDo");
  }

  int BndID(const Point &z) const override {
    if (z[1] == 0) return 2;
    return 3;
  }
};


class Marmousi2AcousticConstantSource : public WaveProblemFromImages2d_RHS_VA {
public:
  Point Pmid;
  double width;
  double pulse_duration;
  double source_factor;

  Marmousi2AcousticConstantSource()
      : WaveProblemFromImages2d_RHS_VA(0, Point(-4.0, 0), Point(13.0, 3.5)),
        Pmid(0.0, 0.0, 0.0, 0.0), pulse_duration(0.1), source_factor(100), width(0.1) {
    Config::Get("ProblemMid", Pmid);
    Config::Get("width", width);
    Config::Get("source_factor", source_factor);
    Config::Get("source_duration", pulse_duration);

    tau = pulse_duration / (2 * M_PI);
  }

  std::string Name() const override {
    return "Marmousi2ConstSource";
  }

  bool HasExactSolution() const override {
    return false;
  }

  double TimeFunction(double t_diff) const {
    if (t_diff > pulse_duration / 2) return 0;
    return 1;
  }

  double SpaceFunction(const double &dist) const {
    if (dist >= width) return 0.0;
    return 1;
  }

  bool HasRHS() const override { return true; }

  double F(double t, const Cell &c, const Point &x, COMPONENT comp) const override {
    if (comp == COMPONENT::P0) {
      double tf = TimeFunction(x.t() - Pmid.t());
      double sf = SpaceFunction(distInfty(Pmid, x));
      return source_factor * tf * sf;
    }
    return 0.0;
  }

  bool hasPointSource(Point &P) const override {
    P = Pmid;
    return false;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override { return 0.0; }


  double Tau_i(const Cell &c, const Point &x, int i) const override {
    if (numL <= 1) return pulse_duration / (2 * M_PI);
    else if (numL == 3) {
      if (i == 1) return 1 / (0.151 * 2 * M_PI);
      if (i == 2) return 1 / (1.93 * 2 * M_PI);
      if (i == 3) return 1 / (18.9 * 2 * M_PI);
      Exit("ToDo");
    }
    Exit("ToDo");
  }

  int BndID(const Point &z) const override {
    if (z[1] == 0) return 2;
    return 3;
  }
};

class Marmousi2AcousticConstantSourceIMAGE : public WaveProblemFromImages2d_RHS_VA2 {
public:
  Point Pmid;
  double width;
  double pulse_duration;
  double source_factor;

  Marmousi2AcousticConstantSourceIMAGE()
      : WaveProblemFromImages2d_RHS_VA2(0, Point(-4.0, 0), Point(13.0, 3.5)),
        Pmid(0.0, 0.0, 0.0, 0.0), pulse_duration(0.1), source_factor(100), width(0.1) {
    Config::Get("ProblemMid", Pmid);
    Config::Get("width", width);
    Config::Get("source_factor", source_factor);
    Config::Get("source_duration", pulse_duration);

    tau = pulse_duration / (2 * M_PI);
  }

  std::string Name() const override {
    return "Marmousi2ConstSource";
  }

  bool HasExactSolution() const override {
    return false;
  }

  double TimeFunction(double t_diff) const {
    if (t_diff > pulse_duration / 2) return 0;
    return 1;
  }

  double SpaceFunction(const double &dist) const {
    if (dist >= width) return 0.0;
    return 1;
  }

  bool HasRHS() const override { return true; }

  double F(double t, const Cell &c, const Point &x, COMPONENT comp) const override {
    if (comp == COMPONENT::P0) {
      double tf = TimeFunction(x.t() - Pmid.t());
      double sf = SpaceFunction(distInfty(Pmid, x));
      return source_factor * tf * sf;
    }
    return 0.0;
  }

  bool hasPointSource(Point &P) const override {
    P = Pmid;
    return false;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override { return 0.0; }


  double Tau_i(const Cell &c, const Point &x, int i) const override {
    if (numL <= 1) return pulse_duration / (2 * M_PI);
    else if (numL == 3) {
      if (i == 1) return 1 / (0.151 * 2 * M_PI);
      if (i == 2) return 1 / (1.93 * 2 * M_PI);
      if (i == 3) return 1 / (18.9 * 2 * M_PI);
      return 0.0;
    }
    Exit("ToDo");
  }

  int BndID(const Point &z) const override {
    if (z[1] == 0) return 2;
    return 3;
  }
};


class Acoustic_BenchmarkSimple_2D : public AcousticProblem {
  double shift;
  double width;
public:
  Acoustic_BenchmarkSimple_2D() : AcousticProblem(2, 0), width(1), shift(-1) {
    Config::Get("width", width);
    Config::Get("shift", shift);
  }

  std::string Name() const override {
    return "Acoustic_BenchmarkSimple_2D: u0 = exp(-1 / (width - s * s) + 1 / width)";
  }

  double Rho(const Cell &c, const Point &x) const override { return 1.0; }

  double Kappa_i(const Cell &c, const Point &x, int i) const override {
    if (i != 0) { THROW("Kappa_i called with i != 0")}
    return 1.0;
  }

  double A0(double s) const {
    if (s * s >= width) return 0;
    return exp(-1 / (width - s * s) + 1 / width);
  }

  double DA0(double s) const {
    if (s * s >= width) return 0;
    double y = width - s * s;
    return -2 * s * (1 / (y * y)) * exp(-1 / y) / exp(-1 / width);
  }

  double A(double s) const { return A0(s - shift); }

  double DA(double s) const { return DA0(s - shift); }

  double u0(const Point &x, COMPONENT comp) const {
    const double tmp = A(x[0]);
    switch (comp) {
      case COMPONENT::V_X:
        return tmp;
      case COMPONENT::V_Y:
        return tmp;
      default:
        return 0.0;
    }
  }

  double ut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    x[0] -= t;
    return u0(x, comp);
  }

  double Mdtut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    x[0] -= x.t();
    const double tmp = DA(x[0]);
    switch (comp) {
      case COMPONENT::V_X:
        return -tmp;
      case COMPONENT::V_Y:
        return -tmp;
      default:
        return 0.0;
    }
  }

  double Aut(double t, const Point &x_, const Cell &c, COMPONENT comp) const override {
    Point x = x_;
    x[0] -= x.t();
    const double tmp = DA(x[0]);
    switch (comp) {
      case COMPONENT::V_X:
        return tmp;
      case COMPONENT::V_Y:
        return tmp;
      default:
        return 0.0;
    }
  }

  bool HasExactSolution() const override { return true; }
};

class RiemannST : public AcousticProblem {
  VectorFieldT<double, 3> y_L{0, 0, 0.9};
  VectorFieldT<double, 3> y_R{0, 0, 1.1};
  VectorFieldT<double, 3> y_L_beta{};
  VectorFieldT<double, 3> y_R_beta{};

public:
  RiemannST() : AcousticProblem(2, 0) {
    VectorFieldT<double, 3> diff = y_R - y_L;
    double beta_R = diff[2] / 2;
    double beta_L = diff[2] / 2;
    y_L_beta = y_L + beta_L * VectorFieldT<double, 3>(-1, 0, 1);
    y_R_beta = y_R + beta_R * VectorFieldT<double, 3>(1, 0, -1);
  }

  bool hasInitialData() const override { return true; }

  std::string Name() const override {
    return "Acoustic_Riemann: x<0 then 0.9 else 1.1";
  }

  bool HasExactSolution() const override { return true; }

  double Mdtut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

  double Aut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

  double F(double t, const Cell &c, const Point &p, COMPONENT comp) const override {
    return 0.0;
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    int i = 0;
    for (COMPONENT co: GetComponents()) {
      if (co == comp) {
        break;
      }
      i++;
    }
    if (x.t() == 0) {
      if (c()[0] < 0) {
        return y_L[i];
      } else {
        return y_R[i];
      }
    }
    if (i >= 3) {
      return 0;
    }
    for (int tt = 0; tt < 5; tt++) {
      if ((0 < c()[0] && x[0] >= abs(x.t() - tt))) {
        return tt % 2 == 0 ? y_R[i] : y_L[i];
      }
      if ((0 > c()[0] && x[0] <= -abs(x.t() - tt))) {
        return tt % 2 == 0 ? y_L[i] : y_R[i];
      }
    }
    if (comp == COMPONENT::V_X) {
      if ((0 <= x.t() && x.t() < 1) || (2 <= x.t() && x.t() < 3)) {
        return 0.1;
      } else if ((1 <= x.t() && x.t() < 2) || (3 <= x.t() && x.t() < 4)) {
        return -0.1;
      }
    }
    if (c()[0] < 0) {
      return y_L_beta[i];
    } else {
      return y_R_beta[i];
    }
  }
};

class Acoustic_Plane_Wave_2D : public AcousticProblem {
  Point Pmid;
  int linear_weight;
public:
  Acoustic_Plane_Wave_2D() : AcousticProblem(2, 0), Pmid(0, 0, 0) {
    Config::Get("rho", rho);
    Config::Get("kappa", kappa);
    Config::Get("ProblemMid", Pmid);
    Config::Get("linear_weight", linear_weight);
  }

  std::string Name() const override {
    return "Acoustic_Plane_Wave_2D: exp(-r * r) * (1 - r * r);";
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    if (x.t() > 0) return 0.0;
    double r = 2 * abs(Pmid[1] - x[1]);
    switch (comp) {
      case COMPONENT::V_X:
        if (r < 1) return exp(-r * r) * (1 - r * r);
        else return 0.0;
      case COMPONENT::V_Y:
        if (r < 1) return -exp(-r * r) * (1 - r * r);
        else return 0.0;
      default:
        return 0.0;
    }
  }

  double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const override {
    switch (linear_weight) {
      case 1: {
        Point ba = b - a;
        if (comp == COMPONENT::V_X) return 0;
        if (comp == COMPONENT::V_Y) return 1.0 / (ba[0] * ba[1]);
        if (comp == COMPONENT::P0) return 1.0 / (ba[0] * ba[1]);
      }
      case 2: {
        Point PPmid = a + 0.5 * (b - a);
        double r = abs(PPmid[1] - x[1]) / abs(PPmid[1] - b[1]);
        switch (comp) {
          case COMPONENT::V_Y:
            if (r < 1) return exp(-r * r) * (1 - r * r);
            else return 0.0;
          case COMPONENT::P0:
            if (r < 1) return exp(-r * r) * (1 - r * r);
            else return 0.0;
          default:
            return 0.0;
        }
      }
      default: Exit("Not implemented");
    }
  }
};

class WaveProblemFromImages2d_RHS : public AcousticProblem {
protected:
  cimg_library::CImg<unsigned char> ImgRho;
  cimg_library::CImg<unsigned char> ImgKappa;
  double left, right, top, bottom;
  double min_rho, max_rho, min_kappa, max_kappa;

  Point Pmid;
  double width;
  double pulse_duration;
  double source_factor;
  int problemLevel;
public:
  std::string Name() const override { return "WaveProblemFromImages2d_RHS"; }

  WaveProblemFromImages2d_RHS(const Point &minP, const Point &maxP)
      : AcousticProblem(2, 0),
        left(minP[0]), right(maxP[0]), top(maxP[1]), bottom(minP[1]),
        min_rho(1.0), max_rho(3.0), min_kappa(1.0), max_kappa(3.0),
        Pmid(0.0, 0.0), width(0.005), pulse_duration(0.125), source_factor(100.0),
        problemLevel(-1) {
    vout(1) << "WaveProblemFromImages2d_RHS: mesh bounding box:"
            << "(" << left << ", " << top << ")" << " - "
            << "(" << right << ", " << bottom << ")" << endl;

    string rho_image = "marmousi_rho";
    string kappa_image = "marmousi_kappa";
    Config::Get("ModelImageRho", rho_image);
    Config::Get("ModelImageKappa", kappa_image);
    std::stringstream rho_pp, kappa_pp;
    rho_pp << string(ProjectMppDir) + "/models/" << rho_image << ".png";
    kappa_pp << string(ProjectMppDir) + "/models/" << kappa_image << ".png";
    ImgRho.load_png(rho_pp.str().c_str());
    ImgKappa.load_png(kappa_pp.str().c_str());

    vout(1) << "WaveProblemFromImages2d_RHS:: rho model read   '"
            << rho_image << "' (size: " << ImgRho.width() << "x"
            << ImgRho.height() << ")" << endl;
    vout(1) << "WaveProblemFromImages2d_RHS:: kappa model read '"
            << rho_image << "' (size: " << ImgKappa.width() << "x"
            << ImgKappa.height() << ")" << endl;

    Config::Get("MinRho", min_rho);
    Config::Get("MaxRho", max_rho);
    Config::Get("MinKappa", min_kappa);
    Config::Get("MaxKappa", max_kappa);

    Config::Get("ProblemMid", Pmid);
    Config::Get("width", width);
    Config::Get("source_factor", source_factor);
    Config::Get("source_duration", pulse_duration);

    Config::Get("ProblemLevel", problemLevel);
  }

  int toImgCoordHorizontal(double c, const cimg_library::CImg<unsigned char> &img) const {
    if (problemLevel < 0)
      return int(round((c - left) / (right - left) * img.width()));

    double ganzzahl = int(c);
    double dezimal = c - ganzzahl;
    double pow2L1 = pow(0.5, problemLevel + 1);
    for (int k = 0; k <= pow(2, problemLevel); k++)
      if (abs(dezimal - (1 + 2 * k) * pow2L1) <= pow2L1)
        return int(round((ganzzahl + (1 + 2 * k) * pow2L1 - left) / (right - left)
                         * img.width()));
    Exit("out of bounds");
  }

  int toImgCoordVertical(double c, const cimg_library::CImg<unsigned char> &img) const {
    if (problemLevel < 0)
      return int(round((c - bottom) / (top - bottom) * img.height()));

    double ganzzahl = int(c);
    double dezimal = c - ganzzahl;
    double pow2L1 = pow(0.5, problemLevel + 1);
    for (int k = 0; k <= pow(2, problemLevel); k++)
      if (abs(dezimal - (1 + 2 * k) * pow2L1) < pow2L1)
        return int(round(
            (ganzzahl + (1 + 2 * k) * pow2L1 - bottom) / (top - bottom)
            * img.height()));
    Exit("out of bounds");
  }

  double Rho(const Point &P) const {
    double x = P[0], y = P[1];
    int ix = toImgCoordHorizontal(x, ImgRho);
    int iy = toImgCoordVertical(y, ImgRho);
    return min_rho + ImgRho(ix, iy) * (max_rho - min_rho) / 255;
  }

  double Kappa(const Point &P) const {
    double x = P[0], y = P[1];
    int ix = toImgCoordHorizontal(x, ImgKappa);
    int iy = toImgCoordVertical(y, ImgKappa);
    return min_kappa + ImgKappa(ix, iy) * (max_kappa - min_kappa) / 255;
  }

  int BndID(const Point &z) const override {
    if (abs(z[1]) < 1e-16) return 2;
    else return 1;
  }
};


class PointSource1D : public AcousticProblem {
public:
  PointSource1D() : AcousticProblem(1, 0) {
  }

  bool HasRHS() const override {
    return true;
  }

  bool hasPointSource(Point &P) const override {
    P = {0.5, 0.0};
    return true;
  }

  std::string Name() const override {
    return "PointSource1D";
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

};

class PointSource2D : public AcousticProblem {
public:
  PointSource2D() : AcousticProblem(2, 0) {
  }

  bool HasRHS() const override {
    return true;
  }

  bool hasPointSource(Point &P) const override {
    P = {0.5, 0.5, 0.0};
    return true;
  }

  std::string Name() const override {
    return "PointSource2D";
  }

  double ut(double t, const Point &x, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

};

class CRCWithInitial : public AcousticProblem {
protected:

  double factor = 10.0;

  Point shotLocation;

public:
  explicit CRCWithInitial(Point shotLocation = {2.0, 0.0}) : AcousticProblem(2, 0) {
  }

  double F(double t, const Cell &c, const Point &x, COMPONENT comp) const override {
    if (comp == COMPONENT::P0) {
      return factor * Ricker(t - 0.5, 0.5) * GaussHat(dist(x, shotLocation), 0.3);
    }
    return 0.0;
  }

  bool HasRHS() const override { return true; }

  bool hasInitialData() const override { return true; }

  std::string Name() const override {
    return "CRCWithInitial";
  }
};

class FWIExampleProblem : public AcousticProblem {
protected:

  double factor = 10000;
  double delta = 0.0625;
  Point shotLocation{-0.125, 0, 0.05};

public:
  explicit FWIExampleProblem() : IProblem("FWIExampleGeometry"), AcousticProblem(2, 0) {
    Config::Get("ShotLocation", shotLocation);
    Config::Get("ShotFactor", factor);
    Config::Get("ShotDelta", delta);
  }

  explicit FWIExampleProblem(std::shared_ptr<Meshes> meshes) : IProblem(std::move(meshes)),
                                                               AcousticProblem(2, 0) {
    Config::Get("ShotLocation", shotLocation);
    Config::Get("ShotFactor", factor);
    Config::Get("ShotDelta", delta);
  }

  double F(double t, const Cell &c, const Point &x, COMPONENT comp) const override {
    if (comp == COMPONENT::P0) {
      double dist = normST(shotLocation - x);
      if (dist < delta) {
        return factor * exp(-6 * dist * dist / delta);
      }
    }
    return 0.0;
  }

  double Rho(const Cell &c, const Point &P) const override {
    return 1.0;
  }

  double Kappa_i(const Cell &c, const Point &P, int i) const override {
    if (i != 0) { THROW("Kappa_i called with i != 0")}
    double x2 = P[1];
    if (0.625 <= x2 && x2 <= 0.75) {
      return 1.0 + 0.55 * sin((x2 - 0.625) / 0.125);
    }
    return 1.0;
  }

  bool HasRHS() const override { return true; }

  bool hasInitialData() const override { return false; }

  double ut(double t, const Point &p, const Cell &c, COMPONENT comp) const override {
    return 0.0;
  }

  std::string Name() const override {
    return "FWIExampleProblem";
  }
};

AcousticProblem *CreateAcousticProblem(const string &name) {
  if (name == "AcousticBenchmarkC2")
    return new AcousticBenchmarkC2(AcousticBenchmarkC2::GetDampingFromConfig());
  if (name == "AcousticWaveDamping")
    return new AcousticWaveDamping(AcousticWaveDamping::GetDampingFromConfig());
  if (name == "AcousticWaveProcTest")
      return new AcousticWaveProcTest(AcousticWaveProcTest::GetDampingFromConfig());
  if (name == "GaussHatAndRicker2D")
    return new GaussHatAndRicker2D();
  if (name == "RiemannWave2D")
    return new RiemannWave2D();
  if (name == "Linear")
    return new Linear();
  if (name == "Quadratic")
    return new Quadratic();
  if (name == "CRC")
    return new CRC();


  if (name == "SinCos") return new SinCosST();

  if (name == "SphericalWave") return new SphericalWaveST();
  if (name == "MaterialJump") return new MaterialJump();
  if (name == "CircularInclusion") return new CircularInclusion();

  if (name == "HighDensity") return new HighDensityProblem();

  if (name == "RiemannST") return new RiemannST();
  if (name == "RiemannJump") return new RiemannJump();
  if (name == "DoubleRiemann") return new DoubleRiemann();

  if (name == "Zero") return new ZeroProblem();

  if (name == "Linear") return new LinearST();
  if (name == "Constant") return new ConstantST();
  if (name == "QuadraticST") return new QuadraticST();
  if (name == "QuadraticDirichlet") return new QuadraticDirichlet();
  if (name == "QuadraticV_X") return new QuadraticV_X();
  if (name == "QuadraticV_Y") return new QuadraticV_Y();
  if (name == "QuadraticPressure") return new QuadraticPressure();

  if (name == "Riemann1D") return new Riemann1DST();
  if (name == "Benchmark1D") return new Benchmark1DST();
  if (name == "PointSource1D") return new PointSource1D();

  if (name == "Benchmark2D") return new Benchmark2DST();
  if (name == "PointSource2D") return new PointSource2D();

  if (name == "Polynomial1DDegree0") return new PolynomialProblem1D(0);
  if (name == "Polynomial1DDegree1") return new PolynomialProblem1D(1);
  if (name == "Polynomial1DDegree2") return new PolynomialProblem1D(2);
  if (name == "Polynomial1DDegree3") return new PolynomialProblem1D(3);
  if (name == "Polynomial1DDegree4") return new PolynomialProblem1D(4);
  if (name == "Polynomial1DDegree5") return new PolynomialProblem1D(5);

  if (name == "Polynomial1DDegree0Dirichlet") return new PolynomialProblem1D(0, true);
  if (name == "Polynomial1DDegree1Dirichlet") return new PolynomialProblem1D(1, true);
  if (name == "Polynomial1DDegree2Dirichlet") return new PolynomialProblem1D(2, true);
  if (name == "Polynomial1DDegree3Dirichlet") return new PolynomialProblem1D(3, true);
  if (name == "Polynomial1DDegree4Dirichlet") return new PolynomialProblem1D(4, true);
  if (name == "Polynomial1DDegree5Dirichlet") return new PolynomialProblem1D(5, true);


  if (name == "Polynomial2DDegree0") return new PolynomialProblem2D(0);
  if (name == "Polynomial2DDegree1") return new PolynomialProblem2D(1);
  if (name == "Polynomial2DDegree2") return new PolynomialProblem2D(2);
  if (name == "Polynomial2DDegree3") return new PolynomialProblem2D(3);
  if (name == "Polynomial2DDegree4") return new PolynomialProblem2D(4);
  if (name == "Polynomial2DDegree5") return new PolynomialProblem2D(5);

  if (name == "Polynomial2DDegree0Dirichlet") return new PolynomialProblem2D(0, true);
  if (name == "Polynomial2DDegree1Dirichlet") return new PolynomialProblem2D(1, true);
  if (name == "Polynomial2DDegree2Dirichlet") return new PolynomialProblem2D(2, true);
  if (name == "Polynomial2DDegree3Dirichlet") return new PolynomialProblem2D(3, true);
  if (name == "Polynomial2DDegree4Dirichlet") return new PolynomialProblem2D(4, true);
  if (name == "Polynomial2DDegree5Dirichlet") return new PolynomialProblem2D(5, true);


  if (name == "Polynomial3DDegree0") return new PolynomialProblem3D(0);
  if (name == "Polynomial3DDegree1") return new PolynomialProblem3D(1);
  if (name == "Polynomial3DDegree2") return new PolynomialProblem3D(2);
  if (name == "Polynomial3DDegree3") return new PolynomialProblem3D(3);
  if (name == "Polynomial3DDegree4") return new PolynomialProblem3D(4);
  if (name == "Polynomial3DDegree5") return new PolynomialProblem3D(5);

  if (name == "Polynomial3DDegree0Dirichlet") return new PolynomialProblem3D(0, true);
  if (name == "Polynomial3DDegree1Dirichlet") return new PolynomialProblem3D(1, true);
  if (name == "Polynomial3DDegree2Dirichlet") return new PolynomialProblem3D(2, true);
  if (name == "Polynomial3DDegree3Dirichlet") return new PolynomialProblem3D(3, true);
  if (name == "Polynomial3DDegree4Dirichlet") return new PolynomialProblem3D(4, true);
  if (name == "Polynomial3DDegree5Dirichlet") return new PolynomialProblem3D(5, true);

  if (name == "SphericalProblem3D") return new SphericalProblem3DST();
  if (name == "TestProblem3D") return new TestProblem3DST();

  if (name == "Polynomial2DDegree0numL0") return new PolynomialProblemNumL(0, 0, false);
  if (name == "Polynomial2DDegree0numL1") return new PolynomialProblemNumL(0, 1, false);
  if (name == "Polynomial2DDegree0numL2") return new PolynomialProblemNumL(0, 2, false);
  if (name == "Polynomial2DDegree0numL3") return new PolynomialProblemNumL(0, 3, false);
  if (name == "Polynomial2DDegree1numL0") return new PolynomialProblemNumL(1, 0, false);
  if (name == "Polynomial2DDegree1numL1") return new PolynomialProblemNumL(1, 1, false);
  if (name == "Polynomial2DDegree1numL2") return new PolynomialProblemNumL(1, 2, false);
  if (name == "Polynomial2DDegree1numL3") return new PolynomialProblemNumL(1, 3, false);
  if (name == "Polynomial2DDegree2numL0") return new PolynomialProblemNumL(2, 0, false);
  if (name == "Polynomial2DDegree2numL1") return new PolynomialProblemNumL(2, 1, false);
  if (name == "Polynomial2DDegree2numL2") return new PolynomialProblemNumL(2, 2, false);
  if (name == "Polynomial2DDegree2numL3") return new PolynomialProblemNumL(2, 3, false);
  if (name == "Polynomial2DDegree3numL0") return new PolynomialProblemNumL(3, 0, false);
  if (name == "Polynomial2DDegree3numL1") return new PolynomialProblemNumL(3, 1, false);
  if (name == "Polynomial2DDegree3numL2") return new PolynomialProblemNumL(3, 2, false);
  if (name == "Polynomial2DDegree3numL3") return new PolynomialProblemNumL(3, 3, false);

  if (name == "Polynomial2DDegree0numL0Dirichlet") return new PolynomialProblemNumL(0, 0, true);
  if (name == "Polynomial2DDegree0numL1Dirichlet") return new PolynomialProblemNumL(0, 1, true);
  if (name == "Polynomial2DDegree0numL2Dirichlet") return new PolynomialProblemNumL(0, 2, true);
  if (name == "Polynomial2DDegree0numL3Dirichlet") return new PolynomialProblemNumL(0, 3, true);
  if (name == "Polynomial2DDegree1numL0Dirichlet") return new PolynomialProblemNumL(1, 0, true);
  if (name == "Polynomial2DDegree1numL1Dirichlet") return new PolynomialProblemNumL(1, 1, true);
  if (name == "Polynomial2DDegree1numL2Dirichlet") return new PolynomialProblemNumL(1, 2, true);
  if (name == "Polynomial2DDegree1numL3Dirichlet") return new PolynomialProblemNumL(1, 3, true);
  if (name == "Polynomial2DDegree2numL0Dirichlet") return new PolynomialProblemNumL(2, 0, true);
  if (name == "Polynomial2DDegree2numL1Dirichlet") return new PolynomialProblemNumL(2, 1, true);
  if (name == "Polynomial2DDegree2numL2Dirichlet") return new PolynomialProblemNumL(2, 2, true);
  if (name == "Polynomial2DDegree2numL3Dirichlet") return new PolynomialProblemNumL(2, 3, true);
  if (name == "Polynomial2DDegree3numL0Dirichlet") return new PolynomialProblemNumL(3, 0, true);
  if (name == "Polynomial2DDegree3numL1Dirichlet") return new PolynomialProblemNumL(3, 1, true);
  if (name == "Polynomial2DDegree3numL2Dirichlet") return new PolynomialProblemNumL(3, 2, true);
  if (name == "Polynomial2DDegree3numL3Dirichlet") return new PolynomialProblemNumL(3, 3, true);

  if (name == "CRCWithInitial") return new CRCWithInitial();
  if (name == "Maxwell") return new Maxwell();
  if (name == "Marmousi2") return new Marmousi2_VA();
  if (name == "Marmousi2_rhs") return new Marmousi2_VA_rhs(3);
  if (name == "Marmousi2_rhsIMAGE") return new Marmousi2_VA_rhsIMAGE(3);
  if (name == "WaveProblemFromImages2d_RHS_Ricker") {
    return new WaveProblemFromImages2d_RHS_Ricker_VA(3);
  }
  if (name == "Marmousi2AcousticConstantSource") {
    return new Marmousi2AcousticConstantSource();
  }
  if (name == "Marmousi2AcousticConstantSourceIMAGE") {
    return new Marmousi2AcousticConstantSourceIMAGE();
  }


  Exit(name + " Error: no Problem found!");
}

std::unique_ptr<AcousticProblem> CreateAcousticProblemUnique(const string &name) {
  return std::unique_ptr<AcousticProblem>(CreateAcousticProblem(name));
}

std::shared_ptr<AcousticProblem> CreateAcousticProblemShared(const string &name) {
  return std::shared_ptr<AcousticProblem>(CreateAcousticProblem(name));
}

std::shared_ptr<AcousticProblem> CreateAcousticProblemShared(const string &name,
                                                             std::shared_ptr<Meshes> meshes) {
  if (name == "FWIExampleProblem") {
    return std::make_shared<FWIExampleProblem>(std::move(meshes));
  }
  Exit(name + " Error: no Problem found!");
}
