#include "ViscoElasticProblems.hpp"
#include "CImg.hpp"


class RickerWavelet_VE : public TProblemViscoElastic {
public:
  Point Pmid;
  double width;
  double pulse_duration;
  double source_factor;
  int linear_weight;

  RickerWavelet_VE()
    : TProblemViscoElastic(2),
      Pmid(-1.0, 1.0, 0.0, 0.4), width(0.2), pulse_duration(0.2),
      source_factor(100.0),
      linear_weight(1) {
    Config::Get("width", width);
    Config::Get("source_factor", source_factor);
    Config::Get("source_duration", pulse_duration);
    Config::Get("linear_weight", linear_weight);
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

  virtual double F(const Point &x, COMPONENT comp) const override {
    if (comp == COMPONENT::V_X)
      return source_factor * TimeFunction(x.t() - Pmid.t())
        * SpaceFunction(dist(Pmid, x));
    return 0.0;
  }

  double ut(const Point &x, COMPONENT comp) const override { return 0.0; }

  double Rho(const Point &x) const override {
    if (abs(x[0]) + abs(x[1]) < 0.5) return 2.0;
    return 1.0;
  }

  virtual double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const {
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch(linear_weight) {
      case 1: {
        Point ba = b - a;
        if (comp != COMPONENT::V_X) return 0;
        else return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }
};

class WaveProblemFromImages2d_VE_RHS : public TProblemViscoElastic {
protected:
  double left, right, top, bottom;
  double min_rho, max_rho, min_kappa, max_kappa, min_mu, max_mu;
  int problemLevel;
  double mlLength;
  short *RhoPixel;
  short *KappaPixel;
  short *MuPixel;
  int RhoWidth;
  int RhoHeight;
  int KappaWidth;
  int KappaHeight;
  int MuWidth;
  int MuHeight;
  double pow2L1;
  double pow2L;

public:
  std::string Name() const override { return "WaveProblemFromImages2d_VE_RHS"; }

  WaveProblemFromImages2d_VE_RHS(const Point &minP, const Point &maxP)
    : TProblemViscoElastic(2), left(minP[0]), right(maxP[0]), top(maxP[1]),
      bottom(minP[1]), min_rho(1.0), max_rho(3.0), min_mu(1.0), max_mu(1.0),
      min_kappa(1.0), max_kappa(3.0), problemLevel(-1), mlLength(1.0) {
    mout << "WaveProblemFromImages2d_VE_RHS: mesh bounding box:"
         << "(" << left << ", " << top << ")" << " - "
         << "(" << right << ", " << bottom << ")" << endl;

    string rho_image = "marmousi_rho";
    string mu_image = "marmousi_mu";
    string kappa_image = "marmousi_kappa";
    Config::Get("ModelImageRho", rho_image);
    Config::Get("ModelImageMu", mu_image);
    Config::Get("ModelImageKappa", kappa_image);

    if (PPM->Master(0)) {
      cimg_library::CImg<unsigned char> ImgRho, ImgKappa, ImgMu;
      std::stringstream rho_pp, mu_pp, kappa_pp;
      rho_pp << string(ProjectMppDir) + "/models/" << rho_image << ".png";
      kappa_pp << string(ProjectMppDir) + "/models/" << kappa_image << ".png";
      mu_pp << string(ProjectMppDir) + "/models/" << mu_image << ".png";
      ImgRho.load_png(rho_pp.str().c_str());
      ImgKappa.load_png(kappa_pp.str().c_str());
      ImgMu.load_png(mu_pp.str().c_str());
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

      MuWidth = ImgMu.width();
      MuHeight = ImgKappa.height();
      MuPixel = new short[MuWidth * MuHeight];
      for (int i = 0; i < MuWidth; ++i)
        for (int j = 0; j < MuHeight; ++j)
          MuPixel[i * MuHeight + j] = ImgMu(i, j);

      PPM->BroadcastInt(RhoWidth);
      PPM->BroadcastInt(RhoHeight);
      PPM->Broadcast(RhoPixel, (RhoWidth * RhoHeight) * sizeof(short), 0);
      PPM->BroadcastInt(KappaWidth);
      PPM->BroadcastInt(KappaHeight);
      PPM->Broadcast(KappaPixel, (KappaWidth * KappaHeight) * sizeof(short), 0);
      PPM->BroadcastInt(MuWidth);
      PPM->BroadcastInt(MuHeight);
      PPM->Broadcast(MuPixel, (MuWidth * MuHeight) * sizeof(short), 0);
    } else {
      RhoWidth = PPM->BroadcastInt();
      RhoHeight = PPM->BroadcastInt();
      RhoPixel = new short[RhoWidth * RhoHeight];
      PPM->Broadcast(RhoPixel, (RhoWidth * RhoHeight) * sizeof(short), 0);
      KappaWidth = PPM->BroadcastInt();
      KappaHeight = PPM->BroadcastInt();
      KappaPixel = new short[KappaWidth * KappaHeight];
      PPM->Broadcast(KappaPixel, (KappaWidth * KappaHeight) * sizeof(short), 0);
      MuWidth = PPM->BroadcastInt();
      MuHeight = PPM->BroadcastInt();
      MuPixel = new short[MuWidth * MuHeight];
      PPM->Broadcast(MuPixel, (MuWidth * MuHeight) * sizeof(short), 0);
    }

    mout << "WaveProblemFromImages2d_VE_RHS:: rho model read   '"
         << rho_image << "' (size: " << RhoWidth << "x"
         << RhoHeight << ")" << endl;
    mout << "WaveProblemFromImages2d_VE_RHS:: kappa model read '"
         << kappa_image << "' (size: " << KappaWidth << "x"
         << KappaHeight << ")" << endl;
    mout << "WaveProblemFromImages2d_VE_RHS:: mu model read   '"
         << mu_image << "' (size: " << MuWidth << "x"
         << MuHeight << ")" << endl;

    Config::Get("MinRho", min_rho);
    Config::Get("MaxRho", max_rho);
    Config::Get("MinKappa", min_kappa);
    Config::Get("MaxKappa", max_kappa);
    Config::Get("MinMu", min_mu);
    Config::Get("MaxMu", max_mu);
    Config::Get("ProblemLevel", problemLevel);
    if (problemLevel != -1) Config::Get("mlLength", mlLength);

    pow2L1 = pow(0.5, problemLevel + 1);
    pow2L = pow(2, problemLevel);
  }

  int toImgCoordHorizontal(double c, int height, int width) const {
    if (problemLevel < 0)
      return int(round((c - left) / (right - left) * width));

    c /= mlLength;
    double ganzzahl = int(c);
    double dezimal = c - ganzzahl;
    for (int k = 0; k < pow2L; k++)
      if (abs(dezimal - (1 + 2 * k) * pow2L1) <= pow2L1)
        return int(round(
          ((ganzzahl + (1 + 2 * k) * pow2L1) * mlLength - left) / (right - left)
            * width));
    Exit("OHAAAA")
  }

  int toImgCoordVertical(double c, int height, int width) const {
    if (problemLevel < 0)
      return int(round((c - bottom) / (top - bottom) * height));

    c /= mlLength;
    double ganzzahl = int(c);
    double dezimal = c - ganzzahl;
    for (int k = 0; k < pow2L; k++)
      if (abs(dezimal - (1 + 2 * k) * pow2L1) < pow2L1)
        return int(round(((ganzzahl + (1 + 2 * k) * pow2L1) * mlLength - bottom)
                           / (top - bottom) * height));
    Exit("OHAAAA");
  }

  double Rho(const Point &P) const {
    double x = P[0], y = P[1];
    int ix = toImgCoordHorizontal(x, RhoHeight, RhoWidth);
    int iy = toImgCoordVertical(y, RhoHeight, RhoWidth);
    return min_rho + RhoPixel[ix * RhoHeight + iy] * (max_rho - min_rho) / 255;
  }

  double Mu(const Point &P, int var = -1) const {
    double x = P[0], y = P[1];
    int ix = toImgCoordHorizontal(x, MuHeight, MuWidth);
    int iy = toImgCoordVertical(y, MuHeight, MuWidth);
    const double Vs = min_mu + MuPixel[ix * MuHeight + iy] * (max_mu - min_mu) / 255;
    const double mu = Vs * Vs * Rho(P) / (1 + numL * tau_s);
    if (var == -1) return mu * (1 + numL * tau_s);
    if (var == 1) return mu;
    else return mu * tau_s;
  }

  double Lambda(const Point &P, int var = -1) const {
    double x = P[0], y = P[1];
    int ix = toImgCoordHorizontal(x, KappaHeight, KappaWidth);
    int iy = toImgCoordVertical(y, KappaHeight, KappaWidth);
    const double Vp =
      min_kappa + KappaPixel[ix * KappaHeight + iy] * (max_kappa - min_kappa) / 255;

    const double mu = Mu(P, var);
    const double lambda = (Vp * Vp * Rho(P) - 2 / 3 * mu) / (1 + numL * tau_p);

    if (var == -1) return lambda * (1 + numL * tau_p);
    if (var == 1) return lambda;
    else return lambda * tau_p;
  }

  double Kappa(const Point &P, int var = -1) const {
    Exit("nope");
    double x = P[0], y = P[1];
    int ix = toImgCoordHorizontal(x, KappaHeight, KappaWidth);
    int iy = toImgCoordVertical(y, KappaHeight, KappaWidth);
    const double Vp =
      min_kappa + KappaPixel[ix * KappaHeight + iy] * (max_kappa - min_kappa) / 255;
    const double kappa = Vp * Vp * Rho(P) / (1 + numL * tau_p);
    if (var == -1) return kappa * (1 + numL * tau_p);
    if (var == 1) return kappa;
    else return kappa * tau_p;
  }
};

class Marmousi2_VE_rhs : public WaveProblemFromImages2d_VE_RHS {
public:
  Point Pmid;
  double width;
  double pulse_duration;
  double source_factor;
  int linear_weight;

  Marmousi2_VE_rhs()
    : WaveProblemFromImages2d_VE_RHS(Point(-4000, 0), Point(13000, 3500)),
      Pmid(0.0, 0.0, 0.0, 0.0), pulse_duration(0.1), source_factor(100), width(100),
      linear_weight(1) {
    Config::Get("ProblemMid", Pmid);
    Config::Get("width", width);
    Config::Get("source_factor", source_factor);
    Config::Get("source_duration", pulse_duration);
    Config::Get("linear_weight", linear_weight);

    tau = pulse_duration / (2 * M_PI);
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

  bool HasRHS() const { return true; }

  virtual double F(const Point &x, int i) const {
    if (i >= dim && i < 2 * dim)
      return source_factor * TimeFunction(x.t() - Pmid.t())
        * SpaceFunction(dist(Pmid, x));
    return 0.0;
  }

  virtual double ut(const Point &x, COMPONENT comp) const { return 0.0; }

  virtual double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const {
    Exit("Todo components")
    int i = 0;
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch(linear_weight) {
      case 1: {
        Point ba = b - a;
        if (i < dim) return 0;
        if (i == 4 || i == 7 || i == 10) return 0;
        else return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }

  double C_s(const Point &x) const { return sqrt(Mu(x) / Rho(x)); }

  double C_p(const Point &x) const { return sqrt((2.0 * Mu(x) + Lambda(x)) / Rho(x)); }
};

class ViscoElastic_Spherical_Wave_2D : public TProblemViscoElastic {
  Point Pmid;
  int linear_weight;
  double width;
public:
  ViscoElastic_Spherical_Wave_2D() :
    TProblemViscoElastic(2, 1), Pmid(0, 0, 0), width(0.25) {
    Config::Get("ProblemMid", Pmid);
    Config::Get("linear_weight", linear_weight);
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

  virtual double ut(const Point &x, COMPONENT comp) const override {
    if (x.t() > 0) return 0.0;
    else return u0(x, comp);
  }

  virtual double ut(const Point &x, const Cell &tc, COMPONENT comp) const override {
    return ut(x, comp);
  }

  virtual double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const override {
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch(linear_weight) {
      case 1: {
        Point ba = b - a;
        if (comp != COMPONENT::V_X) return 0.0;
        else return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }
};

class Elasticity_Constant : public TProblemElasticity {
  int linear_weight;
public:
  Elasticity_Constant()
    : TProblemElasticity(2) {
    Config::Get("linear_weight", linear_weight);
  }

  virtual double ut(const Point &x, const cell &tc, int i) const {
    return ut(x, i);
  }

  double ut(const Point &x, int i) const {
    switch(i) {
      case 3 :
        return 2.0;
      case 4 :
        return 3.0;
      default:
        return 0.0;
    }
  }

  virtual double ut_dual(const Point &x, int i, const Point &a, const Point &b) const {
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch(linear_weight) {
      case 1: {
        Point ba = b - a;
        if (i != 3) return 0;
        else return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }

  bool HasExactSolution() const { return true; }
};

class Elasticity_Linear_2D : public TProblemElasticity {
  int linear_weight;
public:
  Elasticity_Linear_2D()
    : TProblemElasticity(2) {
    Config::Get("linear_weight", linear_weight);
  }

  virtual double ut(const Point &x, const cell &tc, int i) const {
    return ut(x, i);
  }

  double ut(const Point &x, int i) const {
    switch(i) {
      case 3 :
        return x[0];
      default:
        return 0.0;
    }
  }

  bool HasRHS() const { return true; }

  double F(const Point &x, int i) const {
    switch(i) {
      case 0 :
        return -1.0;
      default:
        return 0.0;
    }
  }

  virtual double ut_dual(const Point &x, int i, const Point &a, const Point &b) const {
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch(linear_weight) {
      case 1: {
        Point ba = b - a;
        if (i != 3) return 0;
        else return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }

  bool HasExactSolution() const { return true; }
};

class Spherical_Wave_2D : public TProblemElasticity {
  Point Pmid;
  int linear_weight;
  double width;
public:
  Spherical_Wave_2D()
    : TProblemElasticity(2), Pmid(0, 0, 0), width(0.25) {
    Config::Get("ProblemMid", Pmid);
    Config::Get("linear_weight", linear_weight);
  }

  virtual double ut(const Point &x, const cell &tc, int i) const {
    return ut(x, i);
  }

  double A0(double s) const {
    if (s < width) return pow(cos(s * 0.5 * M_PI / width), 6.0);
    else return 0.0;
  }

  double A(const Point &x) const {
    double r = dist(Pmid, x);
    return A0(r * r);
  }

  double u0(const Point &x, int i) const {
    const double tmp = A(x);
    switch(i) {
      case 3 :
        return (x[0] - Pmid[0]) * tmp;
      case 4 :
        return (x[1] - Pmid[1]) * tmp;
      default:
        return 0.0;
    }
  }

  virtual double ut(const Point &x, int i) const {
    if (x.t() > 0) return 0.0;
    else return u0(x, i);
  }

  virtual double ut_dual(const Point &x, int i, const Point &a, const Point &b) const {
    if (x[0] < a[0] || x[1] < a[1] || x[0] > b[0] || x[1] > b[1]) return 0.0;
    switch(linear_weight) {
      case 1: {
        Point ba = b - a;
        if (i != 3) return 0;
        else return 1.0 / (ba[0] * ba[1]);
      }
      default: Exit("Not implemented");
    }
  }
};

class Plane_Wave_2D : public TProblemElasticity {
public:
  Plane_Wave_2D()
    : TProblemElasticity(2) {}

  virtual double ut(const Point &x, const cell &tc, int i) const {
    if (x.t() > 0)
      return 0.0;
    else if (x[0] < 5.0) {
      switch(i) {
        case 0 :
          return -1.0;
        case 3 :
          return 1.0;
        default:
          return 0.0;
      }
    } else return 0.0;
  }
};

class Elasticity_Plane_Wave_2D : public TProblemElasticity {
  double lambda;
  double mu;
  double rho;
  Point Pmid;
  int linear_weight;
public:
  Elasticity_Plane_Wave_2D()
    : TProblemElasticity(2), Pmid(0, 0, 0) {
    Config::Get("lambda", lambda);
    Config::Get("mu", mu);
    Config::Get("rho", rho);
    Config::Get("ProblemMid", Pmid);
    Config::Get("linear_weight", linear_weight);
  }

  double Lambda(const Point &x) const override { return lambda; }

  double Mu(const Point &x) const override { return mu; }

  double Rho(const Point &x) const override { return rho; }

  double C_s(const Point &x) const { return sqrt(Mu(x) / Rho(x)); }

  double C_p(const Point &x) const { return sqrt((2.0 * Mu(x) + Lambda(x)) / Rho(x)); }

  virtual double ut(const Point &x, const Cell &tc, COMPONENT comp) const {
    if (x.t() > 0) return 0.0;

    double r = abs(Pmid[1] - x[1]);
    switch(comp) {
      case COMPONENT::S0_x:
        if (r < 1) return -exp(-r * r) * (1 - r * r);
        else return 0.0;
      case COMPONENT::S0_y:
        if (r < 1) return exp(-r * r) * (1 - r * r);
        else return 0.0;
      default:
        return 0.0;
    }
  }

  virtual double ut_dual(const Point &x, COMPONENT comp, const Point &a, const Point &b) const {
    switch(linear_weight) {
      case 1: {
        Point ba = b - a;
        if (comp == COMPONENT::S0_x) return 0;
        if (comp == COMPONENT::S0_y) return 1.0 / (ba[0] * ba[1]);
        if (comp == COMPONENT::S0_xy) return 1.0 / (ba[0] * ba[1]);
      }
      case 2: {
        Point PPmid = a + 0.5 * (b - a);
        double r = abs(PPmid[1] - x[1]) / abs(PPmid[1] - b[1]);
        switch(comp) {
          case COMPONENT::S0_y:
            if (r < 1)
              return exp(-r * r) * (1 - r * r);
            else
              return 0.0;
          case COMPONENT::S0_xy :
            if (r < 1)
              return exp(-r * r) * (1 - r * r);
            else
              return 0.0;
          default:
            return 0.0;
        }
        Exit("not implemented");
      }
      default: Exit("Not implemented");
    }
  }
};

TProblemElasticity *CreateElasticProblem(const string &name) {
  if (name == "Linear") return new Elasticity_Linear_2D();
  if (name == "Constant") return new Elasticity_Constant;
  if (name == "SphericalWave2D") return new Spherical_Wave_2D();
  if (name == "PlaneWave2D") return new Elasticity_Plane_Wave_2D();
  else Exit("Error: no wave found!");
}

std::unique_ptr<TProblemElasticity> CreateElasticProblemUnique(const string &name) {
  return std::unique_ptr<TProblemElasticity>(CreateElasticProblem(name));
}

TProblemViscoElastic *CreateViscoElasticProblem(const string &name) {
  if (name == "SphericalWave2D") return new ViscoElastic_Spherical_Wave_2D();
  if (name == "RickerWavelet") return new RickerWavelet_VE();
  if (name == "Marmousi2_rhs") return new Marmousi2_VE_rhs();
  else Exit("Error: no wave found!");
}

std::unique_ptr<TProblemViscoElastic> CreateViscoElasticProblemUnique(const string &name) {
  return std::unique_ptr<TProblemViscoElastic>(CreateViscoElasticProblemUnique(name));
}
