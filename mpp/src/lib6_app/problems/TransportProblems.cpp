#include "TransportProblems.hpp"
#include "HybridMain.hpp"
#include "MixedMain.hpp"


TransportPDEProblemBuilder TransportPDEProblemBuilder::RiemannTransport1D() {
  const double c = 1.0;

  auto amplitude = [](double r) -> double {
    if (abs(r) > 0.06251) return 0.0;
    return 1.0;
  };

  return TransportPDEProblemBuilder("Interval")
    .WithName([]() { return "RiemannTransport1D"; })
    .WithSolution([amplitude, c](double t, const Cell &C, const Point &x) -> double {
      return amplitude(c * t - x[0] + 0.5);
    })
    .WithFlux([](const LevelPair &level, const Cell &c, const Point &x) -> VectorField { 
      return {1.0, 0.0};
    })
    .WithFaceNormalFlux([](const LevelPair &level, const Cell &c, int f, const VectorField &N, const Point &x) -> double {
      return VectorField(1.0, 0.0) * N;
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::InflowTest() {
  return TransportPDEProblemBuilder("Square-10x10")
    .WithName([]() { return "InflowTest"; })
    .WithSolution([](double t, const Cell &C, const Point &x) -> double { 
      if (x[0] == 0.0 && x[1] >= 1.0) return 1.0;
      return 0.0;
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::Inflow() {
  return TransportPDEProblemBuilder("Square-10x10")
    .WithName([]() { return "Inflow"; })
    .WithSolution([](double t, const Cell &C, const Point &x) -> double { 
      if (x[1] == 10) return 10.0;
      return 0.0;
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::Hat() {
  return TransportPDEProblemBuilder("Square-10x10")
    .WithName([]() { return "Hat"; })
    .WithSolution([](double t, const Cell &C, const Point &x) -> double { 
      Point midPoint = Point(5, 5);
      double rr;
      rr = pow(0.5, 3) - norm(midPoint - x);
      if (norm(midPoint - x) < 0.05)
        return rr * 40;
      return 0.0;
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::FatHat() {
  return TransportPDEProblemBuilder("Square-10x10")
    .WithName([]() { return "FatHat"; })
    .WithSolution([](double t, const Cell &C, const Point &x) -> double { 
      if (x[0] >= 1.0 && x[0] <= 4.0 && x[1] >= 1.0 && x[1] <= 4.0) return 1.0;
      else return 0.0;
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::TravelingWave() {
  VectorField Q(0.0, -15.0, 0.0);
  VectorField N(0.0, -1.0, 0.0);
  double c = N * Q;

  auto amplitude = [](double r) -> double {
    if (abs(r) < 0.625) return 1.0;
    else return 0.0;
  };

  return TransportPDEProblemBuilder("UnitSquare")
    .WithName([]() { return "TravelingWave"; })
    .WithRHS(true)
    .WithSolution([amplitude, N, c](double t, const Cell &C, const Point &x) -> double { 
      if (x[0] >= 5 && x[0] <= 15) {
        Point Pmid = Point(0.0, 21.0, 0.0);
        Point xd = x - Pmid;
        return amplitude(c * t - N[0] * xd[0] - N[1] * xd[1] - N[2] * xd[2]);
      }
      return 0.0;
    })
    .WithFlux([Q](const LevelPair &level, const Cell &c, const Point &x) -> VectorField {
      return Q;
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::SphericalWave2D() {
  const double aa = 1.0;
  const double cc = 0.5;

  return TransportPDEProblemBuilder("UnitSquare")
    .WithName([]() { return "SphericalWave2D"; })
    .WithHasExactSolution(true)
    .WithSolution([aa, cc](double t, const Cell &C, const Point &x) -> double { 
      Point midPoint = Point(t, t, 0.0);
      double r = dist(midPoint, x);
      if (r < 1 / cc)
        return aa * pow(cos(cc * Pi * r) + 1.0, 2.0);
      return 0.0;
    })
    .WithFlux([](const LevelPair &level, const Cell &c, const Point &x) -> VectorField {
      return {5.0, 5.0, 0.0};
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::CircleWave2D() {
  const double aa = 1.0;
  const double cc = 0.5;

  return TransportPDEProblemBuilder("Square-10x10")
    .WithName([]() { return "CircleWave2D"; })
    .WithHasExactSolution(true)
    .WithSolution([aa, cc](double t, const Cell &C, const Point &x) -> double { 
      Point midPoint = 5.0 * Point(cos(2. * Pi * t), sin(2. * Pi * t));
      double r = dist(midPoint, x);
      if (r < 1 / cc)
        return aa * pow(cos(cc * Pi * r) + 1.0, 2.0);
      return 0.0;
    })
    .WithFlux([](const LevelPair &level, const Cell &c, const Point &x) -> VectorField {
      return 2 * Pi * VectorField(-x[1], x[0], 0.0);
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::SineWave2D() {
  return TransportPDEProblemBuilder("UnitSquare")
    .WithName([]() { return "SineWave2D"; })
    .WithHasExactSolution(true)
    .WithSolution([](double t, const Cell &C, const Point &x) -> double { 
      return sin((x[0] + x[1] - 2.0 * t) * 1. / 5. * Pi);
    })
    .WithFlux([](const LevelPair &level, const Cell &c, const Point &x) -> VectorField {
      return {1., 1., 0};
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::CirclePacman() {
  double aa = 1.0;
  double cc = 0.5;
  double alpha = Pi / 4;

  return TransportPDEProblemBuilder("Square-10x10")
    .WithName([]() { return "CirclePacman"; })
    .WithSolution([aa, cc, alpha](double t, const Cell &C, const Point &x) -> double { 
      Point Pmid = 5.0 * Point(cos(2. * Pi * t), sin(2. * Pi * t));
      Point P = x - Pmid;
      double r = dist(Pmid, x);

      if (r < 1 / cc) {
        double MidAngle = atan2(Pmid[1], Pmid[0]);
        double angle = atan2(P[1], P[0]) - MidAngle + 8 * Pi - Pi/2;
        angle = fmod(angle, 2 * Pi) ;
        if (abs(angle) < alpha || abs(angle) > 2*Pi - alpha) return 0.0;
        return aa * pow(cos(cc * Pi * r) + 1.0, 2.0);
      }
      return 0.0;
    })
    .WithFlux([](const LevelPair &level, const Cell &c, const Point &x) -> VectorField {
      return 2 * Pi * VectorField(-x[1], x[0], 0.0);
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::Riemann1D() {
  const double thickness = 1.0 / 8.0;
  const double centerStart = 1.5 * thickness;

  auto center = [centerStart](double t) -> double {
    return centerStart + t;
  };

  return TransportPDEProblemBuilder("Interval")
    .WithHasExactSolution(true)
    .WithName([]() { return "Riemann1D"; })
    .WithSolution([thickness, center](double t, const Cell &C, const Point &x) -> double { 
      if (abs(center(t) - x[0]) < thickness / 2.0)
        return 1.0 / thickness;
      else
        return 0.0;
    })
    .WithFlux([](const LevelPair &level, const Cell &c, const Point &x) -> VectorField {
      return {1.0, 0.0, 0.0};
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::Riemann2D() {
  const double thickness = 1.0 / 8.0;
  const double centerStart = 1 - 1.5 * thickness;

  auto center = [centerStart](double t) -> double {
    return centerStart - t;
  };

  return TransportPDEProblemBuilder("UnitSquare")
    .WithName([]() { return "Riemann2D"; })
    .WithHasExactSolution(true)
    .WithSolution([thickness, center](double t, const Cell &C, const Point &x) -> double { 
      if (abs(center(t) - x[1]) < thickness / 2.0)
        return 1.0 / thickness;
      else
        return 0.0;
    })
    .WithFlux([](const LevelPair &level, const Cell &c, const Point &x) -> VectorField {
      return {0.0, -1.0, 0.0};
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::Inflow2D() {
  const double thickness = 1.0 / 8.0;
  const double centerStart = 1 - 1.5 * thickness;

  auto center = [centerStart](double t) -> double { return centerStart - t; };

  return TransportPDEProblemBuilder("UnitSquare")
      .WithName([]() { return "Inflow2D"; })
      .WithRHS([]() { return true; })
      .WithSolution(
          [thickness, center](double t, const Cell &C, const Point &x) -> double { return 0.0; })
      .WithFlux([](const LevelPair &level, const Cell &c, const Point &x) -> VectorField {
        return {0.0, -1.0, 0.0};
      })
      .WithInternalInflow([]() { return true; })
      .WithInflow([](double t, const Cell &c, const Point &x) {
        if ((sqrt(pow(x[0] - 0.5, 2) + pow(x[1] - 0.9, 2)) < 0.05))
          return 1.0; // Todo: check for dependency on timestep.
        else return 0.0;
      });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::GaussHat2D() {
  const Tensor sigma(0.005, 0.0, 0.0, 0.005);
  const double factor = 1.0 / sqrt(pow(2.0 * Pi, 2) * sigma.det());

  auto mu = [](double t) -> Point { 
    return {0.25 + t, 0.25 + t}; 
  };

  return TransportPDEProblemBuilder("UnitSquare")
    .WithName([]() { return "GaussHat2D"; })
    .WithSolution([sigma, factor, mu](double t, const Cell &C, const Point &x) -> double { 
      VectorField diff = (x - mu(t));
      VectorField temp = Invert(sigma) * diff;
      Scalar exponent = -0.5 * diff * temp;
      return factor * exp(exponent);
    })
    .WithFlux([](const LevelPair &level, const Cell &c, const Point &x) -> VectorField {
      return {0.5, 0.5};
    });
}

class LinearTransport : public ITransportProblem {
public:
  VectorField q = {1.0, 1.0, 0.0};
  VectorField a = {4.0/3.0, 2.0/3.0, 0.0};

  LinearTransport() : IProblem("UnitSquare") {}

  std::string Name() const override {
    return "LinearTransport";
  }
  VectorField Flux(const LevelPair &level, const Cell &c, const Point &x) const override {
    return q;
  }
  double divFluxVector(const LevelPair &level, const Cell &c, const Point &x) const override {
    return 0;
  }
  double Solution(double t, const Cell &c, const Point &x) const override {
    return (x - t * q) * a;
  }
  double InflowData(const LevelPair &level, const Cell &c, double time, const Point &x, VectorField &N) const override {
    return Solution(time,c,x) * Flux(level, c, x) * N;
  }
  bool RHS() const { return true; };

  bool HasExactSolution() const override { return true; }
};

class Pollution : public ITransportProblem {
  double d = 1.0 / 16.0;

  MixedMain mixedMain;

public:
  Pollution() : IProblem("UnitSquare"), mixedMain("Laplace2D") {}

  double Solution(double t, const Cell &c, const Point &x) const override {
    if ((abs(1 - x[1] - 1.5 * d) < d / 2) && (abs(0.5 - x[0]) < 3 * d))
      return 1.0;
    else return 0.0;
  }

  VectorField Flux(const LevelPair &level, const Cell &c, const Point &x) const override {
    return mixedMain.EvaluateCellFlux(level, c);
  }

  double FaceNormalFlux(const LevelPair &level, const Cell &c, int face,
                        const VectorField &N, const Point &x) const override {
    return mixedMain.EvaluateNormalFlux(level, c, face);
  }

  std::string Name() const override { return "Pollution"; }

  IProblem &GetFluxProblem() override {
    return mixedMain.GetFluxProblem();
  }

  bool HasFluxProblem() const override {
    return true;
  }

};

class PollutionSquare500 : public ITransportProblem {
  double d = 1.0 / 16.0;
  bool smooth = false;
  MixedMain mixedMain;

public:
  PollutionSquare500() : IProblem("Square500"), mixedMain("LaplaceSquare500") {
    Config::Get("smooth", smooth);
  }

  double dist_x(double x) const {
    return abs(0.5 - x) / (3 * d);
  }

  double dist_y(double t, double y) const {
    return 2 * abs(1 - y - t - 1.5 * d) / d;
  }

  double Solution(double t, const Cell &c, const Point &x) const override {
    if (dist_x(x[0]) > 1) return 0.0;
    if (dist_y(t, x[1]) > 1) return 0.0;
    if (smooth) return (1 + cos(Pi * dist_x(x[0]))) * (1 + cos(Pi * dist_y(t, x[1])));
    return 1.0 / 0.0234375;
  }

  VectorField Flux(const LevelPair &level, const Cell &c, const Point &x) const override {
    return mixedMain.EvaluateCellFlux(level, c);
  }

  double FaceNormalFlux(const LevelPair &level,
                        const Cell &c,
                        int face,
                        const VectorField &N, const Point &x) const override {
    return mixedMain.EvaluateNormalFlux(level, c, face);
  }

  std::string Name() const override { return "PollutionSquare500"; }

  IProblem &GetFluxProblem() override {
    return mixedMain.GetFluxProblem();
  }

  bool HasFluxProblem() const override {
    return true;
  }
};

class PollutionSquare501 : public ITransportProblem {
    double d = 1.0 / 16.0;
    bool smooth = false;
    MixedMain mixedMain;

public:
    PollutionSquare501() : IProblem("Square501"), mixedMain("LaplaceSquare501") {
        Config::Get("smooth", smooth);
    }

    double dist_x(double x) const {
        return abs(0.5 - x) / (3 * d);
    }

    double dist_y(double t, double y) const {
        return 2 * abs(1 - y - t - 1.5 * d) / d;
    }

    double Solution(double t, const Cell &c, const Point &x) const override {
        if (dist_x(x[0]) > 1) return 0.0;
        if (dist_y(t, x[1]) > 1) return 0.0;
        if (smooth) return (1 + cos(Pi * dist_x(x[0]))) * (1 + cos(Pi * dist_y(t, x[1])));
        return 1.0 / 0.0234375;
    }

    VectorField Flux(const LevelPair &level, const Cell &c, const Point &x) const override {
        return mixedMain.EvaluateCellFlux(level, c);
    }

    double FaceNormalFlux(const LevelPair &level,
                          const Cell &c,
                          int face,
                          const VectorField &N, const Point &x) const override {
        return mixedMain.EvaluateNormalFlux(level, c, face);
    }

    std::string Name() const override { return "PollutionSquare501"; }

    IProblem &GetFluxProblem() override {
        return mixedMain.GetFluxProblem();
    }

    bool HasFluxProblem() const override {
        return true;
    }
};

class PollutionSquare2 : public ITransportProblem {
  double d = 1.0 / 16.0;
  bool smooth = false;
  MixedMain mixedMain;

public:
  PollutionSquare2() : IProblem("Square2"), mixedMain("LaplaceSquare2") {
    Config::Get("smooth", smooth);
  }

  bool HasExactSolution() const override { return true; }

  double dist_x(double x) const {
    return abs(0.5 - x) / (3 * d);
  }

  double dist_y(double t, double y) const {
    return 2 * abs(1 - y - t - 1.5 * d) / d;
  }

  double Solution(double t, const Cell &c, const Point &x) const override {
    // return -x[0] - t;
    if (dist_x(x[0]) > 1) return 0.0;
    if (dist_y(t, x[1]) > 1) return 0.0;
    if (smooth) return (1 + cos(Pi * dist_x(x[0]))) * (1 + cos(Pi * dist_y(t, x[1])));
    return 1.0 / 0.0234375;
  }

  VectorField Flux(const LevelPair &level, const Cell &c, const Point &x) const override {
    return mixedMain.EvaluateCellFlux(level, c);
  }

  double FaceNormalFlux(const LevelPair &level, const Cell &c, int face,
                        const VectorField &N, const Point &x) const override {
    return mixedMain.EvaluateNormalFlux(level, c, face);
  }

  std::string Name() const override {
    if (smooth) return "PollutionSquare2 (smmoth)";
    return "PollutionSquare2";
  }

  IProblem &GetFluxProblem() override {
    return mixedMain.GetFluxProblem();
  }

  bool HasFluxProblem() const override {
    return true;
  }
};

class GaussHat1D : public ITransportProblem {
private:
  double sigma2;
  double factor;

public:
  GaussHat1D(double sigma = 0.05) : IProblem("Interval"),
                                    sigma2(sigma * sigma) {
    // TODO: What is "T"? Does this have to be here?
    T = 0.5;
    factor = 1.0 / sqrt(2.0 * Pi * sigma2);
  }

  static Point mu(double t) {
    return {0.25 + t, 0.0};
  }

  Scalar Solution(double t, const Cell &c, const Point &x) const override {
    VectorField diff = (x - mu(t));
    Scalar exponent = -diff * diff / (2 * sigma2);
    return factor * exp(exponent);
  }

  static VectorField TransportFlux(const Point &x) {
    return {1.0, 0.0};
  }

  VectorField Flux(const LevelPair &level, const Cell &c, const Point &x) const override {
    return VectorField{1.0, 0.0};
  }

  string Name() const override { return "GaussHat1D"; }
};

TransportPDEProblemBuilder TransportPDEProblemBuilder::CosHat1D() {
  const double amplitude = 1.0;
  const double frequency = 8.0;

  return TransportPDEProblemBuilder("Interval")
    .WithName([]() { return "CosHat1D"; })
    .WithSolution([amplitude, frequency](double t, const Cell &C, const Point &x) -> double { 
      Point midPoint = Point(0.25 + t, 0.0);
      double r = dist(midPoint, x);
      if (r < (Pi / 2.0) / frequency)
        return amplitude * pow(cos(frequency * r), 2.0);
      return 0.0;
    })
    .WithFlux([](const LevelPair &level, const Cell &c, const Point &x) -> VectorField {
      return {1.0, 0.0};
    });
}

TransportPDEProblemBuilder TransportPDEProblemBuilder::CosHat2D() {
  const double amplitude = 1.0;
  const double frequency = 8.0;

  return TransportPDEProblemBuilder("UnitSquare")
    .WithName([]() { return "CosHat2D"; })
    .WithSolution([amplitude, frequency](double t, const Cell &C, const Point &x) -> double { 
      Point midPoint = 0.25 * (Point(cos(2.0 * Pi * t), sin(2.0 * Pi * t)))
                     + Point(0.5, 0.5);
      double r = dist(midPoint, x);
      if (r < (Pi / 2.0) / frequency)
        return amplitude * pow(cos(frequency * r), 2.0);
      return 0.0;
    })
    .WithFlux([](const LevelPair &level, const Cell &c, const Point &x) -> VectorField {
      return 2 * Pi * VectorField(0.5 - x[1], x[0] - 0.5, 0.0);
    });
}

ITransportProblem *CreateTransportProblem(const string &problemName) {
  // TODO how should flux be implemented for the commented out classes?
  if (problemName == "Riemann1D")
    return TransportPDEProblemBuilder::Riemann1D().Build();
  else if (problemName == "Riemann2D" || problemName == "RiemannUnitSquare")
    return TransportPDEProblemBuilder::Riemann2D().Build();
  else if (problemName == "Inflow2D")
    return TransportPDEProblemBuilder::Inflow2D().Build();
  else if (problemName == "CircleWave2D")
    return TransportPDEProblemBuilder::CircleWave2D().Build();
    //else if (problemName == "Inflow")
    //  return new Inflow();
  else if (problemName == "SphericalWave2D")
    return TransportPDEProblemBuilder::SphericalWave2D().Build();
  else if (problemName == "SineWave2D")
    return TransportPDEProblemBuilder::SineWave2D().Build();
    /*else if (problemName == "FatHat")
      return new FatHat();*/
  else if (problemName == "TravelingWave")
    return TransportPDEProblemBuilder::TravelingWave().Build();
  /*else if (problemName == "InflowTest")
    return new InflowTest();
  else if (problemName == "Hat")
    return new Hat();*/
  if (problemName == "CosHat2D")
    return TransportPDEProblemBuilder::CosHat2D().Build();
  else if (problemName == "CirclePacman")
    return TransportPDEProblemBuilder::CirclePacman().Build();
  else if (problemName == "RiemannTransport1D")
    return TransportPDEProblemBuilder::RiemannTransport1D().Build();
  else if (problemName == "GaussHat2D")
    return TransportPDEProblemBuilder::GaussHat2D().Build();
  else if (problemName == "GaussHat1D")
    return new GaussHat1D();
  else if (problemName == "CosHat1D")
    return TransportPDEProblemBuilder::CosHat1D().Build();
  else if (problemName == "Pollution") {
    return new Pollution();
  } else if (problemName == "LinearTransport") {
    return new LinearTransport();
  } else if (problemName == "PollutionSquare500") {
    return new PollutionSquare500();
  } else if (problemName == "PollutionSquare2") {
    return new PollutionSquare2();
  }
  else if (problemName == "PollutionSquare501") {
      return new PollutionSquare501();
  } else Exit(problemName + " not found")
}

std::unique_ptr<ITransportProblem>
CreateTransportProblemUnique(const std::string &problemName) {
  return std::unique_ptr<ITransportProblem>(CreateTransportProblem(problemName));
}

std::shared_ptr<ITransportProblem>
CreateTransportProblemShared(const std::string &problemName) {
  return std::shared_ptr<ITransportProblem>(CreateTransportProblem(problemName));
}
