#include "EllipticProblems.hpp"

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::Laplace1D() {
  return EllipticPDEProblemBuilder()
    .WithMeshName("Interval")
    .WithHasExactSolution([]() { return true; })
    .WithLoad([](const Point &x, const cell &c) { return 0.0; })
    .WithSolution([](const Point &x) { return 1 - x[0]; })
    .WithFlux([](const Point &x) -> VectorField { return { -1.0, 0.0 }; })
    .WithPermeability([](const Point &x) { return One; })
    .WithName([]() { return "Laplace1D";} );
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::Discontinuous1D() {
  double permeability = 10.0;
  Config::Get("Permeability", permeability);

  return EllipticPDEProblemBuilder()
    .WithMeshName("Interval")
    .WithName([]() { return "Discontinuous1D"; })
    .WithLoad([](const Point &x, const cell &c) { return 0.0; })
    .WithSolution([](const Point &x) { return 0.0; })
    .WithFlux([](const Point &x) -> VectorField { return {-1.0, 0.0 }; })
    .WithPermeability([permeability](const Point &x) {
      if (dist(x, Point(0.5, 0.0)) < 0.2)
        return permeability * One;
      return One;
    });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::Laplace() {
  return EllipticPDEProblemBuilder()
    .WithHasExactSolution([]() { return true; })
    .WithLoad([](const Point &x, const cell &c) { return 0.0; })
    .WithSolution([](const Point &x) { return -x[1]; })
    .WithFlux([](const Point &x) -> VectorField { return {0.0, -1.0 }; })
    .WithPermeability([](const Point &x) { return One; })
    .WithName([]() { return "Laplace"; });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::Laplace2D() {
  return EllipticPDEProblemBuilder::Laplace()
    .WithMeshName("UnitSquare")
    .WithName([]() { return "Laplace2D"; });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::LaplaceSquare500() {
  return EllipticPDEProblemBuilder::Laplace()
    .WithMeshName("Square500")
    .WithSolution([](const Point &x) { return -x[1]; })
    .WithFlux([](const Point &x) -> VectorField { return {0.0, -1.0 }; })
    .WithName([]() { return "LaplaceSquare500"; });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::LaplaceSquare501() {
  return EllipticPDEProblemBuilder::Laplace()
    .WithMeshName("Square501")
    .WithSolution([](const Point &x) { return -x[1]; })
    .WithFlux([](const Point &x) -> VectorField { return {0.0, -1.0 }; })
    .WithName([]() { return "LaplaceSquare501"; });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::LaplaceSquare2() {
  return EllipticPDEProblemBuilder::Laplace()
    .WithMeshName("Square2")
    .WithHasExactSolution([]() { return true; })
    .WithSolution([](const Point &x) { return -x[1]; })
    .WithFlux([](const Point &x) -> VectorField { return {0.0, -1.0 }; })
    .WithName([]() { return "LaplaceSquare2"; });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::Discontinuous() {
  double permeability = 10.0;
  Config::Get("Permeability", permeability);

  return EllipticPDEProblemBuilder()
    .WithLoad([](const Point &x, const cell &c) { return 0.0; })
    .WithSolution([](const Point &x) { return 0.0; })
    .WithFlux([](const Point &x) -> VectorField { return { 0.0, -1.0 }; })
    .WithPermeability([permeability](const Point &x) {
      if (dist(x, Point(0.5, 0.5)) < 0.2)
        return permeability * One;
      return One;
    })
    .WithName([]() { return "Discontinuous"; });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::Discontinuous2D() {
  return EllipticPDEProblemBuilder::Discontinuous()
    .WithMeshName("UnitSquare")
    .WithName([]() { return "Discontinuous2D"; });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::DiscontinuousSquare500() {
  return EllipticPDEProblemBuilder::Discontinuous()
    .WithMeshName("Square500")
    .WithName([]() { return "DiscontinuousSquare500"; });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::Kellogg() {
  return EllipticPDEProblemBuilder()
    .WithMeshName("Square-1x1")
    .WithName([]() { return "Kellogg"; })
    .WithLoad([](const Point &x, const cell &c) { return 0; })
    .WithSolution([](const Point &x) {
      Scalar r = sqrt(x[0] * x[0] + x[1] * x[1]);
      double phi = atan2(x[1], x[0]);
      const double alpha = 0.1;
      const double beta = -14.92256510455152;
      if (phi < 0) phi += 2 * Pi;
      if (phi < 0.5 * Pi)
        return pow(r, alpha) * cos((0.5 * Pi - beta) * alpha)
          * cos(phi * alpha);
      else if (phi < Pi)
        return pow(r, alpha) * cos(0.5 * Pi * alpha)
          * cos((phi - Pi + beta) * alpha);
      else if (phi < 1.5 * Pi)
        return pow(r, alpha) * cos(beta * alpha)
          * cos((phi - 1.5 * Pi) * alpha);
      else
        return pow(r, alpha) * cos((phi - 1.5 * Pi - beta) * alpha);
    })
    .WithFlux([](const Point &x) -> VectorField {
      return {0.0, -1.0};
    })
    .WithPermeability([](const Point &x) {
      if (x[0] * x[1] > 0) return 161.4476387975881 * One;
      return One;
    });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::Rock() {
  double permeability = 1.0;
  Config::Get("Permeability", permeability);

  return EllipticPDEProblemBuilder()
    .WithMeshName("Rock")
    .WithName([]() { return "Rock"; })
    .WithLoad([](const Point &x, const cell &c) { return 0.0; })
    .WithSolution([](const Point &x) { return 0.0; })
    .WithFlux([](const Point &x) -> VectorField {
      return {0.0, -1.0};
    })
    .WithPermeability([permeability](const Point &x) {
      if (x[0] <= 1 && x[1] <= 1) return permeability * One;
      return One;
    });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::Divergent() {
  return EllipticPDEProblemBuilder()
    .WithMeshName("UnitSquare")
    .WithName([]() { return "Divergent"; })
    .WithLoad([](const Point &x, const cell &c) { return 0.0; })
    .WithSolution([](const Point &x) { return -x[1]; })
    .WithFlux([](const Point &x) -> VectorField {
      return {0.0, -1.0};
    })
    .WithPermeability([](const Point &x) {
      int p = (int) (10000.0 * x[0] + 7625 * x[1]);
      p = p % 7;
      return pow(7.7, p) * One;
    });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::LaplaceFicheraCube() {
  return EllipticPDEProblemBuilder()
    .WithMeshName("FicheraCube")
    .WithName([]() { return "LaplaceFicheraCube"; })
    .WithHasExactSolution([]() { return true; })
    .WithLoad([](const Point &x, const cell &c) { return -0.75 * pow(x * x, -0.75); })
    .WithSolution([](const Point &x) { return pow(x * x,0.25); })
    .WithFlux([](const Point &x) -> VectorField {
      double q = pow(x * x,-0.75);
      VectorField flux;
      flux[0] = 0.5 * x[0] * q;
      flux[1] = 0.5 * x[1] * q;
      flux[2] = 0.5 * x[2] * q;
      return flux;
    })
    .WithPermeability([](const Point &x) {
      return One;
    });
}

EllipticPDEProblemBuilder EllipticPDEProblemBuilder::LinearAffineInflow() {
  double rate = -1.0;
  double offset = -1.0;
  Config::Get("InflowRate", rate);
  Config::Get("InflowOffset", offset);
  return EllipticPDEProblemBuilder()
    .WithName([]() { return "Linear Affine Inflow Problem"; })
    .WithFlux([rate,offset](const Point &x) -> VectorField {
      return VectorField(0.0, rate * x[0] + offset);
    });
}

std::unique_ptr<IEllipticProblem> CreateEllipticProblem(const string &problemName) {
  if (problemName == "Rock") return EllipticPDEProblemBuilder::Rock().Build();
  else if (problemName == "Kellogg") return EllipticPDEProblemBuilder::Kellogg().Build();
  else if (problemName == "Laplace") return EllipticPDEProblemBuilder::Laplace().Build();
  else if (problemName == "Laplace1D") return EllipticPDEProblemBuilder::Laplace1D().Build();
  else if (problemName == "Divergent") return EllipticPDEProblemBuilder::Divergent().Build();
  else if (problemName == "LinearAffineInflow") return EllipticPDEProblemBuilder::LinearAffineInflow().Build();
  else if (problemName == "Discontinuous") return EllipticPDEProblemBuilder::Discontinuous().Build();
  else if (problemName == "Discontinuous1D") return EllipticPDEProblemBuilder::Discontinuous1D().Build();
  else if (problemName == "P0Test1D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<0, INTERVAL>().Build();
  else if (problemName == "P1Test1D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<1, INTERVAL>().Build();
  else if (problemName == "P2Test1D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<2, INTERVAL>().Build();
  else if (problemName == "P3Test1D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<3, INTERVAL>().Build();
  else if (problemName == "P4Test1D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<4, INTERVAL>().Build();
#if SpaceDimension >= 2
  else if (problemName == "Discontinuous2D") return EllipticPDEProblemBuilder::Discontinuous2D().Build();
  else if (problemName == "DiscontinuousSquare500") return EllipticPDEProblemBuilder::DiscontinuousSquare500().Build();
  else if (problemName == "Laplace2D") return EllipticPDEProblemBuilder::Laplace2D().Build();
  else if (problemName == "LaplaceSquare500") return EllipticPDEProblemBuilder::LaplaceSquare500().Build();
  else if (problemName == "LaplaceSquare501") return EllipticPDEProblemBuilder::LaplaceSquare501().Build();
  else if (problemName == "LaplaceSquare2") return EllipticPDEProblemBuilder::LaplaceSquare2().Build();
  else if (problemName == "P0Test2D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<0, QUADRILATERAL>().Build();
  else if (problemName == "P1Test2D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<1, QUADRILATERAL>().Build();
  else if (problemName == "P2Test2D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<2, QUADRILATERAL>().Build();
  else if (problemName == "P3Test2D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<3, QUADRILATERAL>().Build();
  else if (problemName == "P4Test2D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<4, QUADRILATERAL>().Build();
  else if (problemName == "P0Test2DTet") return EllipticPDEProblemBuilder::EllipticPolynomialTest<0, TRIANGLE>().Build();
  else if (problemName == "P1Test2DTet") return EllipticPDEProblemBuilder::EllipticPolynomialTest<1, TRIANGLE>().Build();
  else if (problemName == "P2Test2DTet") return EllipticPDEProblemBuilder::EllipticPolynomialTest<2, TRIANGLE>().Build();
  else if (problemName == "P3Test2DTet") return EllipticPDEProblemBuilder::EllipticPolynomialTest<3, TRIANGLE>().Build();
  else if (problemName == "P4Test2DTet") return EllipticPDEProblemBuilder::EllipticPolynomialTest<4, TRIANGLE>().Build();
#endif
#if SpaceDimension >= 3
  else if (problemName == "LaplaceFicheraCube") return EllipticPDEProblemBuilder::LaplaceFicheraCube().Build();
  else if (problemName == "P0Test3D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<0, HEXAHEDRON>().Build();
  else if (problemName == "P1Test3D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<1, HEXAHEDRON>().Build();
  else if (problemName == "P2Test3D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<2, HEXAHEDRON>().Build();
  else if (problemName == "P3Test3D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<3, HEXAHEDRON>().Build();
  else if (problemName == "P4Test3D") return EllipticPDEProblemBuilder::EllipticPolynomialTest<4, HEXAHEDRON>().Build();
  else if (problemName == "P0Test3DTet") return EllipticPDEProblemBuilder::EllipticPolynomialTest<0, TETRAHEDRON>().Build();
  else if (problemName == "P1Test3DTet") return EllipticPDEProblemBuilder::EllipticPolynomialTest<1, TETRAHEDRON>().Build();
  else if (problemName == "P2Test3DTet") return EllipticPDEProblemBuilder::EllipticPolynomialTest<2, TETRAHEDRON>().Build();
  else if (problemName == "P3Test3DTet") return EllipticPDEProblemBuilder::EllipticPolynomialTest<3, TETRAHEDRON>().Build();
  else if (problemName == "P4Test3DTet") return EllipticPDEProblemBuilder::EllipticPolynomialTest<4, TETRAHEDRON>().Build();
#endif

  else Exit(problemName + " not found for SpaceDim=" + std::to_string(SpaceDimension))
}

std::unique_ptr<IEllipticProblem>
CreateEllipticProblemUnique(const std::string &problemName) {
  return CreateEllipticProblem(problemName);
}

std::shared_ptr<IEllipticProblem>
CreateEllipticProblemShared(const std::string &problemName) {
  return std::shared_ptr<IEllipticProblem>(CreateEllipticProblem(problemName));
}
