#include "MultiSampleFEM.hpp"

#include "CirculantEmbedding.hpp"
#include "EllipticPDESolver.hpp"
#include "IStochasticProblem.hpp" // Todo remove
#include "TransportPDESolver.hpp"

class StochasticLaplace1D : public MultiSampleFEM {
public:
  explicit StochasticLaplace1D(const MultiSampleFEMConfig &conf) :
      MultiSampleFEM("Interval", conf), pdeSolver(conf.pdeSolverConf) {}
protected:
  CirculantEmbedding1D circEmb;

  EllipticPDESolver pdeSolver;

  SampleSolution emptyProtocol(const SampleID &id) override { return {pdeSolver.GetSharedDisc(), id}; }

  SampleSolution protocol(const SampleID &id) override {
    auto y = circEmb.DrawSample(id, meshes);

    auto pdeProblem = EllipticPDEProblemBuilder() // Todo: GetLaplace1D(meshes)
      .WithMeshes(meshes)
      .WithPermeability([&y](const Point &x) -> Tensor { return y[floor(x[0] * y.size())] * One; })
      .WithSolution([](const Point &x) -> Scalar { return 1 - x[0]; })
      .WithFlux([](const Point &x) -> VectorField {
        return {-1.0, 0.0};
      })
      .WithLoad([](const Point &x, const cell &c) -> Scalar { return 0.0; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto ellipticSolution = pdeSolver.Run(pdeProblem, id.MeshIndex());

    return {ellipticSolution, id};
  }

  std::string Name() const override { return "StochasticLaplace1D"; };
};

class StochasticLaplace2D : public MultiSampleFEM {
public:
  explicit StochasticLaplace2D(const MultiSampleFEMConfig &conf) :
      MultiSampleFEM("Square", conf), pdeSolver(conf.pdeSolverConf) {}
protected:
  EllipticPDESolver pdeSolver;

  CirculantEmbedding2D circEmb;

  SampleSolution emptyProtocol(const SampleID &id) override { return {pdeSolver.GetSharedDisc(), id}; }

  SampleSolution protocol(const SampleID &id) override {
    auto y = circEmb.DrawSample(id,
                                meshes); // Todo finish refactorings CirculantEmbedding

    //    auto const logNormalDiffusionID = id.WithName("LogNormalDiffusion");
    //    // Todo idea: give in every step a name
    auto pdeProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability([&y](const Point &x) -> Tensor {
        return y[floor(x[0] * y.rows())][floor(x[1] * y.cols())] * One;
      })
      .WithSolution([](const Point &x) -> Scalar { return -x[1]; })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, -1.0};
      })
      .WithLoad([](const Point &x, const cell &c) -> Scalar { return 0.0; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto ellipticSolution = pdeSolver.Run(pdeProblem, id.MeshIndex());

    return {ellipticSolution, id};
  }

  std::string Name() const override { return "StochasticLaplace2D"; };
};

class StochasticLaplace3D : public MultiSampleFEM {
public:
  explicit StochasticLaplace3D(const MultiSampleFEMConfig &conf) :
      MultiSampleFEM("Square", conf), pdeSolver(conf.pdeSolverConf) {}
protected:
  CirculantEmbedding3D circEmb;

  EllipticPDESolver pdeSolver;

  SampleSolution emptyProtocol(const SampleID &id) override { return {pdeSolver.GetSharedDisc(), id}; }

  SampleSolution protocol(const SampleID &id) override {
    auto y = circEmb.DrawSample(id, meshes);

    auto pdeProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability([&y](const Point &x) -> Tensor {
        return y(floor(x[0] * y.FirstDimension()), floor(x[1] * y.SecondDimension()),
                  floor(x[2] * y.ThirdDimension()))
                * One;
      })
      .WithSolution([](const Point &x) -> Scalar { return -x[1]; })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, -1.0, 0.0};
      })
      .WithLoad([](const Point &x, const cell &c) -> Scalar { return 0.0; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto ellipticSolution = pdeSolver.Run(pdeProblem, id.MeshIndex());

    return {ellipticSolution, id};
  }

  std::string Name() const override { return "StochasticLaplace3D"; };
};

// Todo: Rewrite Problemclass:

class SCLaplace2D : public MultiSampleFEM {
public:
  explicit SCLaplace2D(const MultiSampleFEMConfig &conf) :
      MultiSampleFEM("Square", conf), pdeSolver(conf.pdeSolverConf) {}
protected:
  EllipticPDESolver pdeSolver;
  KLExpansionGenerator KLGenerator;

  SampleSolution emptyProtocol(const SampleID &id) override {
    return {pdeSolver.GetSharedDisc(), id};
  }

  SampleSolution protocol(const SampleID &id) override {
    auto y = KLGenerator.DrawSample(id, meshes);
    auto w = KLGenerator.GetSampleWeight();

    //    auto const logNormalDiffusionID = id.WithName("LogNormalDiffusion");

    auto pdeProblem = EllipticPDEProblemBuilder()
                          .WithMeshes(meshes)
                          .WithPermeability([&y](const Point &x) -> Tensor {
                            return exp(y[floor(x[0] * sqrt(y.size())) * floor(sqrt(y.size()) - 1)
                                         + floor(x[1] * sqrt(y.size()))])
                                   * One;
                          })
                          .WithSolution([](const Point &x) -> Scalar { return -x[1]; })
                          .WithFlux([](const Point &x) -> VectorField { return {0.0, -1.0}; })
                          .WithLoad([](const Point &x, const cell &c) -> Scalar { return 0.0; })
                          .WithName([&id]() { return id.IdString(); })
                          .BuildShared();

    auto ellipticSolution = pdeSolver.Run(pdeProblem, id.MeshIndex());

    return {ellipticSolution, id, w};
  }

  std::string Name() const override { return "SCLaplace2D"; };
};

class StochasticLaplace2DTest : public MultiSampleFEM {
public:
  explicit StochasticLaplace2DTest(const MultiSampleFEMConfig &conf) :
      MultiSampleFEM("Square", conf), pdeSolver(conf.pdeSolverConf) {}
protected:
  CirculantEmbedding2D circEmb;

  EllipticPDESolver pdeSolver;

  SampleSolution emptyProtocol(const SampleID &id) override { return {pdeSolver.GetSharedDisc(), id}; }

  SampleSolution protocol(const SampleID &id) override {
    auto y = Random::Uniform(id.commSplit, 2, -1.0, 1.0);

    auto pdeProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability([&y](const Point &x) -> Tensor {
        if (x[1] > 0.125 && x[1] < 0.375) {
          if (x[0] > 0.125 && x[0] < 0.375) return (1 + 0.9 * y[0]) * One;
          if (x[0] > 0.625 && x[0] < 0.875) return (1 + 0.9 * y[1]) * One;
        }
        return One;
      })
      .WithSolution([](const Point &x) -> Scalar { return -x[1]; })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, -1.0};
      })
      .WithLoad([](const Point &x, const cell &c) -> Scalar { return 0.0; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto ellipticSolution = pdeSolver.Run(pdeProblem, id.MeshIndex());

    return {ellipticSolution, id};
  }

  std::string Name() const override { return "StochasticLaplace2DTest"; };
};

class UniformDistributionLaplace2D : public MultiSampleFEM {
public:
  explicit UniformDistributionLaplace2D(const MultiSampleFEMConfig &conf) :
      MultiSampleFEM("Square", conf), pdeSolver(conf.pdeSolverConf) {}
protected:
  EllipticPDESolver pdeSolver;

  SampleSolution emptyProtocol(const SampleID &id) override { return {pdeSolver.GetSharedDisc(), id}; }

  SampleSolution protocol(const SampleID &id) override {
    auto y = Random::Uniform(id.commSplit, 4, -sqrt(12) / 2.0, sqrt(12) / 2.0);

    auto pdeProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability([&y](const Point &x) -> Tensor {
        return One
                * (0.01
                  + exp(exp(-1.0 / 8.0)
                        * (cos(2 * M_PI * x[0]) * y[0] + sin(2 * M_PI * x[0]) * y[1]
                            + cos(2 * M_PI * x[1]) * y[2] + sin(2 * M_PI * x[1]) * y[3])));
      })
      .WithSolution([](const Point &x) -> Scalar { return -x[1]; })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, -1.0};
      })
      .WithLoad([](const Point &x, const cell &c) -> Scalar { return 0.0; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto ellipticSolution = pdeSolver.Run(pdeProblem, id.MeshIndex());

    return {ellipticSolution, id};
  }

  string Name() const { return "UniformDistributionLaplace2D"; }
};

class OCStochasticLaplace2D : public MultiSampleFEM {
public:
  double targetCost = 0.0;

  explicit OCStochasticLaplace2D(const MultiSampleFEMConfig &conf) :
      MultiSampleFEM("Square", conf), pdeSolver(conf.pdeSolverConf) {
    Config::Get("targetCost", targetCost);
  }
protected:
  CirculantEmbedding2D circEmb;

  EllipticPDESolver pdeSolver;

  SampleSolution emptyProtocol(const SampleID &id) override { return {pdeSolver.GetSharedDisc(), id}; }

  SampleSolution protocol(const SampleID &id,
                          const SampleSolution &control) override { // Todo: Additional bool for calculating only the
                                                                    // state solution.
    auto y = circEmb.DrawSample(id, meshes); // Todo finish
    // control.solution.values = this->ComputeValues(control.solution.vector)

    auto statePDEProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability(
          [&y](const Point &x) -> Tensor { return y[floor(x[0] * y.rows())][floor(x[1] * y.cols())] * One; })
      .WithSolution([](const Point &x) -> Scalar { return sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]); })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, -1.0};
      })
      .WithLoad([&control](const Point &x, const cell &c) -> Scalar {
        ScalarElement elem(control.solution.vector, *c);
        return elem.Value(c.GlobalToLocal(x), control.solution.vector);
      })
      .WithHasExactSolution([]() -> bool { return true; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto stateSolution = pdeSolver.Run(statePDEProblem, id.MeshIndex());

    auto adjointPDEProblem = EllipticPDEProblemBuilder() // Todo fix me
      .WithMeshes(meshes)
      .WithPermeability([&y](const Point &x) -> Tensor {
        return y[floor(x[0] * y.rows())][floor(x[1] * y.cols())] * One;
      })
      .WithSolution([](const Point &x) -> Scalar { return 0.0; })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, 0.0};
      })
      .WithLoad([&stateSolution](const Point &x, const cell &c) -> Scalar {
        ScalarElement elem(stateSolution.vector, *c);
        return elem.Value(c.GlobalToLocal(x), stateSolution.vector)
              - sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
      })
      .WithHasExactSolution([]() -> bool { return true; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto adjointSolution = pdeSolver.Run(adjointPDEProblem, id.MeshIndex());

    adjointSolution.vector += targetCost * control.solution.vector; // Todo: Check sgn

    return {adjointSolution, id}; // Todo: Also give back the control for the
                                  // override of QoI valueMap.
  }

  // Todo: Need an option to estimate the function value of the control problem
  // somewhere here
  //  j(u) = E[1/2 * |solution - state(\xi)|²_L²] + targetCost * |control|²_L²

  ValueMap ComputeValues(const SampleSolution &control) { return pdeSolver.ComputeValues(control.solution.vector); }

  SampleSolution protocol(const SampleID &id) override {
    auto y = circEmb.DrawSample(id, meshes); // Todo finish

    auto statePDEProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability(
          [&y](const Point &x) -> Tensor { return y[floor(x[0] * y.rows())][floor(x[1] * y.cols())] * One; })
      .WithSolution([](const Point &x) -> Scalar { return sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]); })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, -1.0};
      })
      .WithLoad([](const Point &x, const cell &c) -> Scalar { return 0.0; })
      .WithHasExactSolution([]() -> bool { return true; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto stateSolution = pdeSolver.Run(statePDEProblem, id.MeshIndex());

    auto adjointPDEProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability([&y](const Point &x) -> Tensor {
        return y[floor(x[0] * y.rows())][floor(x[1] * y.cols())] * One;
      })
      .WithSolution([](const Point &x) -> Scalar { return 0.0; })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, 0.0};
      })
      .WithLoad([&stateSolution](const Point &x, const cell &c) -> Scalar {
        ScalarElement elem(stateSolution.vector, *c);
        return elem.Value(c.GlobalToLocal(x), stateSolution.vector)
              - sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
      })
      .WithHasExactSolution([]() -> bool { return true; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto adjointSolution = pdeSolver.Run(adjointPDEProblem, id.MeshIndex());

    return {adjointSolution, id};
  }

  double protocolQOI(const SampleID &id, const SampleSolution &control) override {
    auto y = circEmb.DrawSample(id, meshes); // Todo: Need to calculate L² Norm of control (once at some point)

    auto statePDEProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability(
          [&y](const Point &x) -> Tensor { return y[floor(x[0] * y.rows())][floor(x[1] * y.cols())] * One; })
      .WithSolution([](const Point &x) -> Scalar { return sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]); })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, -1.0};
      })
      .WithLoad([&control](const Point &x, const cell &c) -> Scalar {
        ScalarElement elem(control.solution.vector, *c);
        return elem.Value(c.GlobalToLocal(x), control.solution.vector);
      })
      .WithHasExactSolution([]() -> bool { return true; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto stateSolution = pdeSolver.Run(statePDEProblem, id.MeshIndex());

    return stateSolution.values["L2Error"]; // give back expectation value from norm and norm from expextation field
                                            // (interesting to look on)
  }

  std::string Name() const override { return "OCStochasticLaplace2D"; }

  //  Scalar Load(const Point &x, const cell &c) const override { // Todo
  //  include me
  //    ScalarElement elem(*rhs, *c);
  //    if (!with_data) {
  //      return elem.Value(c.GlobalToLocal(x), *rhs);
  //    } else {
  //      return Solution(x) - elem.Value(c.GlobalToLocal(x), *rhs);
  //    }
  //  }
};

class OCLaplace2D : public MultiSampleFEM {
public:
  double targetCost = 0.0;

  explicit OCLaplace2D(const MultiSampleFEMConfig &conf) :
      MultiSampleFEM("Square", conf), pdeSolver(conf.pdeSolverConf) {
    Config::Get("targetCost", targetCost);
  }
protected:
  EllipticPDESolver pdeSolver;

  SampleSolution emptyProtocol(const SampleID &id) override { return {pdeSolver.GetSharedDisc(), id}; }

  SampleSolution protocol(const SampleID &id, const SampleSolution &control) override {

    auto statePDEProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability([](const Point &x) -> Tensor { return One; })
      .WithSolution([](const Point &x) -> Scalar { return sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]); })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, -1.0};
      })
      .WithLoad([&control](const Point &x, const cell &c) -> Scalar {
        ScalarElement elem(control.solution.vector, *c);
        return elem.Value(c.GlobalToLocal(x), control.solution.vector);
      })
      .WithHasExactSolution([]() -> bool { return true; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto stateSolution = pdeSolver.Run(statePDEProblem, id.MeshIndex());

    auto adjointPDEProblem = EllipticPDEProblemBuilder() // Todo fix me
      .WithMeshes(meshes)
      .WithPermeability([](const Point &x) -> Tensor { return One; })
      .WithSolution([](const Point &x) -> Scalar { return 0.0; })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, 0.0};
      })
      .WithLoad([&stateSolution](const Point &x, const cell &c) -> Scalar {
        ScalarElement elem(stateSolution.vector, *c);
        return elem.Value(c.GlobalToLocal(x), stateSolution.vector)
              - sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
      })
      .WithHasExactSolution([]() -> bool { return true; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto adjointSolution = pdeSolver.Run(adjointPDEProblem, id.MeshIndex());

    adjointSolution.vector += targetCost * control.solution.vector; // Todo: Check sgn

    return {adjointSolution, id}; // Todo: Also give back the control for the
                                  // override of QoI valueMap.
  }

  // Todo: Need an option to estimate the function value of the control problem
  // somewhere here
  //  j(u) = E[1/2 * |solution - state(\xi)|²_L²] + targetCost * |control|²_L²

  ValueMap ComputeValues(const SampleSolution &control) { return pdeSolver.ComputeValues(control.solution.vector); }

  SampleSolution protocol(const SampleID &id) override {
    auto statePDEProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability([](const Point &x) -> Tensor { return One; })
      .WithSolution([](const Point &x) -> Scalar { return sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]); })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, -1.0};
      })
      .WithLoad([](const Point &x, const cell &c) -> Scalar { return 0.0; })
      .WithHasExactSolution([]() -> bool { return true; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto stateSolution = pdeSolver.Run(statePDEProblem, id.MeshIndex());

    auto adjointPDEProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability([](const Point &x) -> Tensor { return One; })
      .WithSolution([](const Point &x) -> Scalar { return 0.0; })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, 0.0};
      })
      .WithLoad([&stateSolution](const Point &x, const cell &c) -> Scalar {
        ScalarElement elem(stateSolution.vector, *c);
        return elem.Value(c.GlobalToLocal(x), stateSolution.vector)
              - sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
      })
      .WithHasExactSolution([]() -> bool { return true; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto adjointSolution = pdeSolver.Run(adjointPDEProblem, id.MeshIndex());

    return {adjointSolution, id};
  }

  double protocolQOI(const SampleID &id, const SampleSolution &control) override {
    auto statePDEProblem = EllipticPDEProblemBuilder()
      .WithMeshes(meshes)
      .WithPermeability([](const Point &x) -> Tensor { return One; })
      .WithSolution([](const Point &x) -> Scalar { return sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]); })
      .WithFlux([](const Point &x) -> VectorField {
        return {0.0, -1.0};
      })
      .WithLoad([&control](const Point &x, const cell &c) -> Scalar {
        ScalarElement elem(control.solution.vector, *c);
        return elem.Value(c.GlobalToLocal(x), control.solution.vector);
      })
      .WithHasExactSolution([]() -> bool { return true; })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    auto stateSolution = pdeSolver.Run(statePDEProblem, id.MeshIndex());

    return std::pow(stateSolution.values["L2Error"], 2)
           / 2; // give back expectation value from norm and norm from expextation field (interesting to look on)
  }

  std::string Name() const override { return "OCLaplace2D"; }
};

// class OCStochasticLaplace3D : public StochasticLaplace3D { // Todo include me
// like above protected:
//   std::shared_ptr<Vector> rhs;
//
// public:
//   OCStochasticLaplace3D() : IProblem("Hexahedron") {}
//
//  bool with_data = false;
//
//  Scalar Solution(const Point &x) const override {
//    return sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) *
//           sin(2.0 * M_PI * x[2]);
//  }
//
//  void LoadRHS(const Vector &u, bool withData) override {
//    with_data = withData;
//    rhs = std::make_shared<Vector>(u);
//  }
//
//  Scalar Load(const Point &x, const cell &c) const override {
//    ScalarElement elem(*rhs, *c);
//    if (!with_data) {
//      return elem.Value(c.GlobalToLocal(x), *rhs);
//    } else {
//      return Solution(x) - elem.Value(c.GlobalToLocal(x), *rhs);
//    }
//  }
//};
//
// class OptimalLaplace2D : public StochasticLaplace2D {
// protected:
//  std::shared_ptr<Vector> rhs;
//
// public:
//  OptimalLaplace2D() : IProblem("Square") {}
//
//  bool with_data = false;
//
//  bool HasExactSolution() const override { return true; }
//
//  // Scalar Load(const Point &x, const cell &c) const override { return 0.0; }
//
//  Scalar Solution(const Point &x) const override {
//    return sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
//  }  // nontrivial example
//
//  VectorField Flux(const Point &x) const override { return {0.0, -1.0}; }
//
//  Tensor Permeability(const Point &x) const override { return One; }
//
//  std::string Name() const override { return "Laplace"; }
//
//  void LoadRHS(const Vector &u, bool withData) override {
//    with_data = withData;
//    rhs = std::make_shared<Vector>(u);
//  }
//
//  Scalar Load(const Point &x, const cell &c) const override {
//    ScalarElement elem(*rhs, *c);
//    if (!with_data) {
//      // mout << with_data << x << elem.Value(x, *rhs) << endl;
//      return elem.Value(c.GlobalToLocal(x), *rhs);
//    } else {
//      //      mout << with_data << x << Solution(x) - elem.Value(x, *rhs) <<
//      //      endl;
//      return Solution(x) - elem.Value(c.GlobalToLocal(x), *rhs);
//    }
//  }
//};

#if 0

class SparseGridLaplace2D : public IStochasticEllipticProblem {
private:
  double factor = exp(-1.0 / 8.0);

  double a = -sqrt(12) / 2.0;

  double b = sqrt(12) / 2.0;

  double a_min = 0.01;

  SparseGridGenerator sparseGridGen;

  void drawSample(const SampleID &id) override { sparseGridGen.DrawSample(id); }
public:
  SparseGridLaplace2D() :
      IProblem("Square"), sparseGridGen(SparseGridGenerator(GridDomain({{a, a, a, a}, {b, b, b, b}}))) {}

  void InitGenerator(int init) override { sparseGridGen.InitGenerator(meshes, init); }

  // Todo why is weight not working
  //  double SumOfWeights() override {
  //    return sparseGridGen.SumOfWeights();
  //  }
  //
  //  double SampleWeight() const override {
  //    return sparseGridGen.SampleWeight(currentID);
  //  }

  int NumOfSamples() const override { return sparseGridGen.GetNumPoints(); }

  Scalar Load(const Point &x, const cell &c) const override { return 0.0; }

  Scalar Solution(const Point &x) const override { return -x[1]; }

  VectorField Flux(const Point &x) const override { return {0.0, -1.0}; }

  Tensor Permeability(const Point &x) const override {
    return (a_min
            + exp(factor
                  * (cos(2 * M_PI * x[0]) * this->sparseGridGen.EvalSample()[0]
                     + sin(2 * M_PI * x[0]) * this->sparseGridGen.EvalSample()[1]
                     + cos(2 * M_PI * x[1]) * this->sparseGridGen.EvalSample()[2]
                     + sin(2 * M_PI * x[1]) * this->sparseGridGen.EvalSample()[3])))
           * One;
  }

  string Name() const override { return "SparseGridLaplace2D"; }
};

#endif

/*
 * Only Test Protocols
 */

class EmptyLagrangeProtocol : public MultiSampleFEM {
public:
  explicit EmptyLagrangeProtocol(const MultiSampleFEMConfig &conf) :
      MultiSampleFEM("Square", conf),
      disc(std::make_shared<const LagrangeDiscretization>(*meshes, conf.pdeSolverConf.degree)) {}
protected:
  std::shared_ptr<const LagrangeDiscretization> disc;

  SampleSolution emptyProtocol(const SampleID &id) override { return {disc, id}; }

  SampleSolution protocol(const SampleID &id) override { return emptyProtocol(id); }

  double selectQoI(const SampleSolution &sample) const override {
    if (conf.quantityOfInterest == "ProcEvaluation") { return PPM->Proc(0); }
    if (conf.quantityOfInterest == "ZeroEvaluation") { return 0.0; }
    if (conf.quantityOfInterest == "ModuloIndexEvaluation") { return (sample.id.number % 2) ? 1.0 : -1.0; }
    Exit("Quantity of Interest: " + conf.quantityOfInterest + " not implemented")
  }

  std::string Name() const override { return "EmptyLagrangeProtocol"; }
};

class NormalLagrangeProtocol : public MultiSampleFEM {
public:
  explicit NormalLagrangeProtocol(const MultiSampleFEMConfig &conf) :
      MultiSampleFEM("Square", conf),
      disc(std::make_unique<LagrangeDiscretization>(*meshes, conf.pdeSolverConf.degree)) {}
protected:
  std::shared_ptr<LagrangeDiscretization> disc;

  SampleSolution emptyProtocol(const SampleID &id) override { return {disc, id}; }

  SampleSolution protocol(const SampleID &id) override {
    auto y = Random::Normal(id.commSplit);

    auto sample = SampleSolution(disc, id);
    if (id.coarse) {
      sample.solution.values[conf.quantityOfInterest] = 0.0;
      sample.solution.vector = 0.0;
    } else {
      sample.solution.values[conf.quantityOfInterest] = y;

      // Todo, think about how this could be solved nicer
      sample.solution.vector = y;
      for (row r = sample.solution.vector.rows(); r != sample.solution.vector.rows_end(); r++) {
        sample.solution.vector(r, 0) = -r()[1] + y;
      }
    }

    return sample;
  }

  std::string Name() const override { return "NormalLagrangeProtocol"; }
};

class StochasticTransportTodo : public MultiSampleFEM {
public:
  // Todo fix constructor here -> Maybe use list of PDESolverConfigs in MSFEMConfig?
  explicit StochasticTransportTodo(const MultiSampleFEMConfig &conf) :
      MultiSampleFEM("Square", conf), ellipticPDESolver(PDESolverConfig()), transportPDESolver(PDESolverConfig()) {}
protected:
  EllipticPDESolver ellipticPDESolver;

  TransportPDESolver transportPDESolver;

  SampleSolution emptyProtocol(const SampleID &id) override { return {transportPDESolver.GetSharedDisc(), id}; }

  SampleSolution protocol(const SampleID &id) override {

    // Todo different elliptic problems can be setup here (1d, 2d, log-normal, etc.),
    //  this should result in different protocols with meaningful names.
    RVector sample = Random::Uniform(id.commSplit, 2, -1.0, 1.0);
    auto stochasticLaplace = EllipticPDEProblemBuilder() // Todo GetRawLaplace(meshes)
      .WithMeshes(meshes)
      .WithPermeability([&sample](const Point &x) -> Tensor {
        if (x[1] > 0.125 && x[1] < 0.375) {
          if (x[0] > 0.125 && x[0] < 0.375) return (1 + 0.9 * sample[0]) * One;
          if (x[0] > 0.625 && x[0] < 0.875) return (1 + 0.9 * sample[1]) * One;
        }
        return One;
      })
      .WithName([&id]() { return id.IdString(); })
      .BuildShared();

    // Todo New struct EllipticSolution struct with fluxVector
    //  (in elliptic PDESolver SetNormalFlux has to be called)
    auto ellipticSolution = ellipticPDESolver.Run(stochasticLaplace);

    // Todo write builder pattern for TransportPDESolver
    //    auto stochasticPollution = TransportPDEProblemBuilder("Square")
    //        .WithFlux([&ellipticSolution](const Cell &c, const Point &x) -> VectorField {
    //          RTLagrangeElementT<> elem(ellipticSolution.fluxVector, c);
    //          return elem.CellFlux(ellipticSolution.fluxVector); // Todo only use this one
    //        })
    //        .WithFaceNormalFlux([&ellipticSolution](int face, const Cell &c) -> Scalar {
    //          RTLagrangeElement elem(ellipticSolution.vector, c);
    //          RTFaceElement faceElem(ellipticSolution.vector, c, face);
    //          return (ellipticSolution.vector)(elem[face], 0) * elem.Sign(face) / faceElem.Area();
    //        })
    //        .WithName([&id]() { return id.IdString(); })
    //        .BuildShared();
    //
    auto transportSolution = transportPDESolver.Run( // Todo put in problem of above
        CreateTransportProblemShared(conf.pdeSolverConf.problemName));
    return {transportSolution, id};
  }

  std::string Name() const override { return "StochasticTransportTodo"; }
};

MultiSampleFEM *CreateMSFEM(const MultiSampleFEMConfig &conf) {

  if (conf.protocolName == "StochasticLaplace1D") return new StochasticLaplace1D(conf);
  if (conf.protocolName == "StochasticLaplace2D") return new StochasticLaplace2D(conf);
  if (conf.protocolName == "StochasticLaplace3D") return new StochasticLaplace3D(conf);
  if (conf.protocolName == "StochasticLaplace2DTest") return new StochasticLaplace2DTest(conf);
  if (conf.protocolName == "SCLaplace2D") return new SCLaplace2D(conf);

  if (conf.protocolName == "UniformDistributionLaplace2D") return new UniformDistributionLaplace2D(conf);

  if (conf.protocolName == "OCStochasticLaplace2D") return new OCStochasticLaplace2D(conf);
  //  if (conf.protocolName == "OCStochasticLaplace3D")
  //    return new OCStochasticLaplace3D(conf);
  //
  if (conf.protocolName == "OCLaplace2D") return new OCLaplace2D(conf);

  if (conf.protocolName == "EmptyLagrangeProtocol") return new EmptyLagrangeProtocol(conf);

  if (conf.protocolName == "NormalLagrangeProtocol") return new NormalLagrangeProtocol(conf);

  if (conf.protocolName == "StochasticTransport123") return new StochasticTransportTodo(conf);

#if 0
  if (conf.protocolName == "SparseGridLaplace2D") return new SparseGridLaplace2D(conf);
#endif
  Exit(conf.protocolName + " not found");
  return nullptr;
}

std::unique_ptr<MultiSampleFEM> CreateUniqueMSFEM(const MultiSampleFEMConfig &conf) {
  return std::unique_ptr<MultiSampleFEM>(CreateMSFEM(conf));
}

std::shared_ptr<MultiSampleFEM> CreateSharedMSFEM(const MultiSampleFEMConfig &conf) {
  return std::shared_ptr<MultiSampleFEM>(CreateMSFEM(conf));
}
