#ifndef TUTORIAL_ELLIPTICPROBLEM_HPP
#define TUTORIAL_ELLIPTICPROBLEM_HPP

#include <utility>

#include "IProblem.hpp"
#include "Algebra.hpp"


struct EllipticOperands {
  using LoadValue = Scalar;
  using SolutionValue = Scalar;
  using FluxValue = VectorField;
  using PermeabilityValue = Tensor;

  std::string_view name{ "Elliptic" };
};

struct VectorValuedOperands {
  using LoadValue = VectorField;
  using SolutionValue = VectorField;
  using FluxValue = Tensor;
  using PermeabilityValue = Tensor;

  std::string_view name{ "VectorValued" };
};

template<class T>
concept EllipticOperand = std::is_same<T, EllipticOperands>::value || std::is_same<T, VectorValuedOperands>::value;

template<EllipticOperand Operands>
class GenericEllipticProblem : virtual public IProblem {
public:
  virtual Operands::LoadValue Load(const Point &x, const cell &c) const = 0;

  virtual Operands::SolutionValue Solution(const Point &x) const = 0;

  virtual Operands::FluxValue Flux(const Point &x) const = 0;

  virtual Operands::PermeabilityValue Permeability(const Point &x) const = 0;

  void Solution(Vector &u) const {
    for (row r = u.rows(); r != u.rows_end(); ++r) {
      auto S = Solution(r());
      for (int d = 0; d < u.dim(); ++d) {
        u(r, d) = S[d];
      }
    }
  }

  void Permeability(Vector &kappa) const {
    for (row r = kappa.rows(); r != kappa.rows_end(); ++r) {
      auto P = Permeability(r());
      for (int k = 0; k < kappa.dim(); ++k) { // TODO anpassen
        for (int d = 0; d < kappa.dim(); ++d) {
          // Adding everything up to scalar value, since vector is of size = 1
          kappa(r(), 0) += P[k][d] / kappa.dim();
        }
      }
    }
  }
};

using IEllipticProblem = GenericEllipticProblem<EllipticOperands>;
using IVectorValuedProblem = GenericEllipticProblem<VectorValuedOperands>;

class EllipticPDEProblem : public IEllipticProblem {
private:
  friend class EllipticPDEProblemBuilder;

  using FluxFunction = std::function<VectorField(const Point &x)>;
  using PermeabilityFunction = std::function<Tensor(const Point &x)>;
  using SolutionFunction = std::function<Scalar(const Point &x)>;
  using LoadFunction = std::function<Scalar(const Point &x, const cell &c)>;
  using NameFunction = std::function<std::string()>;
  using HasExactSolutionFunction = std::function<bool()>;

  FluxFunction fluxFunction;
  PermeabilityFunction permeabilityFunction;
  SolutionFunction solutionFunction;
  LoadFunction loadFunction;
  NameFunction nameFunction;
  HasExactSolutionFunction hasExactSolution;

public:
  explicit EllipticPDEProblem(std::shared_ptr<Meshes> meshes) : IProblem(std::move(meshes)) {}

  explicit EllipticPDEProblem(std::string meshName) : IProblem(std::move(meshName)) {}

  explicit EllipticPDEProblem() : IProblem() {}

  EllipticPDEProblem(EllipticPDEProblem&);

  Tensor Permeability(const Point &x) const override {
    if (permeabilityFunction) return permeabilityFunction(x);
    return One;
  }

  Scalar Solution(const Point &x) const override {
    if (solutionFunction) return solutionFunction(x);
    return 0.0;
  }

  VectorField Flux(const Point &x) const override {
    if (fluxFunction) return fluxFunction(x);
    return 0.0;
  }

  Scalar Load(const Point &x, const cell &c) const override {
    if (loadFunction) return loadFunction(x, c);
    return 0.0;
  }

  bool HasExactSolution() const override {
      if (hasExactSolution) return hasExactSolution();
      return false;
  }

  std::string Name() const override {
    if (nameFunction) return nameFunction();
    return "EllipticPDEProblem";
  }
};

class EllipticPDEProblemBuilder {
private:
  std::string meshName = "";
  std::shared_ptr<Meshes> meshes;
  EllipticPDEProblem::FluxFunction fluxFunction;
  EllipticPDEProblem::PermeabilityFunction permeabilityFunction;
  EllipticPDEProblem::SolutionFunction solutionFunction;
  EllipticPDEProblem::LoadFunction loadFunction;
  EllipticPDEProblem::NameFunction nameFunction;
  EllipticPDEProblem::HasExactSolutionFunction hasExactSolution;

public:
  [[nodiscard]] explicit EllipticPDEProblemBuilder() {}

  [[nodiscard]] EllipticPDEProblemBuilder &
  WithPermeability(std::function<Tensor(const Point &x)> func) {
    permeabilityFunction = std::move(func);
    return *this;
  }

  [[nodiscard]] EllipticPDEProblemBuilder &
  WithSolution(EllipticPDEProblem::SolutionFunction func) {
    solutionFunction = std::move(func);
    return *this;
  }

  [[nodiscard]] EllipticPDEProblemBuilder &
  WithFlux(EllipticPDEProblem::FluxFunction func) {
    fluxFunction = std::move(func);
    return *this;
  }

  [[nodiscard]] EllipticPDEProblemBuilder &
  WithLoad(EllipticPDEProblem::LoadFunction func) {
    loadFunction = std::move(func);
    return *this;
  }

  [[nodiscard]] EllipticPDEProblemBuilder &
  WithHasExactSolution(EllipticPDEProblem::HasExactSolutionFunction func) {
      hasExactSolution = std::move(func);
      return *this;
  }

  [[nodiscard]] EllipticPDEProblemBuilder &
  WithName(EllipticPDEProblem::NameFunction func) {
    nameFunction = std::move(func);
    return *this;
  }

  [[nodiscard]] EllipticPDEProblemBuilder &
  WithMeshName(std::string&& name) {
    meshName = name;
    return *this;
  }

  /*
  * Set a problem mesh
  *
  * Calling this overrides any mesh name that was set before,
  * causing it to be ignored.
  */
  [[nodiscard]] EllipticPDEProblemBuilder &
  WithMeshes(std::shared_ptr<Meshes> sharedMeshes) {
    meshes = sharedMeshes;
    return *this;
  }

  [[nodiscard]] std::unique_ptr<EllipticPDEProblem> Build() {
    std::unique_ptr<EllipticPDEProblem> problem;
    if (meshes) {
      problem = std::make_unique<EllipticPDEProblem>(meshes);
    } else {
      problem = std::make_unique<EllipticPDEProblem>(meshName);
    }

    problem->hasExactSolution = std::move(hasExactSolution);
    problem->nameFunction = std::move(nameFunction);
    problem->loadFunction = std::move(loadFunction);
    problem->solutionFunction = std::move(solutionFunction);
    problem->fluxFunction = std::move(fluxFunction);
    problem->permeabilityFunction = std::move(permeabilityFunction);
    return problem;
  }

  [[nodiscard]] std::shared_ptr<EllipticPDEProblem> BuildShared() {
    return std::shared_ptr<EllipticPDEProblem>(Build());
  }

  [[nodiscard]] static EllipticPDEProblemBuilder Laplace();
  [[nodiscard]] static EllipticPDEProblemBuilder Laplace1D();
  [[nodiscard]] static EllipticPDEProblemBuilder Laplace2D();
  [[nodiscard]] static EllipticPDEProblemBuilder LaplaceSquare500();
  [[nodiscard]] static EllipticPDEProblemBuilder LaplaceSquare501();
  [[nodiscard]] static EllipticPDEProblemBuilder LaplaceSquare2();
  [[nodiscard]] static EllipticPDEProblemBuilder Discontinuous();
  [[nodiscard]] static EllipticPDEProblemBuilder Discontinuous1D();
  [[nodiscard]] static EllipticPDEProblemBuilder Discontinuous2D();
  [[nodiscard]] static EllipticPDEProblemBuilder DiscontinuousSquare500();
  [[nodiscard]] static EllipticPDEProblemBuilder Kellogg();
  [[nodiscard]] static EllipticPDEProblemBuilder Rock();
  [[nodiscard]] static EllipticPDEProblemBuilder Divergent();
  [[nodiscard]] static EllipticPDEProblemBuilder LaplaceFicheraCube();
  [[nodiscard]] static EllipticPDEProblemBuilder LinearAffineInflow();
  
  template<int degree, CELLTYPE type>
  [[nodiscard]] static EllipticPDEProblemBuilder EllipticPolynomialTest() {
    return EllipticPDEProblemBuilder()
      .WithMeshName(to_string(type))
      .WithHasExactSolution([]() { return true; })
      .WithName([]() { 
        return "Elliptic-P" + std::to_string(degree) + "Test" + to_string(type);
      })
      .WithLoad([](const Point &x, const cell &c) {
        return -CellDim(type) * degree * (degree - 1) * pow(x.Sum(), degree - 2);
      })
      .WithSolution([](const Point &x) {
        return pow(x.Sum(), degree);
      })
      .WithFlux([](const Point &x) -> VectorField {
        VectorField F;
        for (int i = 0; i < CellDim(type); ++i) {
          F[i] = degree * pow(x.Sum(), degree - 1);
        }
        return F;
      })
      .WithPermeability([](const Point &x) {
        return One;
      });
  }
};

std::unique_ptr<IEllipticProblem> CreateEllipticProblem(const string &problemName);

[[deprecated("use CreateEllipticProblem instead")]]
std::unique_ptr<IEllipticProblem>
CreateEllipticProblemUnique(const std::string &problemName);

std::shared_ptr<IEllipticProblem>
CreateEllipticProblemShared(const std::string &problemName);


template<int degree, CELLTYPE type>
class VectorValuedPolynomialTest : public IVectorValuedProblem {
public:
  bool HasExactSolution() const override { return true; }

  string Name() const override {
    return "VectorValued-P" + std::to_string(degree) + "Test" + to_string(type);
  }

  VectorField Load(const Point &x, const cell &c) const override {
    double val = -CellDim(type) * degree * (degree - 1) * pow(x.Sum(), degree - 2);
    return {val, CellDim(type) > 1 ? val : 0.0, CellDim(type) > 2 ? val : 0.0};
  }

  VectorField Solution(const Point &x) const override {
    return {pow(x.Sum(), degree), CellDim(type) > 1 ? pow(x.Sum(), degree) : 0.0,
            CellDim(type) > 2 ? pow(x.Sum(), degree) : 0.0};
  }

  Tensor Flux(const Point &x) const override {
    double val = degree * pow(x.Sum(), degree - 1);
    int cell_dim = CellDim(type);
    double a = CellDim(type) > 1 ? val: 0.0;
    double b = CellDim(type) > 2 ? val : 0.0;
    return Tensor(
        val, a, b,
        a, a, b,
        b, b, b);
  }

  Tensor Permeability(const Point &x) const override {
    return One;
  }

  VectorValuedPolynomialTest() : IProblem(to_string(type)) {}
};


IVectorValuedProblem *CreateVectorValuedProblem(const string &problemName);

std::unique_ptr<IVectorValuedProblem>
CreateVectorValuedProblemUnique(const std::string &problemName);

std::shared_ptr<IVectorValuedProblem>
CreateVectorValuedProblemShared(const std::string &problemName);


#endif //TUTORIAL_ELLIPTICPROBLEM_HPP
