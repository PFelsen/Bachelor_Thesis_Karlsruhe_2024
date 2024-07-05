#ifndef REACTIONPROBLEMS_HPP
#define REACTIONPROBLEMS_HPP

#include "IProblem.hpp"
#include "Algebra.hpp"


class IReactionProblem : virtual public IProblem {
protected:
  double diffusion = 0.01;

  double convection = 1.0;

  double reaction = 5.0;

  double CFL = 0.25;

  double t0 = 0.0;

  double dt = 0.0;

  double T = 1.0;

public:
  IReactionProblem() {
    Config::Get("CFL", CFL);
    Config::Get("t0", t0);
    Config::Get("T", T);
    Config::Get("Diffusion", diffusion);
    Config::Get("Convection", convection);
    Config::Get("Reaction", reaction);
  }

  double GetConvection() const { return convection; }

  double GetDiffusion() const { return diffusion; }

  double GetReaction() const { return reaction; }

  virtual Scalar Load(const Point &x) const = 0;

  virtual VectorField Flux(const Point &x) const = 0;

  virtual Tensor Diffusion(const Point &x) const = 0;

  virtual VectorField Convection(const Point &x) const = 0;

  virtual Scalar Reaction(const Point &x, double c) const = 0;

  virtual Scalar Concentration(double t, const Point &x) const = 0;

  virtual Scalar DerivativeReaction(const Point &, double c) const = 0;

  virtual VectorField CellConvection(const Cell &c, const Point &x) const {
    return Convection(x);
  }

  virtual double FaceConvection(const Cell &c, int face,
                                const VectorField &N, const Point &x) const {
    return Convection(x) * N;
  }


  double GetStepSize(double h) const {
    if (dt == 0.0) return CFL * h;
    else return dt;
  }

  double GetStartTime() const { return t0; }

  double GetEndTime() const { return T; }

  double GetCFL() const { return CFL; }

  bool RHS() const { return true; }
};

class ReactionPDEProblem: public IReactionProblem {
  friend class ReactionPDEProblemBuilder;

  using LoadFunction = std::function<Scalar(const Point &x)>;
  using FluxFunction = std::function<VectorField(const Point &x)>;
  using DiffusionFunction = std::function<Tensor(const Point &x)>;
  using ConvectionFunction = std::function<VectorField(const Point &x)>;
  using ReactionFunction = std::function<Scalar(const Point &x, double c)>;
  using ConcentrationFunction = std::function<Scalar(double t, const Point &x)>;
  using DerivativeReactionFunction = std::function<Scalar(const Point &x, double c)>;
  using CellConvectionFunction = std::function<VectorField(const Cell &c, const Point &x)>;
  using FaceConvectionFunction = std::function<double(const Cell &c, int face, const VectorField &N, const Point &x)>;
  using NameFunction = std::function<std::string()>;

private:
  LoadFunction loadFunction = nullptr;
  FluxFunction fluxFunction = nullptr;
  DiffusionFunction diffusionFunction = nullptr;
  ConvectionFunction convectionFunction = nullptr;
  ReactionFunction reactionFunction = nullptr;
  ConcentrationFunction concentrationFunction = nullptr;
  DerivativeReactionFunction derivativeReactionFunction = nullptr;
  CellConvectionFunction cellConvectionFunction = nullptr;
  FaceConvectionFunction faceConvectionFunction = nullptr;
  NameFunction nameFunction = nullptr;


public:
  explicit ReactionPDEProblem(std::shared_ptr<Meshes> meshes) : IProblem(std::move(meshes)) {}

  explicit ReactionPDEProblem(std::string meshName) : IProblem(std::move(meshName)) {}
  
  ReactionPDEProblem(ReactionPDEProblem&);

  Scalar Load(const Point &x) const override {
    return loadFunction(x);
  }

  VectorField Flux(const Point &x) const override {
    return fluxFunction(x);
  }

  Tensor Diffusion(const Point &x) const override {
    return diffusionFunction(x);
  }

  VectorField Convection(const Point &x) const override {
    return convectionFunction(x);
  }

  Scalar Reaction(const Point &x, double c) const override {
    return reactionFunction(x, c);
  }

  Scalar Concentration(double t, const Point &x) const override {
    return concentrationFunction(t, x);
  }

  Scalar DerivativeReaction(const Point &x, double c) const override {
    return derivativeReactionFunction(x, c);
  }

  VectorField CellConvection(const Cell &c, const Point &x) const {
    if (cellConvectionFunction) return cellConvectionFunction(c, x);
    return Convection(x);
  }

  double FaceConvection(const Cell &c, int face,
                                const VectorField &N, const Point &x) const {
    if (faceConvectionFunction) return faceConvectionFunction(c, face, N, x);
    return Convection(x) * N;
  }

  std::string Name() const override {
    if (nameFunction) return nameFunction();
    return "ReactionPDEProblem";
  }
};

class ReactionPDEProblemBuilder {
private:
  ReactionPDEProblem::LoadFunction loadFunction = nullptr;
  ReactionPDEProblem::FluxFunction fluxFunction = nullptr;
  ReactionPDEProblem::DiffusionFunction diffusionFunction = nullptr;
  ReactionPDEProblem::ConvectionFunction convectionFunction = nullptr;
  ReactionPDEProblem::ReactionFunction reactionFunction = nullptr;
  ReactionPDEProblem::ConcentrationFunction concentrationFunction = nullptr;
  ReactionPDEProblem::DerivativeReactionFunction derivativeReactionFunction = nullptr;
  ReactionPDEProblem::CellConvectionFunction cellConvectionFunction = nullptr;
  ReactionPDEProblem::FaceConvectionFunction faceConvectionFunction = nullptr;
  ReactionPDEProblem::NameFunction nameFunction = nullptr;

  std::string meshName = "";
  std::shared_ptr<Meshes> meshes;

public:
  explicit ReactionPDEProblemBuilder(){}

  [[nodiscard]] ReactionPDEProblemBuilder &
  WithLoad(const ReactionPDEProblem::LoadFunction func) {
    loadFunction = func;
    return *this;
  }

  [[nodiscard]] ReactionPDEProblemBuilder &
  WithFlux(const ReactionPDEProblem::FluxFunction func) {
    fluxFunction = func;
    return *this;
  }

  [[nodiscard]] ReactionPDEProblemBuilder &
  WithDiffusion(const ReactionPDEProblem::DiffusionFunction func) {
    diffusionFunction = func;
    return *this;
  }

  [[nodiscard]] ReactionPDEProblemBuilder &
  WithConvection(const ReactionPDEProblem::ConvectionFunction func) {
    convectionFunction = func;
    return *this;
  }

  [[nodiscard]] ReactionPDEProblemBuilder &
  WithReaction(const ReactionPDEProblem::ReactionFunction func) {
    reactionFunction = func;
    return *this;
  }

  [[nodiscard]] ReactionPDEProblemBuilder &
  WithConcentration(const ReactionPDEProblem::ConcentrationFunction func) {
    concentrationFunction = func;
    return *this;
  }

  [[nodiscard]] ReactionPDEProblemBuilder &
  WithDerivativeReaction(const ReactionPDEProblem::DerivativeReactionFunction func) {
    derivativeReactionFunction = func;
    return *this;
  }

  [[nodiscard]] ReactionPDEProblemBuilder &
  WithCellConvection(const ReactionPDEProblem::CellConvectionFunction func) {
    cellConvectionFunction = func;
    return *this;
  }

  [[nodiscard]] ReactionPDEProblemBuilder &
  WithFaceConvection(const ReactionPDEProblem::FaceConvectionFunction func) {
    faceConvectionFunction = func;
    return *this;
  }

  [[nodiscard]] ReactionPDEProblemBuilder &
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
  [[nodiscard]] ReactionPDEProblemBuilder &
  WithMeshes(std::shared_ptr<Meshes> sharedMeshes) {
    meshes = sharedMeshes;
    return *this;
  }

  std::unique_ptr<ReactionPDEProblem> Build() {
    std::unique_ptr<ReactionPDEProblem> problem;
    if (meshes) {
      problem = std::make_unique<ReactionPDEProblem>(meshes);
    } else {
      problem = std::make_unique<ReactionPDEProblem>(meshName);
    }

    problem->loadFunction = std::move(loadFunction);
    problem->fluxFunction = std::move(fluxFunction);
    problem->diffusionFunction = std::move(diffusionFunction);
    problem->convectionFunction = std::move(convectionFunction);
    problem->diffusionFunction = std::move(diffusionFunction);
    problem->reactionFunction = std::move(reactionFunction);
    problem->concentrationFunction = std::move(concentrationFunction);
    problem->derivativeReactionFunction = std::move(derivativeReactionFunction);
    problem->cellConvectionFunction = std::move(cellConvectionFunction);
    problem->faceConvectionFunction = std::move(faceConvectionFunction);
    problem->nameFunction = std::move(nameFunction);

    return problem;
  }
};

IReactionProblem *
CreateReactionProblem(const std::string &problemName);

std::unique_ptr<IReactionProblem>
CreateReactionProblemUnique(const std::string &problemName);

std::shared_ptr<IReactionProblem>
CreateReactionProblemShared(const std::string &problemName);

#endif //REACTIONPROBLEMS_HPP
