#ifndef SPACETIME_STTRANSPORTPROBLEMS_HPP
#define SPACETIME_STTRANSPORTPROBLEMS_HPP

#include "ProblemBase.hpp"

class TransportProblem : public ProblemBase {
protected:
  int dim;
  int verbose = 0;

  const std::vector<COMPONENT> components;

  static std::vector<COMPONENT> createComponents() {
    std::vector<COMPONENT> components{COMPONENT::V_X};
    return components;
  }

public:
  TransportProblem(int dim) : dim(dim), components(createComponents()) {
    Config::Get("ProblemVerbose", verbose);
  }

  virtual VectorField FluxVector(const LevelPair &level, const cell c, const Point &x) const = 0;

  virtual double NormalFlux(const LevelPair& level,
                            const cell c,
                            const Point &x,
                            const VectorField &N,
                            int face) const {
    return N * FluxVector(level, c, x);
  };
    
  virtual double divFluxVector(const LevelPair& level, const cell c, const Point &x) const { return 0; }

  VectorField Flux(const LevelPair& level, const cell c, const Point &x, double u) const {
    return u * FluxVector(level, c, x);
  }

  double divFlux(const LevelPair &level, const cell c, const Point &x, double &u, VectorField &Du) const {
    return u * divFluxVector(level, c, x) + FluxVector(level, c, x) * Du;
  }

  virtual double ut(const cell c, const Point &x) const {
    return 0.0;
  }

  virtual double InflowData(const LevelPair& level, const cell c, const Point &x, VectorField &N) const { return 0.0; }

  bool HasInflow(const LevelPair& level, const cell c, const Point &x, VectorField N) const { return (FluxVector(level, c, x) * N < 0); }

  bool HasOutflow(const LevelPair& level, const cell c, const Point &x, VectorField N) const { return (FluxVector(level, c, x) * N >= 0); }

  virtual bool hasInitialData() const { return false; }

  virtual bool HasRHS() const { return false; }

  virtual bool HasExactSolution() const { return false; }

  int Dim() const { return 1; }

  std::vector<COMPONENT> GetComponents() const override { return components; }

  const Meshes& GetFluxMeshes() const {
    THROW("GetFluxProblem not implemented.");
  }

  virtual const IProblem& GetFluxProblem() const {
    THROW("No FluxProblem implemented.")
  }

};

TransportProblem *CreateTransportProblemST(const string &name);

std::unique_ptr<TransportProblem> CreateTransportProblemUniqueST(const string &name);

std::shared_ptr<TransportProblem> CreateTransportProblemSharedST(const string &name);


#endif //SPACETIME_STTRANSPORTPROBLEMS_HPP
