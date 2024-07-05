#ifndef GELFAND_HPP
#define GELFAND_HPP

#include "Algebra.hpp"
#include "Assemble.hpp"
#include "GelfandProblem.hpp"
#include "LagrangeDiscretization.hpp"
#include "PDESolver.hpp"

class GelfandAssemble : public IAssemble {
protected:
  const GelfandProblem &problem;
  std::shared_ptr <const LagrangeDiscretization> disc;
  
public:
  GelfandAssemble(const GelfandProblem &p, int degree = 1) :
    problem(p),
    disc(std::make_shared<LagrangeDiscretization>(problem.GetMeshes(), degree)) {}

  const IDiscretization &GetDisc() const { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const { return disc; }

  const char *Name() const { return problem.Name().c_str(); }
  
  void Initialize(Vector &u) const;
    
  double Energy(const Vector &u) const;

  void Residual(const cell &c, const Vector &u, Vector &r) const;

  void Jacobi(const cell &c, const Vector &u, Matrix &A) const;
};

GelfandAssemble *CreateGelfandAssemble(const GelfandProblem &problem, const PDESolverConfig &conf);

std::unique_ptr<GelfandAssemble> CreateGelfandAssembleUnique(const GelfandProblem &problem,
                                                             const PDESolverConfig &conf);

std::shared_ptr<GelfandAssemble> CreateGelfandAssembleShared(const GelfandProblem &problem,
                                                             const PDESolverConfig &conf);

#endif //GELFAND_HPP
