#ifndef PLATE_HPP
#define PLATE_HPP

#include "Algebra.hpp"
#include "ArgyrisDiscretization.hpp"
#include "Assemble.hpp"
#include "PDESolver.hpp"
#include "PlateProblem.hpp"

class PlateAssemble : public IAssemble {
protected:
  const PlateProblem &problem;
  std::shared_ptr <const ArgyrisDiscretization> disc;
  
public:
  PlateAssemble(const PlateProblem &p) :
    problem(p),
    disc(std::make_shared<ArgyrisDiscretization>(problem.GetMeshes())) {}

  const IDiscretization &GetDisc() const { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const { return disc; }

  const char *Name() const { return problem.Name().c_str(); }
  
  void Initialize(Vector &u) const;
    
  double Energy(const Vector &u) const;

  void Residual(const cell &c, const Vector &u, Vector &r) const;

  void Jacobi(const cell &c, const Vector &u, Matrix &A) const;

  void Matrices(Matrix &A, Matrix& B) const;  
};

PlateAssemble *CreatePlateAssemble(const PlateProblem &problem, const PDESolverConfig &conf);

std::unique_ptr<PlateAssemble> CreatePlateAssembleUnique(const PlateProblem &problem,
                                                         const PDESolverConfig &conf);

std::shared_ptr<PlateAssemble> CreatePlateAssembleShared(const PlateProblem &problem,
                                                         const PDESolverConfig &conf);

#endif //PLATE_HPP
