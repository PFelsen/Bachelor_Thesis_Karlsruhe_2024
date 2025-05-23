#ifndef LAGRANGE_VECTORVALUED_ASSEMBLE_HPP
#define LAGRANGE_VECTORVALUED_ASSEMBLE_HPP


#include "IVectorValuedAssemble.hpp"
#include "LagrangeDiscretization.hpp"
#include "VectorFieldElement.hpp"


class LagrangeVectorValuedAssemble : public IVectorValuedAssemble {
protected:
  std::shared_ptr<const LagrangeDiscretization> disc;

public:
  LagrangeVectorValuedAssemble(const IVectorValuedProblem &problem, int degree) :
      disc(std::make_shared<const LagrangeDiscretization>(problem.GetMeshes(), degree, problem.GetMeshes().dim())),
      IVectorValuedAssemble(problem) {}

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override { return disc; };

  const IDiscretization &GetDisc() const override { return *disc; }

  const char *Name() const override;

  void Initialize(Vector &u) const override;

  double Energy(const Vector &u) const override;

  void Residual(const cell &c, const Vector &u, Vector &r) const override;

  void Jacobi(const cell &c, const Vector &u, Matrix &A) const override;

  virtual double EnergyError(const Vector &u) const override;

  double L2(const Vector &u) const override;

  double H1(const Vector &u) const override;

  double L2Error(const Vector &u) const override;

  double L2CellAvgError(const Vector &u) const override;

  double MaxError(const Vector &u) const override;

  double FluxError(const Vector &u) const override;

  double FaceError(const Vector &u) const override;

  virtual void SetExactSolution(Vector &uEx) const override;

};

#endif //TUTORIAL_VECTORVALUEDLAGRANGEELLIPTIC_HPP



