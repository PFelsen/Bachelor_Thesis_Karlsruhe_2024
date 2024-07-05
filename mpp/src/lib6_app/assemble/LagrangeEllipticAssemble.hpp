#ifndef _LAPLACE_H_
#define _LAPLACE_H_

#include "IEllipticAssemble.hpp"
#include "EllipticAssemble.hpp"
#include "LagrangeDiscretization.hpp"
#include "ScalarElement.hpp"
#include "ElementPool.hpp"


class LagrangeEllipticAssemble : public EllipticAssemble<ScalarElement, ScalarFaceElement> {
protected:
  std::shared_ptr<const LagrangeDiscretization> disc;

public:
  LagrangeEllipticAssemble(const IEllipticProblem &problem, int degree) :
      disc(std::make_shared<const LagrangeDiscretization>(problem.GetMeshes(), degree)),
      EllipticAssemble(problem) {}

  const IDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override { return disc; }

  const char *Name() const override;

  void Initialize(Vector &u) const override;

  void Residual(const cell &c, const Vector &u, Vector &r) const override;

  void Jacobi(const cell &c, const Vector &u, Matrix &A) const override;

  FluxPair InflowOutflow(const Vector &u) const override;

  FluxPair PrescribedInflowOutflow(const Vector &u) const override;

  FluxPair OutflowLeftRight(const Vector &u) const override;

  virtual double GoalFunctional(const Vector &u) const override;

  virtual void SetExactSolution(Vector &uEx) const override;

  virtual void SetFlux(const Vector &u, Vector &flux) override;
};

#endif

