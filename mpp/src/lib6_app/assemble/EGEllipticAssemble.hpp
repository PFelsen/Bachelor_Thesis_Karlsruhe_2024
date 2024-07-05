#ifndef EGELLIPTICASSEMBLE_HPP
#define EGELLIPTICASSEMBLE_HPP

#include "IEllipticAssemble.hpp"
#include "EllipticAssemble.hpp"
#include "MixedEGDiscretization.hpp"
#include "MixedEGElement.hpp"
#include "ElementPool.hpp"


class EGEllipticAssemble : public EllipticAssemble<MixedEGElement, MixedEGFaceElement> {
protected:
  std::shared_ptr<const MixedEGDiscretization> disc;

  int degree = disc->Degree();

  double penalty = 20.0;

  double sign = -1.0;

public:
  EGEllipticAssemble(const IEllipticProblem &problem, int degree) :
      disc(std::make_shared<const MixedEGDiscretization>(problem.GetMeshes(), degree)),
      EllipticAssemble(problem) {
    Config::Get("sign", sign);
    Config::Get("penalty", penalty);
  }

  const IDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override { return disc; }

  const char *Name() const override;

  void Initialize(Vector &u) const override;

  void Residual(const cell &c, const Vector &u, Vector &r) const override;

  void Jacobi(const cell &c, const Vector &u, Matrix &A) const override;

  FluxPair InflowOutflow(const Vector &u) const override;

  FluxPair PrescribedInflowOutflow(const Vector &u) const override;

  FluxPair OutflowLeftRight(const Vector &u) const override;

  double GoalFunctional(const Vector &u) const override;

  void SetExactSolution(Vector &uEx) const override;

  void SetFlux(const Vector &u, Vector &flux) override;

  void SetPressure(const Vector &u, Vector &p) const override;
  
  Point FindDirichlet(Vector &u) const;
};

#endif //EGELLIPTICASSEMBLE_HPP
