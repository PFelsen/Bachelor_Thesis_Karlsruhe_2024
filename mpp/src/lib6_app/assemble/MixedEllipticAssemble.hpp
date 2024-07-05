#ifndef _MIXED_H_
#define _MIXED_H_

#include "RTLagrangeDiscretization.hpp"
#include "LagrangeDiscretization.hpp"
#include "IEllipticAssemble.hpp"
#include "RTLagrangeElement.hpp"
#include "RTElement.hpp"
#include "ElementPool.hpp"


class MixedEllipticAssemble : public IEllipticAssemble {
protected:
  mutable FaceElementPool<RTFaceElementT<>> faceElementPool;

  mutable ElementPool<RTLagrangeElementT<>> elementPool;

  std::shared_ptr<const RTLagrangeDiscretization> disc;

public:
  MixedEllipticAssemble(const IEllipticProblem &problem) :
      disc(std::make_shared<RTLagrangeDiscretization>(problem.GetMeshes(), 0, 0)),
      IEllipticAssemble(problem) {};

  const IDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override { return disc; }

  const char *Name() const override;

  void Initialize(Vector &u) const override;

  double Energy(const Vector &u) const override;

  void Residual(const cell &c, const Vector &u, Vector &r) const override;

  void Jacobi(const cell &c, const Vector &u, Matrix &A) const override;

  double EnergyError(const Vector &u) const override;

  double L2(const Vector &u) const override;

  double H1(const Vector &u) const;

  double L2Error(const Vector &u) const override;

  double L2CellAvgError(const Vector &u) const override;

  double MaxError(const Vector &u) const override;

  double FluxError(const Vector &u) const override;

  double FaceError(const Vector &u) const override;

  FluxPair InflowOutflow(const Vector &u) const override;

  FluxPair PrescribedInflowOutflow(const Vector &u) const override;

  FluxPair OutflowLeftRight(const Vector &u) const override;

  double GoalFunctional(const Vector &u) const override;

  void SetPressure(const Vector &u, Vector &p) const override;

  void SetExactSolution(Vector &uEx) const override;

  void SetFlux(const Vector &u, Vector &flux) override;

  double DualPrimalError(const Vector &u) const override;

  void SetPressureFlux(const Vector &u, Vector &p, Vector &flux) const;

  VectorField EvaluateCellFlux(const Vector &flux, const Cell &c) const;

  double EvaluateNormalFlux(const Vector &flux, const Cell &c, int i) const;

  void SetNormalFlux(const Vector &u, Vector &flux);

};

#endif
