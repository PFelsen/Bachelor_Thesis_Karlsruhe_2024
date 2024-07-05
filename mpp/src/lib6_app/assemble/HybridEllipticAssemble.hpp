#ifndef _HYBRID_H_
#define _HYBRID_H_

#include <utility>

#include "IEllipticAssemble.hpp"
#include "RTLagrangeElement.hpp"
#include "RTElement.hpp"
#include "ElementPool.hpp"

class HybridEllipticAssemble : public IEllipticAssemble {
protected:
  mutable FaceElementPool<RTFaceElementT<>> faceElementPool;

  mutable ElementPool<RTLagrangeElementT<>> elementPool;

  std::shared_ptr<const HybridRTLagrangeDiscretization> disc;

public:
  HybridEllipticAssemble(const IEllipticProblem &problem) :
      disc(std::make_shared<HybridRTLagrangeDiscretization>(problem.GetMeshes(), 0, 0)),
      IEllipticAssemble(problem) {}

  const IDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override { return disc; }

  const char *Name() const override;

  void Initialize(Vector &u) const override;

  double Energy(const Vector &u) const override;

  void Residual(const cell &c, const Vector &u, Vector &r) const override;

  void Jacobi(const cell &c, const Vector &u, Matrix &A) const override;

  double EnergyError(const Vector &u) const override;

  double L2(const Vector &u) const override;

  double H1(const Vector &u) const override;

  double L2Error(const Vector &u) const override;

  double L2CellAvgError(const Vector &u) const override;

  double MaxError(const Vector &u) const override;

  double FluxError(const Vector &u) const override;

  double FaceError(const Vector &u) const override;
  
  FluxPair InflowOutflow(const Vector &u) const override;

  FluxPair PrescribedInflowOutflow(const Vector &u) const override;

  FluxPair OutflowLeftRight(const Vector &u) const override;

  double GoalFunctional(const Vector &u) const override;

  void SetExactSolution(Vector &uEx) const override;

  void SetFlux(const Vector &u, Vector &flux) override;

  void SetPressureFlux(const Vector &u, Vector &p, Vector &flux) const;

  double DualPrimalError(const Vector &u) const override;

  void LocalProblem(const cell &c, const Vector &u,
                    RMatrix &A, RVector &B, RVector &R) const;

  void SetNormalFlux(const Vector &u, Vector &flux) override;

  VectorField EvaluateCellFlux(const Vector &flux, const Cell &c) const;

  double EvaluateNormalFlux(const Vector &flux, const Cell &c, int i) const;

  VectorField EvaluateNormal(const Vector &flux, const cell &c, int i) const;

  Scalar Value(const cell &c, const Vector &u, RTLagrangeElement &elem) const;

  Scalar FaceValue(const cell &c, const Vector &u, RTFaceElementT<> &faceElem) const;

  VectorField Flux(const cell &c, const Vector &u, RTLagrangeElement &elem, int q) const;

  VectorField FaceFlux(const cell &c, const Vector &u, const RTFaceElementT<> &faceElem,
                       int q, int face) const;
};

#endif
