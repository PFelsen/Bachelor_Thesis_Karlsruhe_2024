#ifndef DG_VECTORVALUED_ASSEMBLE_HPP
#define DG_VECTORVALUED_ASSEMBLE_HPP

#include "IVectorValuedAssemble.hpp"
#include "MixedEGDiscretization.hpp"
#include "MixedEGVectorFieldElement.hpp"
#include "ElementPool.hpp"


class EGVectorValuedAssemble : public IVectorValuedAssemble {
protected:
  mutable ElementPool<MixedEGVectorFieldElement> elementPool;
  mutable FaceElementPool<MixedEGVectorFieldFaceElement> faceElementPool;

  std::shared_ptr<const MixedEGDiscretization> disc;

  int degree = disc->Degree();

  double penalty = 100.0;

  double sign = -1.0;

public:
  EGVectorValuedAssemble(const IVectorValuedProblem &problem, int degree) :
      disc(std::make_shared<const MixedEGDiscretization>(problem.GetMeshes(), degree,
                                                         problem.GetMeshes().dim())),
      IVectorValuedAssemble(problem) {
  }

  const IDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override { return disc; };

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

#endif //EGVECTORVALUEDASSEMBLE_HPP
