//
// Created by lstengel on 27.07.22.
//

#ifndef TUTORIAL_DGVECTORVALUEDASSEMBLE_H
#define TUTORIAL_DGVECTORVALUEDASSEMBLE_H



#include "IVectorValuedAssemble.hpp"
#include "DGDiscretization.hpp"
#include "DGVectorFieldElement.hpp"
#include "ElementPool.hpp"


class DGVectorValuedAssemble : public IVectorValuedAssemble {
protected:
  mutable ElementPool<DGVectorFieldElement> elementPool;
  mutable FaceElementPool<DGVectorFieldFaceElement> faceElementPool;

  std::shared_ptr<const DGDiscretization> disc;

  double penalty = 20.0;

  int sign = 1;

  int degree = disc->Degree();


public:
  DGVectorValuedAssemble(const IVectorValuedProblem &problem, int degree) :
      disc(std::make_shared<const DGDiscretization>(problem.GetMeshes(), degree, problem.GetMeshes().dim())),
      IVectorValuedAssemble(problem) {
    Config::Get("penalty", penalty);
    Config::Get("sign", sign);
  }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override { return disc;};

  const IDiscretization &GetDisc() const override { return *disc; }

  void PrintInfo() const override {
    mout.PrintInfo("Assemble", verbose,
                   PrintInfoEntry("Assemble Name", Name()),
                   PrintInfoEntry("penalty", penalty),
                   PrintInfoEntry("sign", sign),
                   PrintInfoEntry("Problem", problem.Name()));
  }

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


#endif //TUTORIAL_DGVECTORVALUEDASSEMBLE_H
