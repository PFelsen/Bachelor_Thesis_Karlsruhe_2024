//
// Created by lstengel on 21.07.22.
//

#ifndef TUTORIAL_IVECTORVALUEDASSEMBLE_H
#define TUTORIAL_IVECTORVALUEDASSEMBLE_H

#include "EllipticProblems.hpp"
#include "Assemble.hpp" // kann die verwendet werden?
#include "PDESolver.hpp"


class IVectorValuedAssemble : public IAssemble {
protected:
  const IVectorValuedProblem &problem;

public:
  explicit IVectorValuedAssemble(const IVectorValuedProblem &problem) : problem(problem) {};

  const IVectorValuedProblem &GetProblem() const { return problem; };

  virtual std::shared_ptr<const IDiscretization> GetSharedDisc() const = 0;

  virtual const IDiscretization &GetDisc() const = 0;

  const char *Name() const override = 0;

  void Initialize(Vector &u) const override = 0;

  double Energy(const Vector &u) const override = 0;

  void Residual(const cell &c, const Vector &u, Vector &r) const override = 0;

    void Jacobi(const cell &c, const Vector &u, Matrix &A) const override = 0;

    virtual double EnergyError(const Vector &u) const = 0;

    virtual double L2(const Vector &u) const = 0;

    virtual double H1(const Vector &u) const = 0;

    virtual double L2Error(const Vector &u) const = 0;

    virtual double L2CellAvgError(const Vector &u) const = 0;

    virtual double MaxError(const Vector &u) const = 0;

    virtual double FluxError(const Vector &u) const = 0;

    virtual double FaceError(const Vector &u) const = 0;

    virtual void SetExactSolution(Vector &u_ex) const = 0;

    virtual void SetPressure(const Vector &u, Vector &p) const { p = u; }; // soll das drin bleiben?

    void PrintInfo() const override {
      mout.PrintInfo("VectorValued Assemble", verbose,
                     PrintInfoEntry("VectorValued Assemble Name", Name()),
                     PrintInfoEntry("VectorValued Problem", problem.Name()));
    }
  };

//Probleme anpassen
IVectorValuedAssemble *CreateVectorValuedAssemble(IVectorValuedProblem &problem,
                                                  const PDESolverConfig &conf);

std::unique_ptr<IVectorValuedAssemble> CreateVectorValuedAssembleUnique(IVectorValuedProblem &problem,
                                                                        const PDESolverConfig &conf);

std::shared_ptr<IVectorValuedAssemble> CreateVectorValuedAssembleShared( IVectorValuedProblem &problem,
                                                                         const PDESolverConfig &conf);



#endif //TUTORIAL_IVECTORVALUEDASSEMBLE_H


