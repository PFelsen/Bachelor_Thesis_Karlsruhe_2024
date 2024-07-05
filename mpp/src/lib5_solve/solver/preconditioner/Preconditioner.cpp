#include "GaussSeidel.hpp"
#include "Jacobi.hpp"
#include "Multigrid.hpp"
#include "Richardson.hpp"
#include "SSOR.hpp"
#ifdef USE_SUPERLU
#include "ParallelSolver.hpp"
#include "SuperLU.hpp"
#endif

class NoPreconditioner : public Preconditioner {
public:
  NoPreconditioner() {}

  virtual ~NoPreconditioner() {}

  void Construct(const Matrix &A) override {}

  void Destruct() override {}

  void multiply(Vector &u, const Vector &b) const override {
    u = b;
    u.Accumulate();
  }

  string Name() const override { return "NoPreconditioner"; }
};

std::string PCName(const string &prefix) {
  std::string preCondName = "SuperLU";
  Config::Get(prefix + "Preconditioner", preCondName);
  return preCondName;
}

Preconditioner *GetPCByPrefix(const string &prefix) { return GetPC(PCName(prefix)); }

Preconditioner *GetPC(const string &name) {
  if (name == "Richardson") return new Richardson();
  if (name == "Jacobi") return new Jacobi();
  if (name == "PointBlockJacobi") return new PointBlockJacobi();
  if (name == "PointBlockJacobi_dG") {
    WarningOnMaster("Stop using PointBlockJacobi_dG, it is programmed for a very special case. Can "
                    "be removed according to wieners");
    return new PointBlockJacobi_dG();
  }
  if (name == "DampedJacobi") {
    WarningOnMaster(
        "DampedJacobi offers no new functionality, just use Jacobi") return new DampedJacobi();
  }
  if (name == "CellWiseJacobi") return new CellWiseJacobi();
  if (name == "NoPreconditioner") return new NoPreconditioner();
  if (name == "GaussSeidel") return new GaussSeidel();
  if (name == "PointBlockGaussSeidel") return new PointBlockGaussSeidel();
  if (name == "PointBlockGaussSeidelBackwards") return new PointBlockGaussSeidel(false);
  if (name == "BlockGaussSeidel") return new BlockGaussSeidel();
  if (name == "SSOR") return new SSOR();
  if (name == "SGS") { WarningOnMaster("Use SSOR with omega=1 instead") return new SSOR(1.0); }
  if (name == "Multigrid") return new Multigrid();
#ifdef USE_SUPERLU
  if (name == "SuperLU") return new SuperLU();
  if (name == "SuperLU_local") return new SuperLU_local();
  if (name == "PS") return new ParallelSolver();
  if (name == "LIB_PS") {
    WarningOnMaster(
        "Preconditioner LIB_PS is deprecated: Use PS instead!") return new ParallelSolver();
  }
#endif
  THROW("no preconditioner " + name + " implemented")
}

Preconditioner *GetPC() { return GetPC(PCName("")); }
