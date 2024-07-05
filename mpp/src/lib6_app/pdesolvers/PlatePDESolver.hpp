#include "EigenSolverCreator.hpp"
#include "Newton.hpp"
#include "PlateAssemble.hpp"

class PlatePDESolver : public PDESolver<PlateProblem> {
private:
  std::shared_ptr<Newton> newton;

  std::shared_ptr<PlateAssemble> assemble;
protected:
  void run(Solution &solution) const override;

  void plotVtu(Solution &solution) const override;

  void computeValues(Solution &solution) const override;

  void createAssemble(std::shared_ptr<PlateProblem> problem) override;
public:
  explicit PlatePDESolver(const PDESolverConfig &conf);

  std::string Name() const override { return "PlatePDESolver"; }

  // From PDEMain
  ValueMap ComputeValues(const Vector &u) const override;

  void PrintValues(const Solution &solution) override;

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return assemble->GetSharedDisc();
  }
};

// int PlatePDESolver () {
//   Config::PrintInfo();
//   PlateProblem problem;
//   PlateAssemble assemble(problem);
//   assemble.GetDisc().PrintInfo();
//   mout << "Problem " << problem.Name() << endl;
//   Newton newton;
//   Vector u(0.0, assemble.GetSharedDisc());
//   u.PrintInfo();
//   u.GetMesh().PrintInfo();
//   mpp::plot("u0") << u << mpp::endp;
//   newton(assemble, u);
//   mpp::plot("u") << u << mpp::endp;
//   int N = 10;
//   Config::Get("eigenvalues", N);
//   Vectors U(N, assemble.GetSharedDisc());
//   for (int i=0; i<N; ++i)
//     assemble.Initialize(U[i]);
//   Matrix A(u);
//   Matrix B(u,false);
//   assemble.Matrices(A,B);
//   IEigenSolver *esolver = EigenSolverCreator("LOBPCG").Create();
//   Eigenvalues lambda;
//   (*esolver)(U, lambda, A, B);
//   delete esolver;
//   for (int i=0; i<N; ++i)
//     mout << "eigenvalue[" << i << "] = " << lambda[i] << endl;
//   for (int i=0; i<N; ++i) {
//     std::string filename = "U." + to_string(i);
//     mpp::plot(filename) << U[i] << mpp::endp;
//   }
//   return 0;
// }