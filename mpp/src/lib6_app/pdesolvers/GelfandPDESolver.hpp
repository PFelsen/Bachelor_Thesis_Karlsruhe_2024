#include "GelfandAssemble.hpp"
#include "Newton.hpp"

class GelfandPDESolver : public PDESolver<GelfandProblem> {
private:
  std::shared_ptr<Newton> newton;

  std::shared_ptr<GelfandAssemble> assemble;
protected:
  void run(Solution &solution) const override;

  void plotVtu(Solution &solution) const override;

  void computeValues(Solution &solution) const override;

  void createAssemble(std::shared_ptr<GelfandProblem> problem) override;
public:
  explicit GelfandPDESolver(const PDESolverConfig &conf);

  std::string Name() const override { return "GelfandPDESolver"; }

  // From PDEMain
  ValueMap ComputeValues(const Vector &u) const override;

  void PrintValues(const Solution &solution) override;

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return assemble->GetSharedDisc();
  }
};

// int GelfandPDolver () {
//   Config::PrintInfo();
//   GelfandProblem problem;
//   GelfandAssemble assemble(problem,2);
//   Newton newton;
//   Vector u(0.0, assemble.GetSharedDisc());
//   u.PrintInfo();
//   assemble.Initialize(u);
//   mpp::plot("u0") << u << mpp::endp;
//   newton(assemble, u);
//   mpp::plot("u") << u << mpp::endp;
//   mout << "max(u) = " << u.Max() << endl;
//   return 0;
// }
