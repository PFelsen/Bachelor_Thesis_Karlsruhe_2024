#ifndef EIGENVALUEEXAMPLEUTIL_HPP
#define EIGENVALUEEXAMPLEUTIL_HPP

#include "EigenvalueExampleAssemble.hpp"

#include "MeshesCreator.hpp"
#include "IAEigenvalueMethods.hpp"

std::vector<IAInterval> ComputeEigenvaluesBaseProblem(int N_ev) {
  std::vector<IAInterval> eigenvalues{};
  for (int i = 1; i <= N_ev; ++i)
    for (int j = 1; j <= N_ev; ++j)
      eigenvalues.push_back((i * i + j * j) * IAInterval::PiSqr() - 1.0);

  std::sort(eigenvalues.begin(), eigenvalues.end(),
            [](const IAInterval &a, const IAInterval &b) -> bool { return inf(a) < inf(b); });

  return std::vector<IAInterval>(eigenvalues.begin(), eigenvalues.begin() + N_ev);
}

template<bool VERIFIED>
IAUpperEigenvalueBounds ApplyRayleighRitz(const Meshes &meshes, EigenvalueExampleAssemble<VERIFIED> &assemble,
                                          std::shared_ptr<const IDiscretization> disc_ev,
                                          std::shared_ptr<const IAIDiscretization> ia_disc_ev,
                                          int N_ev, double shift, double t = 0.0) {
  Eigenfcts eigenfcts(N_ev, 0.0, disc_ev);
  eigenfcts.SetIADisc(ia_disc_ev);
  IAUpperEigenvalueBounds Lambda;

  IARayleighRitzMethod<VERIFIED> rayleighRitzMethod;
  rayleighRitzMethod(assemble, eigenfcts, Lambda, t);

  mout << "Verified upper eigenvalue bounds (t= " << t << "):" << endl;
  for (int i = 0; i < Lambda.size(); ++i)
    mout << Lambda[i] << "  ";
  mout << endl << endl;

  return Lambda;
}

template<bool VERIFIED>
IALowerEigenvalueBounds ApplyLehmannGoerisch(const Meshes &meshes, EigenvalueExampleAssemble<VERIFIED> &assemble,
                                             std::shared_ptr<const IDiscretization> disc_ev,
                                             std::shared_ptr<const IAIDiscretization> ia_disc_ev,
                                             int N_ev, double shift, double rho, double t = 0.0) {
  Eigenpairs eigenpairs(N_ev, disc_ev);
  eigenpairs.SetIADisc(ia_disc_ev);
  IALowerEigenvalueBounds lambda;

  IALehmannGoerischMethod<VERIFIED> lehmannGoerischMethod;
  lehmannGoerischMethod(assemble, eigenpairs, lambda, rho, t);

  mout << "Verified lower eigenvalue bounds (t= " << t << "):" << endl;
  for (int i = 0; i < lambda.size(); ++i)
    mout << lambda[i] << "  ";
  mout << endl << endl;

  return lambda;
}

template<bool VERIFIED>
double ApplyHomotopyMethod(const Meshes &meshes, EigenvalueExampleAssemble<VERIFIED> &assemble,
                           std::shared_ptr<const IDiscretization> disc_ev,
                           std::shared_ptr<const IAIDiscretization> ia_disc_ev,
                           int N_ev, double shift, double &rho) {
  Eigenpairs eigenpairs(N_ev, disc_ev);
  eigenpairs.SetIADisc(ia_disc_ev);

  IAHomotopyMethod<VERIFIED> homotopyMethod;
  return homotopyMethod(assemble, eigenpairs, rho, true);
}

template<bool VERIFIED>
IAEigenvalueEnclosures ComputeEnclosures(const Meshes &meshes, EigenvalueExampleAssemble<VERIFIED> &assemble,
                                         std::shared_ptr<const IDiscretization> disc_ev,
                                         std::shared_ptr<const IAIDiscretization> ia_disc_ev,
                                         int N_ev, double shift, double &rho, double step) {
  IAUpperEigenvalueBounds Lambda1 = ApplyRayleighRitz<VERIFIED>(meshes, assemble, disc_ev, ia_disc_ev, N_ev, shift, step);
  IALowerEigenvalueBounds lambda1 = ApplyLehmannGoerisch<VERIFIED>(meshes, assemble, disc_ev, ia_disc_ev, N_ev, shift, rho, step);
  return IAEigenvalueEnclosures(lambda1, Lambda1);
}

template<bool VERIFIED>
void RunExample() {
  int N_ev = 6;
  double shift = 2.0;
  double rho = 127.0;

  std::unique_ptr<Meshes> meshes = MeshesCreator("Square2Triangles").
      WithPLevel(1).
      WithLevel(5).CreateUnique();
  meshes->PrintInfo();

  EigenvalueExampleAssemble<VERIFIED> assemble(shift);
  std::shared_ptr<IDiscretization> disc_ev = assemble.CreateEVDisc(*meshes);
  std::shared_ptr<IAIDiscretization> ia_disc_ev = assemble.CreateIAEVDisc(*meshes);

  IAUpperEigenvalueBounds Lambda = ApplyRayleighRitz<VERIFIED>(*meshes, assemble, disc_ev, ia_disc_ev, N_ev, shift);
  IALowerEigenvalueBounds lambda = ApplyLehmannGoerisch<VERIFIED>(*meshes, assemble, disc_ev, ia_disc_ev, N_ev, shift, rho);

  std::vector<IAInterval> eigenvalues =  ComputeEigenvaluesBaseProblem(N_ev);
  IAEigenvalueEnclosures IALambda(lambda, Lambda);
  mout << "Exact eigenvalues and verified enclosures for the base problem (i.e., with t= 0):" << endl;
  for (int i = 0; i < IALambda.size(); ++i)
    mout << eigenvalues[i] << " \\subseteq " << IALambda[i] << endl;
  mout << endl << endl;

  double final_step = ApplyHomotopyMethod<VERIFIED>(*meshes, assemble, disc_ev, ia_disc_ev, N_ev, shift, rho);

  IAEigenvalueEnclosures IALambda1 = ComputeEnclosures<VERIFIED>(*meshes, assemble, disc_ev, ia_disc_ev, N_ev, shift, rho, final_step);
  mout << "Verified eigenvalue enclosures for the eigenvalue problem (with t= " << final_step << "):" << endl;
  for (int i = 0; i < IALambda1.size(); ++i)
    mout << IALambda1[i] << "  ";
  mout << endl << endl;
}

#endif //EIGENVALUEEXAMPLEUTIL_HPP
