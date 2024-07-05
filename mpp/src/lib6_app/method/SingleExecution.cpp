#include "SingleExecution.hpp"
#include "EllipticPDESolver.hpp"
#include "TransportPDESolver.hpp"
#include "ReactionPDESolver.hpp"
#include "VectorValuedMain.hpp"
#include "AcousticPDESolver.hpp"
#include "GelfandPDESolver.hpp"
#include "PlatePDESolver.hpp"

#ifdef USE_SPACETIME

#include "STAcousticPDESolver.hpp"
#include "STAcousticMain.hpp"
#include "AcousticProblems.hpp"
#include "STTransportMain.hpp"

#endif


void SingleExecution::Method() {
  if (conf.modelName.find("Elliptic") != std::string::npos) {
    EllipticPDESolver pdeSolver(conf);
    auto problem = CreateEllipticProblemShared(conf.problemName);
    problem->GetMeshes().PrintInfo();
    auto solution = pdeSolver.Run(problem);
    pdeSolver.PrintValues(solution);
    return;
  }

  if (conf.modelName.find("DGTransport") != std::string::npos) {
    TransportPDESolver pdeSolver(conf);
    auto problem = CreateTransportProblemShared(conf.problemName);
    problem->GetMeshes().PrintInfo();
    auto solution = pdeSolver.Run(problem);
    return;
  }

  if (conf.modelName.find("Reaction") != std::string::npos) {
    ReactionPDESolver pdeSolver(conf);
    auto problem = CreateReactionProblemShared(conf.problemName);
    problem->GetMeshes().PrintInfo();
    auto solution = pdeSolver.Run(problem);
    return;
  }
  
  if (conf.modelName.find("Gelfand") != std::string::npos) {
    GelfandPDESolver pdeSolver(conf);
    auto problem = CreateGelfandProblemShared(conf.problemName);
    problem->GetMeshes().PrintInfo();
    auto solution = pdeSolver.Run(problem);
    pdeSolver.PrintValues(solution);
    return;
  }

  if (conf.modelName.find("Plate") != std::string::npos) {
    PlatePDESolver pdeSolver(conf);
    auto problem = CreatePlateProblemShared(conf.problemName);
    problem->GetMeshes().PrintInfo();
    auto solution = pdeSolver.Run(problem);
    pdeSolver.PrintValues(solution);
    return;
  }

  
#ifdef USE_SPACETIME
  if (conf.modelName.find("STAcoustic") != std::string::npos) {
    auto problem = CreateAcousticProblemShared(conf.problemName);
    problem->CreateMeshes(MeshesCreator(conf.mesh).WithPLevel(conf.pLevel).WithLevel(conf.level));
    problem->GetMeshes().PrintInfo();
    if (conf.modelName == "STAcousticPDESolver") {
      STAcousticPDESolver pdeSolver(conf);
      auto solution = pdeSolver.Run(problem);
      return;
    }
    if (conf.modelName == "STAcoustic") {
      STAcousticMain pdeSolver(conf);
      auto solution = pdeSolver.Run(problem);
      return;
    }
    if (conf.modelName == "STAcousticGL") {
      STAcousticMainGL pdeSolver(conf);
      auto solution = pdeSolver.Run(problem);
      return;
    }
    if (conf.modelName == "STAcousticMF") {
      STAcousticMainMF pdeSolver(conf);
      auto solution = pdeSolver.Run(problem);
      return;
    }
    if (conf.modelName == "STAcousticPG") {
      STAcousticMainPG pdeSolver(conf);
      auto solution = pdeSolver.Run(problem);
      return;
    }
  }
  if (conf.modelName == "STTransport") {
    auto problem = CreateTransportProblemShared(conf.problemName);
    STTransportMain pdeSolver(conf);
    auto solution = pdeSolver.Run(problem);
    return;
  }
#endif
  if (conf.modelName.find("VectorValued") != std::string::npos) {
    VectorValuedMain pdeSolver(conf);
    auto problem = CreateVectorValuedProblemShared(conf.problemName);
    auto solution = pdeSolver.Run(problem);
    return;
  }

  if (conf.modelName.find("Acoustic") != std::string::npos) {
    AcousticPDESolver pdeSolver(conf);
    auto problem = CreateAcousticProblemShared(conf.problemName);
    problem->GetMeshes().PrintInfo();
    auto solution = pdeSolver.Run(problem);
    return;
  }

  Exit(conf.modelName + " not found")
}
