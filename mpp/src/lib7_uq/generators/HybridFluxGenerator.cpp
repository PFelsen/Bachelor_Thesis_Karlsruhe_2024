#include "HybridEllipticAssemble.hpp"
#include "HybridFluxGenerator.hpp"
#include "EllipticPDESolver.hpp"


//void HybridFaceNormalFluxGenerator::InitGenerator(
//  std::shared_ptr<Meshes> meshes, int init) {
//  SampleGenerator<Scalar>::InitGenerator(meshes, init);
//  string problemName;
//  if (meshes->dim() == 1) problemName = "StochasticLaplace1D";
//  if (meshes->dim() == 2) problemName = "StochasticLaplace2D";
//
//  PDESolverConfig conf(
//      meshes->Level(), 1, "HybridElliptic", problemName
//  );
//
//  pdeSolver = new EllipticPDESolver(conf);
//}

void HybridFaceNormalFluxGenerator::drawSample(const SampleID &id) {
  mout.StartBlock(Name());
  vout(1) << id << endl;
  solutionFaceValues = std::make_unique<SampleSolution>(pdeSolver->GetSharedDisc(), id);
  solutionFaceFlux = std::make_unique<SampleSolution>(pdeSolver->GetSharedDisc(), id);
  vout(1) << "Run internal PDESolver" << endl;
  Solution solution = pdeSolver->Run(problem);
  vout(1) << "Set Normal Flux" << endl;
  pdeSolver->SetNormalFlux(solutionFaceValues->solution.vector, solutionFaceFlux->solution.vector);
  mout.EndBlock(verbose == 0);
}

Scalar HybridFaceNormalFluxGenerator::EvalSample(int face, const Cell &c) const {
  RTLagrangeElement elem(solutionFaceFlux->solution.vector, c);
  RTFaceElement faceElem(solutionFaceFlux->solution.vector, c, face);
  return (solutionFaceFlux->solution.vector)(elem[face], 0) * elem.Sign(face) / faceElem.Area();
}

HybridFaceNormalFluxGenerator::~HybridFaceNormalFluxGenerator() {
  delete pdeSolver;
}

VectorField HybridCellFluxGenerator::EvalSample(const Cell &c) const {
  RTLagrangeElement elem(generator.solutionFaceFlux->solution.vector, c);
  return elem.CellFlux(generator.solutionFaceFlux->solution.vector);
}

