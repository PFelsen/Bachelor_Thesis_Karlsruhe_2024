#include "STSingleRunMethod.hpp"
#include "Plotting.hpp"

STSingleRunMethod::STSingleRunMethod(STMainBuilder builder) : STMethod(builder) {}

void STSingleRunMethod::run() {
  Vector U(0.0, assemble->GetDisc(), levels);
  U.PrintInfo();
  U.GetMatrixGraph().PrintInfo();
  U.GetMatrixGraph().PrintMatrixMemoryInfo();
  U.GetMesh().PrintInfo();

  // Todo -> All this stuff should happen in LinearSolver
  U.SetAccumulateFlag(true);
  Matrix B(U);
  Vector RHS(0.0, U);
  assemble->System(B, RHS);

  //mout << B << endl;

  SetStartingValue(B, U, RHS);
  (*solver)(B, true);
  RHS.SetAccumulateFlag(false);
  U += (*solver) * RHS;

  printVTK(U, *assemble);
  plotVtu(U);
  computeNorms(U, B, RHS);

  results.PrintInfo();
}
