#include "STPGViscoAcousticAssemble.hpp"
#include "STPGViscoElasticAssemble.hpp"
#include "STDGElasticity.hpp"
#include "STDGViscoAcousticAssemble.hpp"
#include "STDGTransportAssemble.hpp"
#include "STGLViscoAcousticAssemble.hpp"
#include "STDGMatrixFreeViscoAcousticAssemble.hpp"

#include "SpaceTimePlotting.hpp"

double STAssemble::MhalfLuMinusF(const Vector &U, const std::string &filename) const {
  Vector v(0.0, U);
  for(cell c = U.cells(); c != U.cells_end(); c++){
    v(c(),0) = MminushalfLuMinusF(c, U);
  }

  //mout.PrintEntries(filename, 3, PrintInfoEntry("max", max(abs(v.Max()), abs(v.Min()))));

  //printVTK_Eta(v, *this, filename+"_L"+std::to_string(U.Level().space));
  return v.norm();
}

void STAssemble::PlotErrorEstimator(const Vector &Eta, std::string filename) {
  VtuPlot plot(filename, {.parallelPlotting=false,
                      .plotBackendCreator=std::make_unique<DefaultSpaceTimeBackend>});
  plot.AddData("Eta", Eta, 0);
  plot.PlotFile();
}


// TODO: use same name for Assemble and for model string
std::unique_ptr<STAssemble> CreateSTAssemble(
  const string &modelName, const Meshes &meshes, DegreePair degree, const string &probName) {

  if (modelName == "STPGViscoAcousticAssemble")
    return std::make_unique<STPGViscoAcousticAssemble>(meshes, degree, probName);
  else  if (modelName == "STDGViscoAcousticAssemble")
    return std::make_unique<STDGViscoAcousticAssemble>(meshes, degree, probName);
  else if (modelName == "STGLViscoAcousticAssemble" || modelName == "STGLGLViscoAcousticAssemble")
    return std::make_unique<STGLViscoAcousticAssemble>(meshes, degree, probName);
  else if (modelName == "STPGViscoElasticAssemble")
    return std::make_unique<STPGViscoElasticAssemble>(meshes, degree, probName);
  else if (modelName == "STDGMatrixFreeViscoAcousticAssemble")
    return std::make_unique<STDGMatrixFreeViscoAcousticAssemble>(meshes, degree, probName);
    // Todo: remove case below
  else if (modelName == "DGElasticity_adaptive")
    return std::make_unique<TDGElasticityAssemble>(meshes, degree, probName);
  Exit(modelName + " Model not found")
}
