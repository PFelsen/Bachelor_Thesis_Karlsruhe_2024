#ifndef _SPACETIME_H_
#define _SPACETIME_H_

#include "m++.hpp"
#include "MeshesCreator.hpp"
#include "SpaceTimePreconditioner.hpp"
#include "STPathManager.hpp"
#include "CyclicPreconditioner.hpp"
#include "Results.hpp"
#include "assemble/STAssemble.hpp"

/*
 * Todo Clean up this file ! ! !
 *
 */

class Errors {
  Errors() = delete;

public:

  std::map<std::string, vector<double>> errors;

  RVector GF_Value;
  RVector GF;
  vector<int> NCells;
  vector<int> NDoFs;
  RVector Int;
  RVector IntGN;
  RVector L1;
  RVector L2;
  RVector LInf;
  RVector L1_int;
  RVector L2_int;
  RVector LInf_int;
  RVector L2_proj;
  RVector DGSemiNormExactSolution;
  RVector DGSemiNormSolution;

  RVector ReconstructionErrorBound_L1;
  RVector ReconstructionErrorBound_L2;
  RVector ReconstructionErrorBound_E;

  RVector GN;
  RVector W;
  RVector V;

  Errors(int size)
    :
    GF_Value(0.0, size),
    GF(0.0, size),
    NCells(vector<int>(size)),
    NDoFs(vector<int>(size)),
    Int(0.0, size),
    IntGN(0.0, size),
    L1(0.0, size),
    L2(0.0, size),
    LInf(0.0, size),
    L1_int(0.0, size),
    L2_int(0.0, size),
    LInf_int(0.0, size),
    L2_proj(0.0, size),
    DGSemiNormExactSolution(0.0, size),
    DGSemiNormSolution(0.0, size),
    GN(0.0, size),
    W(0.0, size),
    ReconstructionErrorBound_L1(0.0, size),
    ReconstructionErrorBound_L2(0.0, size),
    ReconstructionErrorBound_E(0.0, size),
    V(0.0, size) {
  }

  vector<double> &operator[](const std::string &key) {
    auto elem = errors.find(key);
    if (elem == errors.end()) {
      elem = errors.insert({key, vector<double>()}).first;
    }
    return elem->second;
  }
};

class SpaceTimeMain {
private:
  int verbose = 0;
  bool exact;
  bool restore = 0;
  int sDeg = 0;
  int tDeg = 0;
  string probName = "SphericalWave2D";
  double theta = 1e-5;
  double theta_min = 1e-10;

  string functional_string = "none";
  int refinement = 0;
  string refine_by = "abs_value";
  int load_balancing = 0;
  Point roi_min;
  Point roi_max;
  double dual_red = 1e-2;
  double dual_eps = 1e-6;
  int vtkplot = 1;
  int coarse_high_poly = 0;
  bool useSparseMatrix = false;
  bool doMeasure = true;
  double timeScaleFactor = 0.5;
  int time_refinements = 0;
  std::string modelName;
  string pathChoice = "single";

  bool calc_conf_rec = false;

  std::shared_ptr<Meshes> meshes;
  std::unique_ptr<STAssemble> assemble;

  Results results;

public:
  Results& getResults(){
    return results;
  }

  void run();

  SpaceTimeMain() {
    Config::Get("SpacetimeVerbose", verbose);
    Config::Get("restoreAdaptiveRun", restore);
    Config::Get("degree", sDeg);
    Config::Get("time_degree", tDeg);
    Config::Get("Problem", probName);
    Config::Get("theta", theta);
    Config::Get("theta_min", theta_min);
    Config::Get("refine_by", refine_by);
    Config::Get("GoalFunctional", functional_string);
    Config::Get("refinement_steps", refinement);
    Config::Get("load_balancing", load_balancing);
    Config::Get("roi_min", roi_min);
    Config::Get("roi_max", roi_max);
    Config::Get("dual_red", dual_red);
    Config::Get("dual_eps", dual_eps);
    Config::Get("vtkplot", vtkplot);
    Config::Get("high_poly_deg_on_coarse_grids", coarse_high_poly);
    Config::Get("UseSparseMatrix", useSparseMatrix);
    Config::Get("doMeasure", doMeasure);
    Config::Get("ScaleFactor", timeScaleFactor);
    Config::Get("Scales", time_refinements);
    Config::Get("PathChoice", pathChoice);
    Config::Get("ConformingReconstruction", calc_conf_rec);
    Config::Get("Model", modelName);

    meshes =
        MeshesCreator().WithTimeRefinement(time_refinements).CreateShared();
  }

protected:

  void handleExactSolution(Vector &U,
                           Matrix &B,
                           const int run1,
                           Results &);

  void checkForRestore(const LevelPair level, int &run);
};


#endif
