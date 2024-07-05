#ifndef SPACETIME_STMAIN_HPP
#define SPACETIME_STMAIN_HPP

#include "SpaceTimePlotting.hpp"
#include "Results.hpp"
#include "MeshesCreator.hpp"
#include "ErrorEstimator.hpp"

class STMainBuilder {
public:
  std::string meshName = "";
  std::string problemName = "";
  std::string modelName = "";
  std::string errorEstimator = "";
  std::vector<std::string> callbacks;
  int space_deg = -1;
  int time_deg = -1;

  int level = -1;
  int plevel = -1;

  int refinement_steps = -1;

  STMainBuilder &WithMesh(const std::string &mesh) {
    this->meshName = mesh;
    return *this;
  }

  STMainBuilder &WithProblem(const std::string &problem) {
    this->problemName = problem;
    return *this;
  }

  STMainBuilder &WithModel(const std::string &model) {
    this->modelName = model;
    return *this;
  }

  STMainBuilder &WithLevel(int level) {
    this->level = level;
    return *this;
  }

  STMainBuilder &WithPLevel(int plevel) {
    this->plevel = plevel;
    return *this;
  }

  STMainBuilder &WithSpaceDegree(int sDeg) {
    this->space_deg = sDeg;
    return *this;
  }

  STMainBuilder &WithTimeDegree(int tDeg) {
    this->time_deg = tDeg;
    return *this;
  }

  STMainBuilder &WithRefinements(int refinements) {
    this->refinement_steps = refinements;
    return *this;
  }

  STMainBuilder &WithErrorEstimator(const std::string &errorEstimator) {
    this->errorEstimator = errorEstimator;
    return *this;
  }

  STMainBuilder &WithCallbacks(const std::vector<std::string> &callbacks) {
    this->callbacks = callbacks;
    return *this;
  }

  STMainBuilder &ReadConfig() {
    if (meshName.empty()) Config::Get("Mesh", meshName);
    if (problemName.empty()) Config::Get("Problem", problemName);
    if (modelName.empty()) Config::Get("Model", modelName);
    if (errorEstimator.empty()) Config::Get("ErrorEstimator", errorEstimator);
    if (callbacks.empty()) Config::Get("Callbacks", callbacks);
    if (space_deg < 0) Config::Get("degree", space_deg);
    if (time_deg < 0) Config::Get("time_degree", time_deg);
    if (level < 0) Config::Get("level", level);
    if (plevel < 0) Config::Get("plevel", plevel);
    if (refinement_steps < 0) {
      Config::Get("refinement_steps", refinement_steps);
      if (refinement_steps < 0) refinement_steps = 0;
    }
    return *this;
  }
};

class STMethod {
protected:
  int verbose = -1;

  int plotting = -1;

  Results results;
  std::unique_ptr<Meshes> meshes;
  std::unique_ptr<STAssemble> assemble;
  std::unique_ptr<LinearSolver> solver;

  std::unique_ptr<Vector> oldSolution = nullptr;

  LevelPair levels;

  std::unique_ptr<ErrorEstimator> ee;

  std::unique_ptr<Preconditioner> createPC();

  void computeNorms(const Vector &U, const Matrix &M, const Vector &RHS);

  void plotParams(const Vector &U);

  void plotDegrees(const Vector& U);

  void plotSolution(const Vector& U);

  void plotVtu(const Vector &U);

  void checkForHighDegree(LevelPair levels);

  virtual void SetStartingValue(Matrix &M, Vector &U, Vector &RHS);

public:
  explicit STMethod(STMainBuilder builder);

  STMethod() : STMethod(STMainBuilder().ReadConfig()) {
    
  }

  virtual void run() = 0;

  Results GetResults();
};

#endif //SPACETIME_STMAIN_HPP
