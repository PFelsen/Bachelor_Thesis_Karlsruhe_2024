#ifndef PDESOLVER_HPP
#define PDESOLVER_HPP

#include "IProblem.hpp"
#include "MemoryLogger.hpp"
#include "SeismogramData.hpp"
#include "IDiscretization.hpp"
#include "Vector.hpp"

#include <typeinfo>
#include <utility>

#if BUILD_UQ
// Todo remove
struct PriorSamplerConfig{
  std::string generator = "Normal";
  double mean = 0;
  double variance = 1;

  PriorSamplerConfig();
  PriorSamplerConfig WithGenerator(std::string generator);
  PriorSamplerConfig WithMean(double mean);
  PriorSamplerConfig WithVariance(double stdv);
};

// Todo remove
struct ProposalConfig {
  std::string proposal_type = "RandomWalk";
  double step_length = 0.2;
  PriorSamplerConfig prior_sampler = PriorSamplerConfig();
  PriorSamplerConfig initial_sampler = PriorSamplerConfig().WithVariance(-1);
  string field_name = "Simple2DKL";
  bool adaptive = false;

  ProposalConfig();

  ProposalConfig WithStepLength(double step_length);
  ProposalConfig WithProposalType(std::string type);
  ProposalConfig WithPriorConfig(PriorSamplerConfig conf);
  ProposalConfig WithInitialSampler(PriorSamplerConfig conf);
  ProposalConfig WithRandomField(std::string name);
  ProposalConfig WithAdaptive();
};
#endif

// Todo: This one should be leaner
struct PDESolverConfig {
  double T{};

  double dt{};

  int level{};

  int pLevel{};

  int degree{};

  int rkorder{};

  int minLevel{};

  int timeDegree{};

  std::string mesh{};

  std::string modelName{};

  std::string problemName{};

  std::string preconditioner = "SuperLU";

  bool dualPrimal = false;

  std::string distribution = "none";

  std::string linearSolver = "GMRES";

  bool usePrevious = false;

  PDESolverConfig(int level, int degree, std::string modelName, std::string problemName) :
      level(level), pLevel(level - 1), degree(degree), modelName(std::move(modelName)),
      problemName(std::move(problemName)) {}

  PDESolverConfig(const std::string &modelName, const std::string &problemName) :
      PDESolverConfig(0, 0, modelName, problemName) {}

  PDESolverConfig() {
    if (!Config::IsInitialized()) return;

    Config::Get("Level", level);
    Config::Get("level", level);
    pLevel = max(level - 1, 0);
    Config::Get("plevel", pLevel);

    Config::Get("degree", degree);

    Config::Get("Model", modelName);
    Config::Get("MinLevel", minLevel);
    Config::Get("Problem", problemName);

    Config::Get("Preconditioner", preconditioner);
    Config::Get("time_degree", timeDegree);
    Config::Get("Mesh", mesh);
    Config::Get("Distribution", distribution);

    Config::Get("T", T);
    Config::Get("dt", dt);

    Config::Get("rkorder", rkorder);

    Config::Get("DualPrimal", dualPrimal);
    Config::Get("LinearSolver", linearSolver);

    Config::Get("UsePrevious", usePrevious);
  }

  PDESolverConfig WithLevel(int _level) {
    pLevel = (preconditioner == "Multigrid") ? minLevel : _level - 1;
    this->level = _level;
    return *this;
  }

  PDESolverConfig WithPLevel(int _plevel) {
    this->pLevel = _plevel;
    return *this;
  }

  PDESolverConfig WithDegree(const DegreePair &deg) {
    this->degree = deg.space;
    this->timeDegree = deg.time;
    return *this;
  }

  PDESolverConfig WithDegree(int _degree) {
    this->degree = _degree;
    return *this;
  }

  PDESolverConfig WithMinLevel(int minimumLevel) {
    minLevel = minimumLevel;
    return *this;
  }

  PDESolverConfig WithModel(std::string model) {
    modelName = std::move(model);
    return *this;
  }

  PDESolverConfig WithProblem(std::string problem) {
    problemName = std::move(problem);
    return *this;
  }

  PDESolverConfig WithMesh(std::string meshName) {
    mesh = meshName;
    return *this;
  }

#if BUILD_UQ
  // Todo remove again
  ProposalConfig proposal_config = ProposalConfig();

  std::vector<std::pair<std::string, std::vector<double>>> observations = {{"Point", {0.5, 0.5}}};

  PDESolverConfig WithObservations(std::vector<std::pair<std::string, std::vector<double>>> obs) {
    this->observations = obs;
    return *this;
  }

  PDESolverConfig WithProposalConfig(ProposalConfig conf) {
    this->proposal_config = conf;
    return *this;
  }
#endif
};

typedef std::map<std::string, double> ValueMap;
typedef std::map<std::string, std::vector<double>> MultiValueMap;

struct Solution {
  explicit Solution(std::shared_ptr<const IDiscretization> disc)
      : vector(Vector(0.0, std::move(disc))) {}

  explicit Solution(std::shared_ptr<const IDiscretization> disc, LevelPair meshIndex) :
      vector(Vector(0.0, std::move(disc), meshIndex)) {}

  explicit Solution(const Vector &vector) : vector(Vector(vector)) {}

  double computingTime{};

  double initTime{};

  std::string name{};

  ValueMap values{};

  // TODO: find a better way to store time-dependent data or maybe this is fine?
  MultiValueMap multiValues{};

  bool converged{};

  Vector vector;

  std::shared_ptr<SeismogramData> seismograms;
};

template <typename Problem, typename SolutionType>
class GenericPDESolver {
protected:
  int verbose = 1;

  int plotting = 0;

  PDESolverConfig conf;

  virtual void run(SolutionType &solution) const {};

  virtual void runAdjoint(SolutionType &solution) const {};

  virtual void plotVtu(SolutionType &solution) const {};

  virtual void computeValues(SolutionType &solution) const {};

  virtual void computeAdjointValues(SolutionType &solution) const {};

  virtual void createAssemble(std::shared_ptr<Problem> problem) {};

  virtual void createAdjointAssemble(std::shared_ptr<Problem> problem) {};

  virtual void printValues(SolutionType &solution) const {};

public:
  explicit GenericPDESolver(PDESolverConfig conf) : conf(std::move(conf)) {
    Config::Get("PDESolverVerbose", verbose);
    Config::Get("VtuPlot", plotting);
  }

  // Todo remove this one, or discuss other usage.
  //  Maybe Problem should be solved on all Meshes passed with Problem, and then gives back a Solution Collection?
  SolutionType Run(std::shared_ptr<Problem> problem, int commSplit = 0) {
    LevelPair level;
    
    if (problem->IsMeshInitialized()){
      level = problem->GetMeshes().FineLevel().WithCommSplit(commSplit);
    } else {
      level = {conf.level, conf.level, 0, commSplit};
    }
    return Run(problem, level);
  }

  // Todo, only use this one
  // Todo rename LevelPair to MeshIndex
  // Todo, make problem const&
  // Maybe: Template based solution
  SolutionType Run(std::shared_ptr<Problem> problem, LevelPair meshIndex) {
    mout.StartBlock(Name());
    vout(1) << "Solve Problem " << problem->Name() << endl;
    double start = MPI_Wtime();
    createAssemble(problem);
    SolutionType solution(GetSharedDisc(), meshIndex);
    solution.initTime = (MPI_Wtime() - start);
    run(solution);
    if (plotting) plotVtu(solution);
    computeValues(solution);
    PrintValues(solution);
    solution.computingTime = (MPI_Wtime() - start);
    mout.EndBlock(verbose == 0);
    return solution;
  }

  SolutionType RunAdjoint(std::shared_ptr<Problem> problem, LevelPair meshIndex) {
    mout.StartBlock(Name());
    vout(1) << "Solve AdjointProblem " << problem->Name() << endl;
    double start = MPI_Wtime();
    createAdjointAssemble(std::move(problem));
    SolutionType solution(GetSharedDisc(), meshIndex);
    solution.initTime = (MPI_Wtime() - start);
    runAdjoint(solution);
    if (plotting) plotVtu(solution);
    computeAdjointValues(solution);
    solution.computingTime = (MPI_Wtime() - start);
    mout.EndBlock(verbose == 0);
    return solution;
  }

  virtual std::string Name() const { return "PDESolver"; };

  virtual ~GenericPDESolver() = default;

  // TODO: Delete this ComputeValues and use computeValues(Solution &)
  virtual std::map<std::string, double> ComputeValues(const Vector &u) const { return {}; };

  // Todo merge this with virtual void plotVtu(Solution &solution) const {}
  virtual void PlotVtu(const Vector &u) const {};

  // Todo Call this one if verbose is high enough and protect -> Prints ValuesMap or Maybe introduce a PrintInfo for Solution
  virtual void PrintValues(const SolutionType &solution) {};

  const PDESolverConfig &GetConfig() const { return conf; }

  virtual std::shared_ptr<const IDiscretization> GetSharedDisc() const = 0;
};

template<typename Problem>
using PDESolver = GenericPDESolver<Problem, Solution>;

#endif //PDESOLVER_HPP
