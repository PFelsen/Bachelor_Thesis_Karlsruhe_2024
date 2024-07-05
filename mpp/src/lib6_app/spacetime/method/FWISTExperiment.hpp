#ifndef MPP_FWISTEXPERIMENT_HPP
#define MPP_FWISTEXPERIMENT_HPP

#include "PDESolver.hpp"
#include "STAcousticMain.hpp"

class FWISTExperiment {
private:
  std::unique_ptr<STAcousticMain> pdeSolver = nullptr;

  std::shared_ptr<Meshes> meshes;

  size_t ITERATIONS = 1;
  double ALPHA = 0.000001;

public:
  FWISTExperiment(PDESolverConfig conf) : pdeSolver(std::make_unique<STAcousticMain>(conf)) {
    Config::Get("FWI_ITERATIONS", ITERATIONS);
    Config::Get("FWI_ALPHA", ALPHA);

    meshes = MeshesCreator("FWIExampleGeometry").WithPLevel(conf.pLevel).WithLevel(conf.level).CreateShared();

  }

  void Method() {
    mout.StartBlock("FWI");
    std::shared_ptr<AcousticProblem> problem = CreateAcousticProblemShared(pdeSolver->GetConfig().problemName, meshes);

    Solution solution = pdeSolver->Run(problem);
    SeismogramData &observations = *solution.seismograms;

    pdeSolver->clearOldSolution();

    auto materialDisc = std::make_shared<STDiscretizationT_DGDG<>>(*meshes, DegreePair{0, 0}, 2, false);
    Vector material(1.0, materialDisc);


    for (int i = 0; i < ITERATIONS; i++) {
      FlexibleAcousticProblemData builder{.meshes = meshes,
                                          .problem = problem,
                                          .kappa = [&material](const Point &cellMidPoint){
                                                              return material(cellMidPoint, 0); },
                                          .rho = [&material](const Point &cellMidPoint){
                                                              return material(cellMidPoint, 1); },
                                          .FPressure = [&problem](double t, const Cell &c, const Point &x){
                                                              return problem->F(t, c, x, COMPONENT::P0); },
                                          .name = "ForwardProblem with current material and same source as " + problem->Name()};

      auto forwardProblem = std::make_shared<FlexibleAcousticProblem>(builder);

      Solution forwardSolution = pdeSolver->Run(forwardProblem);

      SeismogramData differenceObservations((*forwardSolution.seismograms) - observations);
      
      double norm = differenceObservations.GetData().normAccumulated();

      FlexibleAcousticProblemData adjointBuilder{.meshes = meshes,
                                                 .problem = problem,
                                                 .kappa = [&material](const Point &cellMidPoint){
                                                                  return material(cellMidPoint, 0); },
                                                 .rho = [&material](const Point &cellMidPoint){
                                                                  return material(cellMidPoint, 1); },
                                                 .sourceData = std::make_shared<SeismogramData>(differenceObservations),
                                                 .name = "AdjointProblem with current material and source is difference of seismograms"};

      auto adjointProblem = std::make_shared<FlexibleAcousticProblem>(adjointBuilder);

      Solution backwardSolution = pdeSolver->RunAdjoint(adjointProblem, material.Level());

      Vector materialUpdate = pdeSolver->CalculateMaterialUpdate(material, forwardSolution.vector, backwardSolution.vector);

      material -= ALPHA * materialUpdate;

      for(cell c=material.cells(); c!= material.cells_end(); c++){
        material(c(), 1) = 1.0;
        if (c()[1] < 0.2) material(c(), 0) = 1.0;
      }
      material.MakeAdditive();
      material.Accumulate();

      double min2 = 1e100;
      double max2 = -1e100;
      for(const double &d : material.GetData().Data()){
        min2 = std::min(min2, d);
        max2 = std::max(max2, d);
      }
      mout.PrintIteration(1,
                          PrintIterEntry<int>{"n", i, 0, 1},
                          PrintIterEntry<double>{"seismogram misfit", norm, 0, 1},
                          PrintIterEntry<double>{"min(material)", PPM->Min(min2), 0, 1},
                          PrintIterEntry<double>{"max(material)", PPM->Max(max2), 0, 1});
    }
  }
};

#endif