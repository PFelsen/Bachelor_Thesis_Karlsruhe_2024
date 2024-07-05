#ifndef MLUQ_TESTSEQUENTIALMONTECARLO_H
#define MLUQ_TESTSEQUENTIALMONTECARLO_H

using ConfigMap = std::map<std::string, std::string>;

const ConfigMap SMCConfig = {
  // ----- Metropolis Kernel -----
  {"SampleGenerator", "Normal"},
  {"ProposalGenerator", "RandomWalkMetropolis"},
  {"step_length", "0.2"},  
  
  // ----- Prior -----
  {"StochasticField", "Gauss"},
  {"Mean", "0.0"}, // TODO: should prior not always be mean free? Instead possibly: put measurement?
  {"sigma", "1.0"},
  {"lambda", "[0.15, 0.15]"},
  {"smoothing", "1.8"},

  // ----- Problem -----
  {"ProblemName", "DummyTest"},

  // TODO: Combine into Sequential Monte Carlo
  // ----- Particle -----
  {"NMarkovSteps", "20"},
  {"Measurement", "{1}"},
  {"Precision", "{1}"},
  // ----- Particle set -----
  {"NParticles", "100"},
  {"ResamplingScheme", "Multinomial"},
  // ----- Sequential Monte Carlo -----
  {"RelativeESSmin", "0.5"},
  {"Adaptive", "false"},


  // TODO: Not checked, just taken from template:

  // Verbose
  {"SLEstimatorVerbose", "1"}, // ??
  {"MainVerbose", "1"},
  {"MeshVerbose", "0"}, // ??
  {"ConfigVerbose", "1"},
  {"LinearVerbose", "0"}, // ??
  {"NewtonVerbose", "0"}, // ??
  {"PDESolverVerbose", "0"}, // ??
  {"GeneratorVerbose", "0"},

  // ----- Plotting -----
  {"VtuPlot", "0"}, // ??
};

#endif //MLUQ_TESTSEQUENCTIALMONTECARLO_H

