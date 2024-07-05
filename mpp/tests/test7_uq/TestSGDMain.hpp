#ifndef TESTSGDMAIN_H
#define TESTSGDMAIN_H

#include <map>
#include <string>

using ConfigMap = std::map<std::string, std::string>;

const ConfigMap DefaultSGDConfigMap = {{"Protocol", "OCLaplace2D"},
                                       {"Model", "LagrangeElliptic"},
                                       {"degree", "1"},
                                       {"Quantity", "L2Error"},
                                       {"maxSteps", "100"},
                                       {"level", "4"},
                                       {"targetCost", "0.0"},
                                       {"stepsizeRule", "constant"},
                                       {"alpha", "1.0"},
                                       {"secondaryQOI", "value"},
                                       {"SGDVerbose", "1"},
                                       {"est_sample_size", "1"},
                                       {"M", "1"},
                                       {"overkill", "true"},
                                       {"theta", "55"},
                                       {"nu", "10"},
                                       {"gamma", "1000"},
                                       {"StochasticField", "LogNormal"},
                                       {"Mean", "0.0"},
                                       {"sigma", "1.0"},
                                       {"lambda", "[0.15, 0.15]"},
                                       {"smoothing", "1.0"},
                                       {"VtuPlot", "0"},
                                       {"SLEstimatorVerbose", "1"},
                                       {"MLMCVerbose", "1"},
                                       {"MainVerbose", "1"},
                                       {"MeshVerbose", "0"},
                                       {"MeshesVerbose", "0"},
                                       {"ConfigVerbose", "1"},
                                       {"LinearVerbose", "0"},
                                       {"NewtonVerbose", "0"},
                                       {"AssembleVerbose", "0"},
                                       {"PDESolverVerbose", "0"},
                                       {"GeneratorVerbose", "0"},
                                       {"ReferenceValue", "0.00551712"},
                                       {"ReferenceValueStepsize", "100"},
                                       {"u_a", "-1000"},
                                       {"u_b", "1000"},
                                       {"descentType", "SGD"},
                                       {"GradEst", "SGD"},
                                       {"breakpointCriteria", "iterations"}};

const ConfigMap Dummy = {};


#endif // TESTSGDMAIN_H
