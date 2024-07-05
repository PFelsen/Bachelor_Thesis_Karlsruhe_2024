#ifndef MLMC_STOCHASTICTRANSPORTPROBLEM_HPP
#define MLMC_STOCHASTICTRANSPORTPROBLEM_HPP

#include "IStochasticProblem.hpp"
#include "TransportProblems.hpp"
#include "TimeSeries.hpp"


class IStochasticTransportProblem : public IStochasticProblem, public ITransportProblem {};

IStochasticTransportProblem *
CreateStochasticTransportProblem(const string &problemName);

std::unique_ptr<IStochasticTransportProblem>
CreateStochasticTransportProblemUnique(const std::string &problemName);

std::shared_ptr<IStochasticTransportProblem>
CreateStochasticTransportProblemShared(const std::string &problemName);

#endif