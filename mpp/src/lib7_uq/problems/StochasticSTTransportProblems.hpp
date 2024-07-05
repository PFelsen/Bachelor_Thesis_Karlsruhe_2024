#ifndef MPP_STOCHASTICSTTRANSPORTPROBLEMS_HPP
#define MPP_STOCHASTICSTTRANSPORTPROBLEMS_HPP

#include "TransportProblems.hpp"
#include "IStochasticProblem.hpp"


class IStochasticSTTransportProblem : public IStochasticProblem, public ITransportProblem {
public:
  IStochasticSTTransportProblem(int dim) : ITransportProblem() {}
};

IStochasticSTTransportProblem *
CreateTransportStochasticProblem(const string &problemName);

std::unique_ptr<IStochasticSTTransportProblem>
CreateTransportStochasticProblemUnique(const std::string &problemName);

std::shared_ptr<IStochasticSTTransportProblem>
CreateTransportStochasticProblemShared(const std::string &problemName);

#endif //STOCHASTICSTACOUSTICPROBLEMS_HPP
