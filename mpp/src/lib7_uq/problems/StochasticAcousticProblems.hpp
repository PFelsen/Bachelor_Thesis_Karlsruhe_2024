#ifndef STOCHASTICACOUSTICPROBLEMS_HPP
#define STOCHASTICACOUSTICPROBLEMS_HPP

#include "AcousticProblems.hpp"
#include "IStochasticProblem.hpp"

class IStochasticAcousticProblem : public IStochasticProblem, public AcousticProblem {
protected:
  IStochasticAcousticProblem(int dim, int numL) : AcousticProblem(dim, numL) {}
};

IStochasticAcousticProblem *
CreateStochasticAcousticProblem(const string &problemName);

std::unique_ptr<IStochasticAcousticProblem>
CreateStochasticAcousticProblemUnique(const string &problemName);

std::shared_ptr<IStochasticAcousticProblem>
CreateStochasticAcousticProblemShared(const std::string &problemName);

#endif //STOCHASTICACOUSTICPROBLEMS_HPP
