#ifndef STOCHASTICSTACOUSTICPROBLEMS_HPP
#define STOCHASTICSTACOUSTICPROBLEMS_HPP

#include "AcousticProblems.hpp"
#include "IStochasticProblem.hpp"


class IStochasticSTAcousticProblem : public IStochasticProblem, public AcousticProblem {
public:
  IStochasticSTAcousticProblem(int dim, int nL) : AcousticProblem(dim, nL) {}
};

IStochasticSTAcousticProblem *
CreateAcousticStochasticProblem(const string &problemName);

std::unique_ptr<IStochasticSTAcousticProblem>
CreateViscoAcousticStochasticProblemUnique(const std::string &problemName);

std::shared_ptr<IStochasticSTAcousticProblem>
CreateViscoAcousticStochasticProblemShared(const std::string &problemName);

#endif //STOCHASTICSTACOUSTICPROBLEMS_HPP
