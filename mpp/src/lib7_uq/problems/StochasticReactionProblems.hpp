#ifndef STOCHASTICREACTIONPROBLEMS_HPP
#define STOCHASTICREACTIONPROBLEMS_HPP

#include "IStochasticProblem.hpp"
#include "ReactionProblems.hpp"


class IStochasticReactionProblem : public IStochasticProblem, public IReactionProblem {};

IStochasticReactionProblem *
CreateStochasticReactionProblem(const std::string &problemName);

std::unique_ptr<IStochasticReactionProblem>
CreateStochasticReactionProblemUnique(const std::string &problemName);

std::shared_ptr<IStochasticReactionProblem>
CreateStochasticReactionProblemShared(const std::string &problemName);

#endif //STOCHASTICREACTIONPROBLEMS_HPP
