#ifndef PROPOSAL_GENERATOR_H
#define PROPOSAL_GENERATOR_H

#include "RVector.hpp"
#include "SamplerSelection.hpp"


class ProposalGenerator {
protected:
  int n_proposed = 0;
  int n_accepted = 0;
  double step_length = 0.5;
  ProposalGenerator(ProposalConfig conf, RVector &randomInput);
public:
  std::unique_ptr<Sampler> prior_sampler;
  RVector state;
  RVector &proposal;
  virtual void Propose() = 0;
  void Accept();
  RVector &GetProposal();
  RVector &GetState();
  double GetStepSize();
  double GetAcceptanceProbability();
  void SetState(RVector new_state);
  void AdaptStepSize();
  virtual std::string Name() const = 0;
  virtual double AcceptanceProbability(double old_potential, double new_potential,
                                       double temperature = 1) = 0;
};

class RandomWalkGenerator : public ProposalGenerator {
public:
  RandomWalkGenerator(ProposalConfig conf, RVector &expansion_parameters);
  void Propose();

  std::string Name() const override { return "Random Walk"; }

  double AcceptanceProbability(double old_potential, double new_potential, double temperature = 1);
};

class PCN_Generator : public ProposalGenerator {
public:
  PCN_Generator(ProposalConfig conf, RVector &expansion_parameters);
  void Propose();

  std::string Name() const override { return "Preconditioned Crank Nicolson"; }

  double AcceptanceProbability(double old_potential, double new_potential, double temperature = 1);
};

class PriorDrawGenerator : public ProposalGenerator {
public:
  PriorDrawGenerator(ProposalConfig conf, RVector &expansion_parameters) :
      ProposalGenerator(conf, expansion_parameters) {}

  void Propose();

  std::string Name() const override { return "Prior Draw"; }

  double AcceptanceProbability(double old_potential, double new_potential, double temperature = 1);
};

ProposalGenerator *CreateProposalGenerator(ProposalConfig conf, RVector &randomInput);
std::unique_ptr<ProposalGenerator> CreateProposalGeneratorUnique(ProposalConfig conf,
                                                                 RVector &randomInput);

#endif // PROPOSAL_GENERATOR_H
