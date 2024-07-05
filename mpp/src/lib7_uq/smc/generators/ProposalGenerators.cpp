#include "ProposalGenerators.hpp"

ProposalGenerator::ProposalGenerator(ProposalConfig conf, RVector &expansion_parameters) :
    proposal(expansion_parameters) {
  prior_sampler = CreateSamplerUnique(conf.prior_sampler);

  if (conf.initial_sampler.variance == -1) proposal = prior_sampler->DrawSample(proposal.size());
  else {
    std::unique_ptr<Sampler> initial_sampler = CreateSamplerUnique(conf.initial_sampler);
    proposal = initial_sampler->DrawSample(proposal.size());
  }
  state = proposal;
}

RVector &ProposalGenerator::GetProposal() { return proposal; }

RVector &ProposalGenerator::GetState() { return state; }

void ProposalGenerator::SetState(RVector new_state) { state = new_state; }

void ProposalGenerator::Accept() {
  state = proposal;
  n_accepted++;
}

void ProposalGenerator::AdaptStepSize() {
  if (n_proposed >= 8) {
    double acceptance_rate = double(n_accepted) / n_proposed;
    if (acceptance_rate > 0.3 && step_length * 1.3 < 1) {
      step_length *= 1.3;
      n_accepted = 0;
      n_proposed = 0;
    }
    if (acceptance_rate < 0.15) {
      step_length *= 0.7;
      n_accepted = 0;
      n_proposed = 0;
    }
  }
}

double ProposalGenerator::GetStepSize() { return step_length; }

double ProposalGenerator::GetAcceptanceProbability() {
  if (n_proposed > 0) return double(n_accepted) / n_proposed;
  return -1;
}

RandomWalkGenerator::RandomWalkGenerator(ProposalConfig conf, RVector &expansion_parameters) :
    ProposalGenerator(conf, expansion_parameters) {
  step_length = conf.step_length;
}

void RandomWalkGenerator::Propose() {
  RVector xi = prior_sampler->DrawSample(state.size()); // TODO COMM_SPLIT
  proposal = state + step_length * xi;
  n_proposed++;
}

double RandomWalkGenerator::AcceptanceProbability(double old_potential, double new_potential,
                                                  double temperature) {
  return std::min<double>(1, std::exp(temperature * (old_potential - new_potential))
                                 * prior_sampler->Prob(proposal) / prior_sampler->Prob(state));
}

PCN_Generator::PCN_Generator(ProposalConfig conf, RVector &randomInput) :
    ProposalGenerator(conf, randomInput) {
  step_length = conf.step_length;
}

void PCN_Generator::Propose() {
  RVector xi = prior_sampler->DrawSample(state.size()); // TODO COMM_SPLIT
  proposal = std::sqrt(1 - step_length * step_length) * state + step_length * xi;
  n_proposed++;
}

double PCN_Generator::AcceptanceProbability(double old_potential, double new_potential,
                                            double temperature) {
  return std::min<double>(1, std::exp(temperature * (old_potential - new_potential)));
}

void PriorDrawGenerator::Propose() {
  proposal = prior_sampler->DrawSample(state.size()); // TODO COMM_SPLIT
  n_proposed++;
}

double PriorDrawGenerator::AcceptanceProbability(double old_potential, double new_potential,
                                                 double temperature) {
  return std::min<double>(1, std::exp(temperature * (old_potential - new_potential)));
}

ProposalGenerator *CreateProposalGenerator(ProposalConfig conf, RVector &expansion_parameters) {
  if (conf.proposal_type == "RandomWalk")
    return new RandomWalkGenerator(conf, expansion_parameters);
  if (conf.proposal_type == "pCN") return new PCN_Generator(conf, expansion_parameters);
  if (conf.proposal_type == "PriorDraw") return new PriorDrawGenerator(conf, expansion_parameters);
  Exit(conf.proposal_type + " not found, available proposals are RandomWalk, pCN and PriorDraw.")
};

std::unique_ptr<ProposalGenerator> CreateProposalGeneratorUnique(ProposalConfig conf,
                                                                 RVector &expansion_parameters) {
  return std::unique_ptr<ProposalGenerator>(CreateProposalGenerator(conf, expansion_parameters));
};

ProposalConfig::ProposalConfig() {
  if (!Config::IsInitialized()) return;

  Config::Get("ProposalType", proposal_type);
  Config::Get("step_length", step_length);
  Config::Get("RandomField", field_name);
}

ProposalConfig ProposalConfig::WithStepLength(double step_length) {
  this->step_length = step_length;
  return *this;
}

ProposalConfig ProposalConfig::WithProposalType(std::string type) {
  this->proposal_type = type;
  return *this;
}

ProposalConfig ProposalConfig::WithPriorConfig(PriorSamplerConfig conf) {
  this->prior_sampler = conf;
  return *this;
}

ProposalConfig ProposalConfig::WithRandomField(std::string name) {
  this->field_name = name;
  return *this;
}

ProposalConfig ProposalConfig::WithInitialSampler(PriorSamplerConfig conf) {
  this->initial_sampler = conf;
  return *this;
}

ProposalConfig ProposalConfig::WithAdaptive() {
  this->adaptive = true;
  return *this;
}
