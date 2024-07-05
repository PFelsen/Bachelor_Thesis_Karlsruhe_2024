#include "TestEnvironment.hpp"
#include "ProposalGenerators.hpp"
#include "RVector.hpp"

void TestProposalBehavior(const ProposalConfig& conf, double expectedProposedState) {
    RVector proposal(1);
    std::unique_ptr<ProposalGenerator> generator = CreateProposalGeneratorUnique(conf, proposal);
    EXPECT_EQ(generator->state[0], 1);

    generator->Propose();
    RVector proposedState = generator->GetProposal();
    EXPECT_EQ(proposedState[0], expectedProposedState);
}

TEST(ProposalCorrect, RandomWalk) {
    ProposalConfig conf = ProposalConfig().WithProposalType("RandomWalk")
      .WithStepLength(1)
      .WithPriorConfig(PriorSamplerConfig().WithGenerator("Test"));
    TestProposalBehavior(conf, 2);
}
TEST(ProposalCorrect, pCN) {
    ProposalConfig conf = ProposalConfig().WithProposalType("pCN")
      .WithStepLength(0.5)
      .WithPriorConfig(PriorSamplerConfig().WithGenerator("Test"));
    TestProposalBehavior(conf, std::sqrt(0.75) + 0.5);
}
TEST(ProposalCorrect, PriorDraw) {
    ProposalConfig conf = ProposalConfig().WithProposalType("PriorDraw")
      .WithStepLength(1)
      .WithPriorConfig(PriorSamplerConfig().WithGenerator("Test"));
    TestProposalBehavior(conf, 1);
}

TEST(AcceptanceProbabilityTest, RandomWalkNormal) {
  ProposalConfig conf = ProposalConfig()
    .WithProposalType("RandomWalk")
    .WithPriorConfig(PriorSamplerConfig().WithGenerator("Normal").WithVariance(1));
  RVector proposal(1);
  std::unique_ptr<ProposalGenerator> generator = CreateProposalGeneratorUnique(conf, proposal);
  generator->GetProposal()[0] = 1;
  generator->GetState()[0] = 0;
  EXPECT_DOUBLE_EQ(generator->AcceptanceProbability(0, 1), std::exp(-1.5));
}
TEST(AcceptanceProbabilityTest, RandomWalkUniform) {
  ProposalConfig conf = ProposalConfig()
    .WithProposalType("RandomWalk")
    .WithPriorConfig(PriorSamplerConfig().WithGenerator("Uniform").WithVariance(1));
  RVector proposal(1);
  std::unique_ptr<ProposalGenerator> generator = CreateProposalGeneratorUnique(conf, proposal);
  generator->GetProposal()[0] = 0.5;
  generator->GetState()[0] = 0;
  EXPECT_DOUBLE_EQ(generator->AcceptanceProbability(0, 1), std::exp(-1));
  generator->GetProposal()[0] = 1.02;
  generator->GetState()[0] = 0;
  EXPECT_DOUBLE_EQ(generator->AcceptanceProbability(0, 1), 0);
}
TEST(AcceptanceProbabilityTest, pCN) {
  ProposalConfig conf = ProposalConfig().WithProposalType("pCN");
  RVector proposal(1);
  std::unique_ptr<ProposalGenerator> generator = CreateProposalGeneratorUnique(conf, proposal);
  EXPECT_DOUBLE_EQ(generator->AcceptanceProbability(0, 1), std::exp(-1));
  EXPECT_DOUBLE_EQ(generator->AcceptanceProbability(1.02, 1), 1);
}
TEST(AcceptanceProbabilityTest, PriorDraw) {
  ProposalConfig conf = ProposalConfig().WithProposalType("PriorDraw");
  RVector proposal(1);
  std::unique_ptr<ProposalGenerator> generator = CreateProposalGeneratorUnique(conf, proposal);
  EXPECT_DOUBLE_EQ(generator->AcceptanceProbability(0, 1), std::exp(-1));
  EXPECT_DOUBLE_EQ(generator->AcceptanceProbability(1.02, 1), 1);
}
TEST(AcceptanceProbabilityTest, GaussianProbability) {
  NormalSampler sampler(PriorSamplerConfig().WithMean(1).WithVariance(2));
  EXPECT_DOUBLE_EQ(sampler.Prob(RVector(2.0, 1)), std::exp(-0.25));
}

TEST(AcceptanceProbabilityTest, AdaptiveProposal) {
  ProposalConfig conf = ProposalConfig()
    .WithProposalType("RandomWalk")
    .WithStepLength(0.5)
    .WithPriorConfig(PriorSamplerConfig().WithGenerator("Uniform").WithVariance(1));
  RVector proposal(1);
  std::unique_ptr<ProposalGenerator> generator = CreateProposalGeneratorUnique(conf, proposal);
  for (int i = 0; i < 8; i++) {
    generator->Propose();
    EXPECT_DOUBLE_EQ(generator->GetAcceptanceProbability(), 0);
    generator->AdaptStepSize();
  }
  EXPECT_DOUBLE_EQ(generator->GetStepSize(), 0.5*0.7);
  for (int i = 0; i < 8; i++) {
    generator->Propose();
    generator->Accept();
    EXPECT_DOUBLE_EQ(generator->GetAcceptanceProbability(), 1);
    generator->AdaptStepSize();
  }
  EXPECT_DOUBLE_EQ(generator->GetStepSize(), 0.5*0.7*1.3);
}

int main(int argc, char **argv) {
  return MppTest(
      MppTestBuilder(argc, argv).
          WithConfigEntry("ConfigVerbose", "0").
          WithConfigEntry("NewtonVerbose", "0").
          WithConfigEntry("LinearVerbose", "0").
          WithConfigEntry("MeshVerbose", "0").
          WithConfigEntry("AssembleVerbose", "0").
          WithRandomInitialized().
          WithPPM()
  ).RUN_ALL_MPP_TESTS();
}