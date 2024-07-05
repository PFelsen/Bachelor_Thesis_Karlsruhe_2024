#include "TestEnvironment.hpp"
#include "SequentialMonteCarlo.hpp"
#include "OtherTools.hpp"

TEST(MCMCHyperParametersWithoutPDE, BurnIn) {
  const int num_iterations = 50;
  std::vector<std::vector<double>> no_burn_in_errors(num_iterations);
  std::vector<std::vector<double>> only_chain_errors(num_iterations);
  std::vector<std::vector<double>> increasing_burn_in100(num_iterations);
  std::vector<std::vector<double>> increasing_burn_in1000(num_iterations);

  std::vector<double> no_burn_in_lengths;
  std::vector<double> chain_length;
  std::vector<double> burn_in100;
  std::vector<double> burn_in1000;

  double precision = 10;
  double variance = 1000;
  double posterior_mean = (1/precision)/((1/precision) + 1/variance);

  for (int i = 0; i < num_iterations; ++i) {
    SequentialMonteCarlo mcmc{
      SequentialMonteCarloConfig()
        .WithVersion("MCMC")
        .WithNIntermediate(10000)
        .WithNStepsAtOnce(1)
        .WithParticleSetConfig(ParticleSetConfig()
          .WithNParticles(1)
          .WithParticleConfig(ParticleConfig()
            .WithMeasurement({1.0})
            .WithPrecision({precision})
            .WithPDEModelConfig(PDESolverConfig()
              .WithModel("Test")
              .WithProposalConfig(ProposalConfig()
                .WithProposalType("RandomWalk")
                .WithRandomField("Simple2DKL")
                .WithStepLength(0.02)
                .WithAdaptive()
                .WithPriorConfig(PriorSamplerConfig()
                  .WithGenerator("Normal")
                  .WithVariance(variance))
                .WithInitialSampler(PriorSamplerConfig()
                  .WithMean(0.5)
                  .WithVariance(0.5)
                  .WithGenerator("Normal")
                  )))))};
      mcmc.Method();
      std::vector<RVector> markov_chain = mcmc.target_dist.EvaluateParticle(0, 0);
      for( int m = 9; m < 1000; m += 10) {
        if (i == 0) no_burn_in_lengths.push_back(m);
        std::vector<RVector> data(markov_chain.begin(), markov_chain.begin() + m);
        RVector estimate = Mean(data);
        no_burn_in_errors[i].push_back( RelativeError1(estimate, RVector(posterior_mean, 1), RVector(1.0, 1)) );
      }
      for( int m = 1010; m < 10000; m += 10) {
        if (i == 0) chain_length.push_back(m-1000);
        std::vector<RVector> data(markov_chain.begin()+1000, markov_chain.begin() + m);
        RVector estimate = Mean(data);
        only_chain_errors[i].push_back( RelativeError1(estimate, RVector(posterior_mean, 1), RVector(1.0, 1)) );
      }
      for( int m_0 = 0; m_0 < 100; m_0 += 2) {
        if (i == 0) burn_in100.push_back(m_0);
        std::vector<RVector> data(markov_chain.begin() + m_0, markov_chain.begin() + 100);
        RVector estimate = Mean(data);
        increasing_burn_in100[i].push_back( RelativeError1(estimate, RVector(posterior_mean, 1), RVector(1.0, 1)) );
      }
      for( int m_0 = 0; m_0 < 1000; m_0 += 4) {
        if (i == 0) burn_in1000.push_back(m_0);
        std::vector<RVector> data(markov_chain.begin() + m_0, markov_chain.begin() + 1000);
        RVector estimate = Mean(data);
        increasing_burn_in1000[i].push_back( RelativeError1(estimate, RVector(1.0, 1), RVector(1.0, 1)) );
      }
  }

  std::vector<double> no_burn_in_errors_mean = mean_in_column(no_burn_in_errors);
  std::vector<double> no_burn_in_errors_var = variance_in_column(no_burn_in_errors, no_burn_in_errors_mean);
  std::vector<double> only_chain_errors_mean = mean_in_column(only_chain_errors);
  std::vector<double> only_chain_errors_var = variance_in_column(only_chain_errors, only_chain_errors_mean);
  std::vector<double> increasing_burn_in100_mean = mean_in_column(increasing_burn_in100);
  std::vector<double> increasing_burn_in100_var = variance_in_column(increasing_burn_in100, increasing_burn_in100_mean);
  std::vector<double> increasing_burn_in1000_mean = mean_in_column(increasing_burn_in1000);
  std::vector<double> increasing_burn_in1000_var = variance_in_column(increasing_burn_in1000, increasing_burn_in1000_mean);

  mout.PrintInfo("SMCHyperparameterData", 1, 
    // PrintInfoEntry("Study", "BurnInAndMarkovChainLength"),
    PrintInfoEntry("NoBurnInMean", no_burn_in_errors_mean), 
    PrintInfoEntry("NoBurnInVar", no_burn_in_errors_var), 
    PrintInfoEntry("NoBurnInChainLength", no_burn_in_lengths),
    PrintInfoEntry("OnlyChainErrorsMean", only_chain_errors_mean), 
    PrintInfoEntry("OnlyChainErrorsVar", only_chain_errors_var), 
    PrintInfoEntry("OnlyChainLength", chain_length),
    PrintInfoEntry("IncreasingBurnIn100Mean", increasing_burn_in100_mean), 
    PrintInfoEntry("IncreasingBurnIn100Var", increasing_burn_in100_var), 
    PrintInfoEntry("BurnIn100Lenth", burn_in100),
    PrintInfoEntry("IncreasingBurnIn1000Mean", increasing_burn_in1000_mean), 
    PrintInfoEntry("IncreasingBurnIn1000Var", increasing_burn_in1000_var),
    PrintInfoEntry("BurnIn1000Lenth", burn_in1000));
}


int main(int argc, char **argv) {
  return MppTest(
    MppTestBuilder(argc, argv).
      WithConfigEntry("ConfigVerbose", "1").
      WithConfigEntry("NewtonVerbose", "0").
      WithConfigEntry("LinearVerbose", "0").
      WithRandomInitialized().
      WithPPM()
  ).RUN_ALL_MPP_TESTS();
}
