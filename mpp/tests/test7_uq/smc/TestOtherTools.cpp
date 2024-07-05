#include "TestEnvironment.hpp"
#include "OtherTools.hpp"
#include "json.hpp"

using json = nlohmann::json;

// Todo ise filesystem header, this test is dangerous in parallel,
//  PPM guards added, this should be removed again and written as ctest!

TEST(OtherTools, MeanInColumn) {
  if (!PPM->Master(0)) return;
  std::vector<std::vector<double>> data = {{0, 1}, {2, 3}};
  EXPECT_DOUBLE_EQ(mean_in_column(data)[0], 1);
  EXPECT_DOUBLE_EQ(mean_in_column(data)[1], 2);
}

TEST(OtherTools, VarianceInColumn) {
  if (!PPM->Master(0)) return;
  std::vector<std::vector<double>> data = {{0, 1}, {2, 3}};
  std::vector<double> mean = {1, 2};
  EXPECT_DOUBLE_EQ(variance_in_column(data, mean)[0], 2);
  EXPECT_DOUBLE_EQ(variance_in_column(data, mean)[1], 2);
}

TEST(TestJSON, WriteJson)  {
  if (!PPM->Master(0)) return;
  std::vector<double> reference_distribution = {1, 2, 3, 4, 5};
  json jsonData;
  jsonData["reference_distribution"] = reference_distribution;
  std::ofstream outputFile;
  outputFile.open("test.json");
  if (!outputFile.is_open()) {
    std::invalid_argument("Json not found");
  }
  outputFile << jsonData.dump(4);
  outputFile.close();
}

TEST(TestJSON, ReadJson) {
  if (!PPM->Master(0)) return;
  std::vector<double> reference_distribution;
  std::ifstream inputFile;
  inputFile.open("test.json");
  if (!inputFile.is_open()) {
    std::invalid_argument("Json not found");
  }
  json jsonData;
  inputFile >> jsonData;
  if (jsonData.contains("reference_distribution")) {
      reference_distribution = jsonData["reference_distribution"].get<std::vector<double>>();
  }
  inputFile.close();

  ASSERT_EQ(reference_distribution.size(), 5);
  ASSERT_EQ(reference_distribution[0], 1);
  ASSERT_EQ(reference_distribution[1], 2);
  ASSERT_EQ(reference_distribution[2], 3);
  ASSERT_EQ(reference_distribution[3], 4);
  ASSERT_EQ(reference_distribution[4], 5);
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