#include <fstream>
#include <string>
#include <vector>

#include "Logging.hpp"
#include "TestEnvironment.hpp"

#ifdef BUILD_IA
#include <interval.hpp>

#include "IACInterval.hpp"
#include "IAInterval.hpp"
#endif

namespace IncludeEmptyLines {
enum __IncludeEmptyLines {
  Yes,
  No,
};
};

std::vector<std::string> ReadLines(int includeEmptyLines = IncludeEmptyLines::No) {
  std::ifstream infile(Config::GetBuildPath() + "/tests/test0_basic/logfile.log");
  std::string line;
  std::vector<std::string> lines;
  while (std::getline(infile, line)) {
    if (includeEmptyLines == IncludeEmptyLines::Yes || !line.empty()) lines.push_back(line);
  }
  return lines;
}

nlohmann::json ReadJson() {
  std::ifstream infile(Config::GetBuildPath() + "/tests/test0_basic/logfile.json");
  if (!infile) throw std::runtime_error("could not open File");
  nlohmann::json j;
  infile >> j;
  return j;
}

class DummyClass {
  std::string name = "";
public:
  DummyClass(std::string name) : name(name) {}

  template<typename S>
  friend LogTextStream<S> &operator<<(LogTextStream<S> &stream, const DummyClass &dc) {
    stream << "DummyClass with name <" << dc.name << ">" << std::endl;
    return stream;
  }
};

template<typename LOGGING_TYPE>
class TestLogging : public ::testing::Test {
protected:
  LOGGING_TYPE &getLogger() {
    if constexpr (std::is_same_v<LOGGING_TYPE, MasterLogging>) return mout;
    else if constexpr ((std::is_same_v<LOGGING_TYPE, LogTextStream<ProcessLogging>>)) return pout;
    else THROW("Not implemented")
  }

  std::string getIndent() {
    if constexpr (std::is_same_v<LOGGING_TYPE, MasterLogging>) return "   ";
    else if constexpr ((std::is_same_v<LOGGING_TYPE, LogTextStream<ProcessLogging>>)) return "0: ";
    else THROW("Not implemented")
  }

  template<typename T>
  void test(T value, std::string expectedLine) {
    expectedLine = getIndent() + expectedLine;
    getLogger() << value << std::endl;
    std::string logLine = ReadLines().back();
    EXPECT_EQ(logLine, expectedLine);
    // TODO: Check why this is not working!
    //    getLogger() << value << "\n" << std::flush;
    //    logLine = ReadLines().back();
    //    EXPECT_EQ(logLine, expectedLine);
  }

  void testDoubleBlock(double value, std::string expectedLine) {
    getLogger() << std::beginD;
    test<double>(value, expectedLine);
    getLogger() << std::endD;
  }
};

#define TEST_LOGGING(loggingClass)                                                                 \
                                                                                                   \
  TEST_F(loggingClass, ShortTest) {                                                                \
    test<short>(257, "257");                                                                       \
    test<short>(-4830, "-4830");                                                                   \
  }                                                                                                \
                                                                                                   \
  TEST_F(loggingClass, IntegerTest) {                                                              \
    test<int>(10542, "10542");                                                                     \
    test<int>(-34861, "-34861");                                                                   \
  }                                                                                                \
                                                                                                   \
  TEST_F(loggingClass, DoubleTest) {                                                               \
    test<double>(1.8754, "1.8754");                                                                \
    test<double>(-43.7462, "-43.7462");                                                            \
  }                                                                                                \
                                                                                                   \
  TEST_F(loggingClass, DoubleBlockTest) {                                                          \
    testDoubleBlock(0.0, "   0.00000000");                                                         \
    testDoubleBlock(-0.0, "  -0.00000000");                                                        \
    testDoubleBlock(infty, "    infty    ");                                                       \
    testDoubleBlock(-infty, "   -infty    ");                                                      \
    testDoubleBlock(1e-2, "   0.01000000");                                                        \
    testDoubleBlock(-1e-2, "  -0.01000000");                                                       \
    testDoubleBlock(1e2, " 100.00000000");                                                         \
    testDoubleBlock(-1e2, "-100.00000000");                                                        \
  }                                                                                                \
                                                                                                   \
  TEST_F(loggingClass, MultipleLineDoubleBlockTest) {                                              \
    getLogger() << std::beginD << 1.0045 << std::endl                                              \
                << -40.758 << std::endl                                                            \
                << 5.098 << std::endl                                                              \
                << std::endD;                                                                      \
    std::vector<std::string> lines = ReadLines();                                                  \
    std::string logLine = lines.back();                                                            \
    lines.pop_back();                                                                              \
    EXPECT_EQ(logLine, getIndent() + "   5.09800000");                                             \
    logLine = lines.back();                                                                        \
    lines.pop_back();                                                                              \
    EXPECT_EQ(logLine, getIndent() + " -40.75800000");                                             \
    logLine = lines.back();                                                                        \
    EXPECT_EQ(logLine, getIndent() + "   1.00450000");                                             \
  }                                                                                                \
                                                                                                   \
  TEST_F(loggingClass, StringTest) { test<std::string>("Test string", "Test string"); }            \
                                                                                                   \
  TEST_F(loggingClass, MultipleLineStringTest) {                                                   \
    getLogger() << std::string("Multiple line\ntest string") << std::endl;                         \
    std::vector<std::string> lines = ReadLines();                                                  \
    std::string logLine = lines.back();                                                            \
    lines.pop_back();                                                                              \
    EXPECT_EQ(logLine, getIndent() + "test string");                                               \
    logLine = lines.back();                                                                        \
    EXPECT_EQ(logLine, getIndent() + "Multiple line");                                             \
  }                                                                                                \
                                                                                                   \
  TEST_F(loggingClass, charArrayTest) {                                                            \
    test<const char *>("Test char array", "Test char array");                                      \
  }                                                                                                \
                                                                                                   \
  TEST_F(loggingClass, MultipleLineCharArrayTest) {                                                \
    getLogger() << "Multiple line\ntest char array" << std::endl;                                  \
    std::vector<std::string> lines = ReadLines();                                                  \
    std::string logLine = lines.back();                                                            \
    lines.pop_back();                                                                              \
    EXPECT_EQ(logLine, getIndent() + "test char array");                                           \
    logLine = lines.back();                                                                        \
    EXPECT_EQ(logLine, getIndent() + "Multiple line");                                             \
  }                                                                                                \
                                                                                                   \
  TEST_F(loggingClass, NewLineTest) {                                                              \
    getLogger() << "first" << std::endl << std::endl << "second" << std::endl;                     \
    std::vector<std::string> lines = ReadLines(IncludeEmptyLines::Yes);                            \
    std::string logLine = lines.back();                                                            \
    lines.pop_back();                                                                              \
    EXPECT_EQ(logLine, getIndent() + "second");                                                    \
    logLine = lines.back();                                                                        \
    lines.pop_back();                                                                              \
    EXPECT_TRUE(logLine.empty());                                                                  \
    logLine = lines.back();                                                                        \
    EXPECT_EQ(logLine, getIndent() + "first");                                                     \
                                                                                                   \
    getLogger() << "first" << std::endl                                                            \
                << "\n"                                                                            \
                << "second" << std::endl;                                                          \
    lines = ReadLines(IncludeEmptyLines::Yes);                                                     \
    logLine = lines.back();                                                                        \
    lines.pop_back();                                                                              \
    EXPECT_EQ(logLine, getIndent() + "second");                                                    \
    logLine = lines.back();                                                                        \
    lines.pop_back();                                                                              \
    EXPECT_EQ(logLine, getIndent());                                                               \
    logLine = lines.back();                                                                        \
    EXPECT_EQ(logLine, getIndent() + "first");                                                     \
                                                                                                   \
    getLogger() << "first" << std::endl                                                            \
                << "\n\n"                                                                          \
                << "second" << std::endl;                                                          \
    lines = ReadLines(IncludeEmptyLines::Yes);                                                     \
    logLine = lines.back();                                                                        \
    lines.pop_back();                                                                              \
    EXPECT_EQ(logLine, getIndent() + "second");                                                    \
    logLine = lines.back();                                                                        \
    lines.pop_back();                                                                              \
    EXPECT_EQ(logLine, getIndent());                                                               \
    logLine = lines.back();                                                                        \
    lines.pop_back();                                                                              \
    EXPECT_EQ(logLine, getIndent());                                                               \
    logLine = lines.back();                                                                        \
    EXPECT_EQ(logLine, getIndent() + "first");                                                     \
  }                                                                                                \
                                                                                                   \
  TEST_F(loggingClass, OperatorTest) {                                                             \
    getLogger() << DummyClass("dummy1");                                                           \
    getLogger() << DummyClass("dummy2");                                                           \
    std::vector<std::string> lines = ReadLines();                                                  \
    std::string logLine = lines[lines.size() - 2];                                                 \
    EXPECT_EQ(logLine, getIndent() + "DummyClass with name <dummy1>");                             \
    logLine = lines[lines.size() - 1];                                                             \
    EXPECT_EQ(logLine, getIndent() + "DummyClass with name <dummy2>");                             \
  }

//========================================================================================
// Mout tests
//========================================================================================
using MoutTest = TestLogging<MasterLogging>;

TEST_LOGGING(MoutTest);

TEST_F(MoutTest, BeginlTest) {
  mout << "first" << std::beginl << "second" << std::endl;
  std::vector<std::string> lines = ReadLines();
  std::string logLine = lines.back();
  lines.pop_back();
  EXPECT_EQ(logLine, getIndent() + "second");
  logLine = lines.back();
  EXPECT_EQ(logLine, getIndent() + "first");
}

TEST_F(MoutTest, BlockTest) {
  mout.StartBlock("TestBlock");
  mout << "output" << std::endl;
  std::string logLine = ReadLines().back();
  EXPECT_EQ(logLine, getIndent() + "TestBlock: output");
  mout.EndBlock();
  logLine = ReadLines().back();
  EXPECT_TRUE(logLine.starts_with(getIndent() + "TestBlock:"));
  EXPECT_TRUE(logLine.ends_with("seconds"));
}

TEST_F(MoutTest, MultipleBlockTest) {
  mout.StartBlock("TestBlock");
  mout << "output" << std::endl;
  std::string logLine = ReadLines().back();
  EXPECT_EQ(logLine, getIndent() + "TestBlock: output");

  mout.StartBlock("IndentedBlock1");
  mout << "output" << std::endl;
  logLine = ReadLines().back();
  EXPECT_EQ(logLine, getIndent() + "  IndentedBlock1: output");
  mout.EndBlock();
  logLine = ReadLines().back();
  EXPECT_TRUE(logLine.starts_with(getIndent() + "  IndentedBlock1:"));
  EXPECT_TRUE(logLine.ends_with("seconds"));

  mout.StartBlock("IndentedBlock2");
  mout << "output" << std::endl;
  logLine = ReadLines().back();
  EXPECT_EQ(logLine, getIndent() + "  IndentedBlock2: output");
  mout.EndBlock();
  logLine = ReadLines().back();
  EXPECT_TRUE(logLine.starts_with(getIndent() + "  IndentedBlock2:"));
  EXPECT_TRUE(logLine.ends_with("seconds"));

  mout << "output" << std::endl;
  logLine = ReadLines().back();
  EXPECT_EQ(logLine, getIndent() + "TestBlock: output");
  mout.EndBlock();
  logLine = ReadLines().back();
  EXPECT_TRUE(logLine.starts_with(getIndent() + "TestBlock:"));
  EXPECT_TRUE(logLine.ends_with("seconds"));
}

// JSON Tests

using JsonTest = TestLogging<MasterLogging>;

TEST_F(JsonTest, printEntry) {
  std::string testValue = "TestValueRaw";
  std::string testKey = "TestKeyRaw";
  mout.printEntryRaw(testKey, testValue);
  mout.flush();
  nlohmann::json j = ReadJson();
  EXPECT_EQ(j[testKey], testValue);
}

TEST_F(JsonTest, printInfo) {
  std::string testValue = "testValueInfo";
  std::string testKey = "testKeyInfo";
  mout.PrintInfo(testKey, 3, PrintInfoEntry(testKey + "1", testValue + "1"),
                 PrintInfoEntry(testKey + "2", testValue + "2"),
                 PrintInfoEntry(testKey + "3", testValue + "3"),
                 PrintInfoEntry(testKey + "4", std::vector<double>{0, 1, 2, 3}),
                 PrintInfoEntry(testKey + "5", std::vector<std::string>{"0", "1", "2", "3"}));

  mout.flush();
  nlohmann::json j = ReadJson();
  EXPECT_EQ(j[testKey + " Info"][testKey + "1"], testValue + "1");
  EXPECT_EQ(j[testKey + " Info"][testKey + "2"], testValue + "2");
  EXPECT_EQ(j[testKey + " Info"][testKey + "3"], testValue + "3");

  const auto &doubleArray = j[testKey + " Info"][testKey + "4"];
  const auto &stringArray = j[testKey + " Info"][testKey + "5"];
  EXPECT_EQ(doubleArray.size(), 4);
  EXPECT_EQ(stringArray.size(), 4);
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(doubleArray[i], i);
    EXPECT_EQ(stringArray[i], std::to_string(i));
  }
}

TEST_F(JsonTest, groups) {
  std::string testValue = "testValueGroup";
  std::string testKey = "testKeyGroup";

  mout.StartBlock(testKey);
  mout.StartBlock();
  mout.StartBlock(testKey);

  mout.printEntryRaw(testKey, testValue);
  mout.EndBlock();
  mout.EndBlock();
  mout.EndBlock();
  mout.flush();
  nlohmann::json j = ReadJson();

  EXPECT_EQ(j[testKey][testKey][testKey], testValue);
}

// The JSON Logger is expected to not store blocks
// that do not contain any data
TEST_F(JsonTest, dropEmptyGroups) {
  // Baseline test: We should not drop groups *with* data
  mout.StartBlock("doNotDropMe");
  mout.printEntryRaw("key", "value");
  mout.EndBlock();

  // But if the group is empty then it should *not*
  // appear in the final json
  mout.StartBlock("dropMe");
  mout.EndBlock();

  mout.flush();
  nlohmann::json j = ReadJson();

  EXPECT_TRUE(j.contains("doNotDropMe"));
  EXPECT_FALSE(j.contains("dropMe"));
}

// TEST_F(JsonTest, PointInJson){
//   PointT<double, 3, 1>(0.0, 0.0, 0.0, 0.0)
// }

#ifdef BUILD_IA
TEST_F(JsonTest, IAIntervalTest) {
  constexpr static double inf = 10;
  constexpr static double sup = 15;
  IAInterval testValue(inf, sup);
  std::string testKey = "testKeyIAInterval";
  mout.printEntryRaw(testKey, testValue);
  mout.flush();
  nlohmann::json j = ReadJson();

  EXPECT_EQ(j[testKey]["inf"], inf);
  EXPECT_EQ(j[testKey]["sup"], sup);
}

TEST_F(JsonTest, IACIntervalTest) {
  constexpr static double rel = 10;
  constexpr static double im = 15;
  IACInterval testValue(rel, im);
  std::string testKey = "testKeyIACInterval";
  mout.printEntryRaw(testKey, testValue);
  mout.flush();
  nlohmann::json j = ReadJson();

  EXPECT_EQ(j[testKey]["rel"]["inf"], rel);
  EXPECT_EQ(j[testKey]["rel"]["sup"], rel);
  EXPECT_EQ(j[testKey]["im"]["inf"], im);
  EXPECT_EQ(j[testKey]["im"]["sup"], im);
}
#endif

// TEST_F(MoutTest, WarningTest) {
//   Warning("Test warning")
//   int lineNumber = __LINE__ - 1;
//   std::vector<std::string> lines = ReadLines();
//   std::string logLine = lines[lines.size() - 3];
//   EXPECT_EQ(logLine, "Warning: Test warning");
//   logLine = lines[lines.size() - 2];
//   EXPECT_TRUE(logLine.ends_with("tests/test0_basic/utility/TestLogging.cpp:"
//   + std::to_string(lineNumber))); logLine = lines[lines.size() - 1];
//   EXPECT_EQ(logLine, "         on line " + std::to_string(lineNumber));
//
//   WarningOnProc("Test parallel warning")
//   //TODO: to check at end of file
// }

//========================================================================================
// Pout tests
//========================================================================================
using PoutTest = TestLogging<LogTextStream<ProcessLogging>>;

TEST_LOGGING(PoutTest);

int main(int argc, char **argv) {
  Config::SetLogPath(Config::GetBuildPath() + "/tests/test0_basic/");
  Config::SetJsonPath(Config::GetBuildPath() + "/tests/test0_basic/");
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithFileLogging().WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}
