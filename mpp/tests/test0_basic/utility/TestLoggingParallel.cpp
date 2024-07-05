#include "Logging.hpp"
#include "TestEnvironment.hpp"

#include <fstream>

std::vector<std::string> ReadLines() {
  std::ifstream infile(Config::GetBuildPath() + "/tests/test0_basic/logfile.log");
  std::string line;
  std::vector<std::string> lines;
  while (std::getline(infile, line)) {
    if (!line.empty()) lines.push_back(line);
  }
  return lines;
}

std::vector<std::string> ReadLinesWithEmpty() {
  std::ifstream infile(Config::GetBuildPath() + "/tests/test0_basic/logfile.log");
  std::string line;
  std::vector<std::string> lines;
  while (std::getline(infile, line)) {
    lines.push_back(line);
  }
  return lines;
}

//========================================================================================
// Mout tests
//========================================================================================

class ParallelMoutTest : public ::testing::Test {
protected:
  std::string indent = "";

  ParallelMoutTest() {
    if (PPM->Size(0) > 1) {
      int emptyLength = (int)(std::log10(PPM->Size(0) - 1)) + 1;
      indent = std::string(emptyLength, ' ') + "  ";
    }
  }

  template<typename T>
  void test(T value, std::string expectedLine) {
    mout << value << std::endl;
    std::string logLine;
    if (PPM->Master()) {
      logLine = ReadLines().back();
      EXPECT_EQ(logLine, indent + expectedLine);
    }
    //    mout << value << "\n" << std::flush;
    //    if (PPM->Master()) {
    //      logLine = ReadLines().back();
    //      EXPECT_EQ(logLine, indent + expectedLine);
    //    }
  }

  void testDoubleBlock(double value, std::string expectedLine) {
    mout << std::beginD;
    test<double>(value, expectedLine);
    mout << std::endD;
  }
};

TEST_F(ParallelMoutTest, ShortTest) {
  test<short>(257, "257");
  test<short>(-4830, "-4830");
}

TEST_F(ParallelMoutTest, IntegerTest) {
  test<int>(10542, "10542");
  test<int>(-34861, "-34861");
}

TEST_F(ParallelMoutTest, DoubleTest) {
  test<double>(1.8754, "1.8754");
  test<double>(-43.7462, "-43.7462");
}

TEST_F(ParallelMoutTest, DoubleBlockTest) {
  testDoubleBlock(0.0, "   0.00000000");
  testDoubleBlock(-0.0, "  -0.00000000");
  testDoubleBlock(infty, "    infty    ");
  testDoubleBlock(-infty, "   -infty    ");
  testDoubleBlock(1e-2, "   0.01000000");
  testDoubleBlock(-1e-2, "  -0.01000000");
  testDoubleBlock(1e2, " 100.00000000");
  testDoubleBlock(-1e2, "-100.00000000");
}

TEST_F(ParallelMoutTest, MultipleLineDoubleBlockTest) {
  mout << std::beginD << 1.0045 << std::endl << -40.758 << "\n" << 5.098 << std::endl << std::endD;
  if (PPM->Master()) {
    std::vector<std::string> lines = ReadLines();
    std::string logLine = lines.back();
    lines.pop_back();
    EXPECT_EQ(logLine, indent + "   5.09800000");
    logLine = lines.back();
    lines.pop_back();
    EXPECT_EQ(logLine, indent + " -40.75800000");
    logLine = lines.back();
    EXPECT_EQ(logLine, indent + "   1.00450000");
  }
}

TEST_F(ParallelMoutTest, StringTest) { test<std::string>("Test string", "Test string"); }

TEST_F(ParallelMoutTest, MultipleLineStringTest) {
  mout << std::string("Multiple line\ntest string") << std::endl;
  if (PPM->Master()) {
    std::vector<std::string> lines = ReadLines();
    std::string logLine = lines.back();
    lines.pop_back();
    EXPECT_EQ(logLine, indent + "test string");
    logLine = lines.back();
    EXPECT_EQ(logLine, indent + "Multiple line");
  }
}

TEST_F(ParallelMoutTest, charArrayTest) {
  test<const char *>("Test char array", "Test char array");
}

TEST_F(ParallelMoutTest, MultipleLineCharArrayTest) {
  mout << "Multiple line\ntest char array" << std::endl;
  if (PPM->Master()) {
    std::vector<std::string> lines = ReadLines();
    std::string logLine = lines.back();
    lines.pop_back();
    EXPECT_EQ(logLine, indent + "test char array");
    logLine = lines.back();
    EXPECT_EQ(logLine, indent + "Multiple line");
  }
}

TEST_F(ParallelMoutTest, BlockTest) {
  mout.StartBlock("TestBlock");
  mout << "output" << std::endl;
  std::string logLine;
  if (PPM->Master()) {
    logLine = ReadLines().back();
    EXPECT_EQ(logLine, indent + "TestBlock: output");
  }
  mout.EndBlock();
  if (PPM->Master()) {
    logLine = ReadLines().back();
    EXPECT_TRUE(logLine.starts_with(indent + "TestBlock:"));
    EXPECT_TRUE(logLine.ends_with("seconds"));
  }
}

TEST_F(ParallelMoutTest, MultipleBlockTest) {
  mout.StartBlock("TestBlock");
  mout << "output" << std::endl;
  std::string logLine;
  if (PPM->Master()) {
    logLine = ReadLines().back();
    EXPECT_EQ(logLine, indent + "TestBlock: output");
  }

  mout.StartBlock("IndentedBlock1");
  mout << "output" << std::endl;
  if (PPM->Master()) {
    logLine = ReadLines().back();
    EXPECT_EQ(logLine, indent + "  IndentedBlock1: output");
  }
  mout.EndBlock();
  if (PPM->Master()) {
    logLine = ReadLines().back();
    EXPECT_TRUE(logLine.starts_with(indent + "  IndentedBlock1:"));
    EXPECT_TRUE(logLine.ends_with("seconds"));
  }

  mout.StartBlock("IndentedBlock2");
  mout << "output" << std::endl;
  if (PPM->Master()) {
    logLine = ReadLines().back();
    EXPECT_EQ(logLine, indent + "  IndentedBlock2: output");
  }
  mout.EndBlock();
  if (PPM->Master()) {
    logLine = ReadLines().back();
    EXPECT_TRUE(logLine.starts_with(indent + "  IndentedBlock2:"));
    EXPECT_TRUE(logLine.ends_with("seconds"));
  }

  mout << "output" << std::endl;
  if (PPM->Master()) {
    logLine = ReadLines().back();
    EXPECT_EQ(logLine, indent + "TestBlock: output");
  }
  mout.EndBlock();
  if (PPM->Master()) {
    logLine = ReadLines().back();
    EXPECT_TRUE(logLine.starts_with(indent + "TestBlock:"));
    EXPECT_TRUE(logLine.ends_with("seconds"));
  }
}

TEST_F(ParallelMoutTest, BeginlTest) {
  mout << "first" << std::beginl << "second" << std::endl;
  if (PPM->Master()) {
    std::vector<std::string> lines = ReadLines();
    std::string logLine = lines.back();
    lines.pop_back();
    EXPECT_EQ(logLine, indent + "second");
    logLine = lines.back();
    EXPECT_EQ(logLine, indent + "first");
  }
}

// TEST_F(ParallelMoutTest, WarningTest) {
//   Warning("Test warning")
//   int lineNumber = __LINE__ - 1;
//   if (PPM->Master()) {
//     std::vector<std::string> lines = ReadLines();
//     std::string logLine = lines[lines.size() - 3];
//     EXPECT_EQ(logLine, "Warning: Test warning");
//     logLine = lines[lines.size() - 2];
//     EXPECT_TRUE(logLine.ends_with(
//         "tests/test0_basic/utility/TestParallelLogging.cpp:" + std::to_string(lineNumber)));
//     logLine = lines[lines.size() - 1];
//     EXPECT_EQ(logLine, "         on line " + std::to_string(lineNumber));
//   }
//
//   WarningOnProc("Parallel warning")
//   //TODO: to check at end of file
// }

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

TEST_F(ParallelMoutTest, OperatorTest) {
  mout << DummyClass("dummy1");
  mout << DummyClass("dummy2");
  if (PPM->Master()) {
    std::vector<std::string> lines = ReadLines();
    std::string logLine = lines[lines.size() - 2];
    EXPECT_EQ(logLine, indent + "DummyClass with name <dummy1>");
    logLine = lines[lines.size() - 1];
    EXPECT_EQ(logLine, indent + "DummyClass with name <dummy2>");
  }
}

//========================================================================================
// Pout tests
//========================================================================================

class ParallelPoutTest : public ::testing::Test {
protected:
  int emptyLength = 0;

  ParallelPoutTest() {
    if (PPM->Size(0) > 1) { emptyLength = (int)(std::log10(PPM->Size(0) - 1)) + 1; }
  }

  std::string getIndent(int proc) {
    int pL = 1;
    if (proc > 0) pL = (int)(std::log10(proc)) + 1;
    return std::to_string(proc) + ":" + std::string(emptyLength - pL, ' ') + " ";
  }

  void startsWithOnProcs(const std::string &expectedLine) {
    if (!PPM->Master()) return;
    std::vector<std::string> lines = ReadLines();
    for (int p = 0; p < PPM->Size(0); ++p) {
      std::string logLine = lines[lines.size() - PPM->Size() + p];
      EXPECT_TRUE(logLine.starts_with(getIndent(p) + expectedLine));
    }
  }

  void endsWithOnProcs(const std::string &expectedLine) {
    if (!PPM->Master()) return;
    std::vector<std::string> lines = ReadLines();
    for (int p = 0; p < PPM->Size(0); ++p) {
      std::string logLine = lines[lines.size() - PPM->Size() + p];
      EXPECT_TRUE(logLine.ends_with(expectedLine));
    }
  }

  void testOnProcs(const std::string &expectedLine) {
    if (!PPM->Master()) return;
    std::vector<std::string> lines = ReadLines();
    for (int p = 0; p < PPM->Size(0); ++p) {
      std::string logLine = lines[lines.size() - PPM->Size() + p];
      EXPECT_EQ(logLine, getIndent(p) + expectedLine);
    }
  }

  void testOnProcs2Lines(const std::string &expectedLine1, const std::string &expectedLine2,
                         bool singleBlock) {
    if (!PPM->Master()) return;
    std::vector<std::string> lines = ReadLines();
    if (singleBlock) {
      for (int p = 0; p < PPM->Size(0); ++p) {
        std::string logLine = lines[lines.size() - 2 * PPM->Size() + 2 * p];
        EXPECT_EQ(logLine, getIndent(p) + expectedLine1);
        logLine = lines[lines.size() - 2 * PPM->Size() + 2 * p + 1];
        EXPECT_EQ(logLine, getIndent(p) + expectedLine2);
      }
    } else {
      for (int p = 0; p < PPM->Size(0); ++p) {
        std::string logLine = lines[lines.size() - 2 * PPM->Size() + p];
        EXPECT_EQ(logLine, getIndent(p) + expectedLine1);
        logLine = lines[lines.size() - PPM->Size() + p];
        EXPECT_EQ(logLine, getIndent(p) + expectedLine2);
      }
    }
  }

  void testOnProcs2Lines(const std::string &expectedLine, bool singleBlock) {
    testOnProcs2Lines(expectedLine, "second line", singleBlock);
  }

  template<typename T>
  void test(T value, std::string expectedLine) {
    pout << value << std::endl;
    testOnProcs(expectedLine);
    pout << value << "\n";
    pout << "second line" << std::endl;
    testOnProcs2Lines(expectedLine, true);
    pout << value << std::endl;
    pout << "second line" << std::endl;
    testOnProcs2Lines(expectedLine, false);
  }

  void testDoubleBlock(double value, std::string expectedLine) {
    pout << std::beginD;
    test<double>(value, expectedLine);
    pout << std::endD;
  }
};

TEST_F(ParallelPoutTest, ShortTest) {
  test<short>(257, "257");
  test<short>(-4830, "-4830");
}

TEST_F(ParallelPoutTest, IntegerTest) {
  test<int>(10542, "10542");
  test<int>(-34861, "-34861");
}

TEST_F(ParallelPoutTest, DoubleTest) {
  test<double>(1.8754, "1.8754");
  test<double>(-43.7462, "-43.7462");
}

TEST_F(ParallelPoutTest, DoubleBlockTest) {
  testDoubleBlock(0.0, "   0.00000000");
  testDoubleBlock(-0.0, "  -0.00000000");
  testDoubleBlock(infty, "    infty    ");
  testDoubleBlock(-infty, "   -infty    ");
  testDoubleBlock(1e-2, "   0.01000000");
  testDoubleBlock(-1e-2, "  -0.01000000");
  testDoubleBlock(1e2, " 100.00000000");
  testDoubleBlock(-1e2, "-100.00000000");
}

TEST_F(ParallelPoutTest, MultipleLineDoubleBlockTest) {
  pout << std::beginD << 1.0045 << std::endl << -40.758 << "\n" << 5.098 << std::endl << std::endD;
  if (PPM->Master()) {
    std::vector<std::string> lines = ReadLines();
    for (int p = 0; p < PPM->Size(0); ++p) {
      std::string logLine = lines[lines.size() - 3 * PPM->Size() + p];
      EXPECT_EQ(logLine, getIndent(p) + "   1.00450000");
      logLine = lines[lines.size() - 2 * PPM->Size() + 2 * p];
      EXPECT_EQ(logLine, getIndent(p) + " -40.75800000");
      logLine = lines[lines.size() - 2 * PPM->Size() + 2 * p + 1];
      EXPECT_EQ(logLine, getIndent(p) + "   5.09800000");
    }
  }
}

TEST_F(ParallelPoutTest, StringTest) { test<std::string>("Test string", "Test string"); }

TEST_F(ParallelPoutTest, MultipleLineStringTest) {
  pout << std::string("Multiple line\ntest string") << std::endl;
  testOnProcs2Lines("Multiple line", "test string", true);
}

TEST_F(ParallelPoutTest, charArrayTest) {
  test<const char *>("Test char array", "Test char array");
}

TEST_F(ParallelPoutTest, MultipleLineCharArrayTest) {
  pout << "Multiple line\ntest char array" << std::endl;
  testOnProcs2Lines("Multiple line", "test char array", true);
}

// TODO: Blocks not implemented for ProcessLogging
// TEST_F(ParallelPoutTest, BlockTest) {
//   pout.StartBlock("TestBlock");
//   pout << "output" << std::endl;
//   testOnProcs("TestBlock: output");
//   pout.EndBlock();
//   startsWithOnProcs("TestBlock:");
//   endsWithOnProcs("seconds");
// }

// TEST_F(ParallelPoutTest, MultipleBlockTest) {
//   pout.StartBlock("TestBlock");
//   pout << "output" << std::endl;
//   testOnProcs("TestBlock: output");
//
//   pout.StartBlock("IndentedBlock1");
//   pout << "output" << std::endl;
//   testOnProcs("  IndentedBlock1: output");
//   pout.EndBlock();
//   startsWithOnProcs("  IndentedBlock1:");
//   endsWithOnProcs("seconds");
//
//   pout.StartBlock("IndentedBlock2");
//   pout << "output" << std::endl;
//   testOnProcs("  IndentedBlock2: output");
//   pout.EndBlock();
//   startsWithOnProcs("  IndentedBlock2:");
//   endsWithOnProcs("seconds");
//
//   pout << "output" << std::endl;
//   testOnProcs("TestBlock: output");
//   pout.EndBlock();
//   startsWithOnProcs("TestBlock:");
//   endsWithOnProcs("seconds");
// }
// TODO: beginl not implemented for ProcessLogging
// TEST_F(ParallelPoutTest, BeginlTest) {
//   pout << "first" << std::beginl << "second" << std::endl;
//   testOnProcs2Lines("first", "second", true);
// }

TEST_F(ParallelPoutTest, NewLineTest) {
  pout << "first" << std::endl << std::endl << "second" << std::endl;
  testOnProcs("second");
  if (PPM->Master(0)) {
    std::vector<std::string> lines = ReadLinesWithEmpty();
    EXPECT_TRUE(lines[lines.size() - PPM->Size(0) - 1].empty());
  }
}

int main(int argc, char **argv) {
  Config::SetLogPath(Config::GetBuildPath() + "/tests/test0_basic/");
  MppTest mppTest = MppTestBuilder(argc, argv).WithPPM().WithFileLogging().WithScreenLogging();
  return mppTest.RUN_ALL_MPP_TESTS();
}
