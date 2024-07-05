#include "MemoryLogger.hpp"
#include "TestEnvironment.hpp"

// Tolerance of 5 %
const double MEMORY_TOLERANCE = 0.05;

TEST(TestMemoryLogger, TestSetup) {
  MemoryLogger::Tare();
  MemoryLogger::LogMemory("Setup Memory");
  // One int is 4 Byte, hence the line below adds 4 MByte * 16 per process
  int amountOfInt = 1024 * 1024 * 16; // This gives 64 MByte per process
  std::vector<int> myVector(amountOfInt, 0);
  MemoryLogger::LogMemory("added std::vector<int>");
  EXPECT_LE(MemoryLogger::GetMemoryStampVector().back(),
            PPM->Size(0) * 64 * (1 + MEMORY_TOLERANCE));

  std::vector<int> mySecVector(amountOfInt, 0);
  MemoryLogger::LogMemory("added std::vector<int>");
  // Now it should be 128 MByte per process
  EXPECT_LE(MemoryLogger::GetMemoryStampVector().back(),
            PPM->Size(0) * 128 * (1 + MEMORY_TOLERANCE));

  MemoryLogger::Tare();
  // Now its 0 MByte per process again
  std::vector<int> myThirdVector(amountOfInt, 0);
  MemoryLogger::LogMemory("added std::vector<int>");
  // And finally 64 MByte per process
  EXPECT_LE(MemoryLogger::GetMemoryStampVector().back(),
            PPM->Size(0) * 64 * (1 + MEMORY_TOLERANCE));

  MemoryLogger::PrintInfo();
}

int main(int argc, char **argv) {
  return MppTest(MppTestBuilder(argc, argv)
                     .WithConfigEntry("MemoryLoggerVerbose", 2)
                     .WithScreenLogging()
                     .WithPPM())
      .RUN_ALL_MPP_TESTS();
}