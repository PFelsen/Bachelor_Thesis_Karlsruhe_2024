#ifndef MPP_MEMORYLOGGER_H
#define MPP_MEMORYLOGGER_H

#include <utility>
#include <vector>

#include "Config.hpp"
#include "Parallel.hpp"

struct MemoryStamp {
  MemoryStamp(std::string description, double memoryInMByte, double timeSinceInitInSeconds) :
      description(std::move(description)), memorySinceTareInMByte(memoryInMByte),
      timeSinceTareInSeconds(timeSinceInitInSeconds) {}

  double memorySinceTareInMByte;

  double timeSinceTareInSeconds;

  std::string description;
};

class MemoryLogger {
public:
  MemoryLogger() = delete;

  MemoryLogger(const MemoryLogger &random) = delete;

  MemoryLogger &operator=(const MemoryLogger &random) = delete;

  // returns the physical memory im MB for all processes
  // sends to master before logging
  static double MemoryUsage(bool perProcess = false);

  // returns the physical memory im MB for all processes
  // uses pout for logging
  static double TotalMemoryUsage(bool perProcess = false);

  static void LogMemory(const std::string &description = "");

  static void Tare();

  static void PrintInfo();

  static std::vector<double> GetTimeStampVector();

  static std::vector<double> GetMemoryStampVector();

  static std::vector<std::string> GetDescriptionVector();
private:
  static int verbose;

  static double initTime;

  static double offSetMemory;

  static std::vector<MemoryStamp> memoryStamps;
};

#endif // MPP_MEMORYLOGGER_H
