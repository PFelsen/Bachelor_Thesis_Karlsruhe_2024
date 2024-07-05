#include "Logging.hpp"

#include <stack>

#include "Config.hpp"
#include "ExchangeBuffer.hpp"
#include "Parallel.hpp"
#include "TimeDate.hpp"

std::unique_ptr<MasterLogging> MasterLogging::instance;
LogTextStream<ProcessLogging> *ProcessLogging::instance = nullptr;

MasterLogging::MasterLogging(EnabledLoggers enabledLoggers) {
  onMaster = PPM->Master(0);
  this->setEnabledLoggers(enabledLoggers);

  // Build info is only show in json, so we temporarily
  // disable the text logger (regardless of whether or not it
  // was previously enabled)
  EnabledLoggers oldConfiguration = enabledLoggers;
  enabledLoggers[LoggingBackend::Text] = false;
  this->setEnabledLoggers(enabledLoggers);
  PrintBuildInfo();
  this->setEnabledLoggers(oldConfiguration);
}

JsonLogger::JsonLogger() {
  // FIXME: If we are not the master process, we shouldn't initialize
  //        the logging backends in the first place.
  onMaster = PPM->Master(0);
  if (onMaster) {
    printIterationContext = {};
    json internalStorage = json(json::value_t::object);
    openBlocks.push_back(JsonBlock{internalStorage, 0});
    fileOut.open((Config::GetJsonFileNameWithPath() + ".json").c_str(), "json");
  }
}

void JsonLogger::StartBlock(const std::string &blockName = "") {
  const auto start = std::chrono::system_clock::now();
  if (blockName == "") {
    openBlocks.back().unnamedOpenGroups++;
  } else {
    json block(json::value_t::object);
    openBlocks.push_back(JsonBlock{block, 0, 0, blockName, start});
  }
}

void JsonLogger::EndBlock() {
  if (openBlocks.back().unnamedOpenGroups > 0) {
    openBlocks.back().unnamedOpenGroups--;
  } else {
    // NOTE: The top-level block is implicitly created. The user-created
    //       groups are placed inside it, but it should *never* be closed
    Assert(openBlocks.size() != 1);

    auto closedBlock = openBlocks.back();
    openBlocks.pop_back();

    if (closedBlock.terminatesPrintIteration) endPrintIterationContext();

    if (openBlocks.size() == 1) { dirty = true; }

    // Empty groups are dropped
    if (!closedBlock.data.empty()) {

      // check if it is a PoutBlock
      if (!closedBlock.data.is_array()) {
        closedBlock.data["time"] =
            (std::chrono::system_clock::now() - closedBlock.startTime).count();
      }
      openBlocks.back().data[closedBlock.blockName] = closedBlock.data;
    }
  }
}

void JsonLogger::flush() {
  if (onMaster && dirty) {
    fileOut.seekp(0, std::ios_base::beg);
    fileOut << std::setw(2) << openBlocks.front().data << std::endl;
    dirty = false;
  }
}

JsonLogger::~JsonLogger() {
  if (onMaster) endPrintIterationContext();
  flush();
}

TextLogger::TextLogger() {
  // FIXME: If we are not the master process, we shouldn't initialize
  //        the logging backends in the first place.
  onMaster = PPM->Master();
  if (onMaster) {
    // Load settings from config
    Config::Get("enableFileLogging", fileEnabled, true);
    Config::Get("enableScreenLogging", screenEnabled, true);

    int precision = 8;
    Config::Get("precision", precision, true);
    setPrecision(precision);

    if (fileEnabled) fileOut.open((Config::GetLogFileNameWithPath() + ".log").c_str(), "log");

    // Log start information
    // + 1 because HOST_NAME_MAX does not include the terminating null byte
    std::string hostname(HOST_NAME_MAX + 1, '\0');
    gethostname(hostname.data(), hostname.size());
    hostname.resize(hostname.find_first_of('\0'));

    logPrefix = "   ";

    startTime = Date();
    numProcs = PPM->Size();
    // std::format would be quite nice to have...
    this->write("start program on ");
    this->write(numProcs);
    this->write(" procs at ");
    this->write(startTime);
    *this << std::endl;
    this->write("Running on: ");
    this->write(hostname);
    *this << std::endl;
  }
}

std::string getMsg(int error_cnt, int warning_cnt) {
  if (error_cnt != 0) {
    if (warning_cnt != 0) {
      return " with " + std::to_string(error_cnt) + " error(s) and " + std::to_string(warning_cnt)
             + " warning(s):\n";
    } else {
      return " with " + std::to_string(error_cnt) + " error(s):\n";
    }
  } else {
    if (warning_cnt != 0) {
      return " with " + std::to_string(warning_cnt) + " warning(s):\n";
    } else {
      return "";
    }
  }
}

TextLogger::~TextLogger() {
  if (!onMaster) return;

  std::string msg = getMsg(error_cnt, warning_cnt);

  if (error_cnt != 0) msg += error_msg;
  if (warning_cnt != 0) msg += warning_msg;

  Date endTime = Date();
  this->write("end program after ");
  this->write(endTime - startTime);
  this->write(" on ");
  this->write(numProcs);
  this->write(" procs at ");
  this->write(endTime);
  this->write(msg);
  this->endl();
}

void TextLogger::StartBlock(const std::string &blockName) {
  const auto start = std::chrono::system_clock::now();
  if (blockName.empty()) openBlocks.push(TextBlock{start, {}});
  else openBlocks.push(TextBlock{start, blockName});

  if (openBlocks.size() > 1) IncreaseIndent();
}

void TextLogger::StartPrintInfoBlock(const std::string &blockName) {
  *this << blockName << ":" << std::endl;
  IncreaseIndent();
}

void TextLogger::EndBlock(bool mute, const std::string &msg) {
  Assert(!openBlocks.empty());

  if (!mute) {
    *this << msg;
    *this << std::chrono::system_clock::now() - openBlocks.top().startTime;
    *this << std::endl;
  }

  if (openBlocks.size() > 1) DecreaseIndent();

  openBlocks.pop();
}

void TextLogger::EndBlock(int classVerbose, int verbose) { EndBlock(classVerbose < verbose, ""); }

void TextLogger::StartBlockWithoutIndent(const std::string &blockName) {
  const auto start = std::chrono::system_clock::now();
  if (blockName.empty()) openBlocks.push(TextBlock{start, {}});
  else openBlocks.push(TextBlock{start, blockName});
}

void TextLogger::EndBlockWithoutIndent(bool mute) {
  Assert(!openBlocks.empty());

  if (!mute) {
    *this << "time ";
    *this << openBlocks.top().startTime << std::endl;
  }

  openBlocks.pop();
}

void TextLogger::AddErrorMsg(const std::string &msg) {
  if (msg.size() == 0) return;
  error_msg += msg;
  error_cnt++;
}

void TextLogger::AddWarningMsg(const std::string &msg) {
  if (msg.size() == 0) return;
  warning_msg += msg;
  warning_cnt++;
}

TextLogger &TextLogger::operator<<(std::ostream &(*f)(std::ostream &)) {
  if (f == (std::basic_ostream<char> & (*)(std::basic_ostream<char> &)) & std::endl) {
    if (fileEnabled) fileOut << std::endl;
    if (screenEnabled) std::cout << std::endl;
    lineEmpty = true;
  } else if (f == (std::basic_ostream<char> & (*)(std::basic_ostream<char> &)) & std::flush) {
    if (fileEnabled) fileOut << std::flush;
    if (screenEnabled) std::cout << std::flush;
  } else if (f == (std::basic_ostream<char> & (*)(std::basic_ostream<char> &)) & std::beginl) {
    if (fileEnabled) fileOut << "\n" << std::flush;
    if (screenEnabled) std::cout << "\n" << std::flush;
    lineEmpty = true;
  } else if (f == (std::basic_ostream<char> & (*)(std::basic_ostream<char> &)) & std::beginD) {
    doubleBlock = true;
  } else if (f == (std::basic_ostream<char> & (*)(std::basic_ostream<char> &)) & std::endD) {
    doubleBlock = false;
  } else {
    THROW("Not implemented in TextLogger!")
  }

  return *this;
}

TextLogger &TextLogger::operator<<(const double &s) {
  if (doubleBlock) {
    double absValue = std::abs(s);
    if (std::abs(s - infty) < 1.0) {
      int l = getPrecision() / 2;
      std::string msg(l, ' ');
      msg.append("infty").append(std::string(getPrecision() - l, ' '));
      (*this) << msg;
    } else if (std::abs(s + infty) < 1.0) {
      int l = getPrecision() / 2;
      std::string msg(l - 1, ' ');
      msg.append("-infty ").append(std::string(std::max(getPrecision() - l - 1, 0), ' '));
      (*this) << msg;
    } else if (s == 0.0 || (0.001 < absValue && absValue < 1000.0)) {
      char buf[32];
      sprintf(buf, doubleOutputFormat.c_str(), s);
      (*this) << buf;
    } else {
      if (s >= 0) (*this) << " ";

      int precision = getPrecision();
      if (screenEnabled) {
        std::cout.precision(precision - 1);
        std::cout << s;
        std::cout.precision(precision);
      }
      if (fileEnabled) {
        fileOut.precision(precision - 1);
        fileOut << s;
        fileOut.precision(precision);
      }
    }
  } else {
    printIndentMsg();
    this->write(s);
  }
  return *this;
}

void TextLogger::printFromPLogging(const std::string &line_p, int proc) {
  if (!lineEmpty) { (*this) << std::endl; }

  if (line_p != "") {
    std::string prefix = std::to_string(proc) + ": " + indent;
    if (!openBlocks.empty()) prefix += openBlocks.top().name.value_or("") + ": ";

    // Log each line individually
    std::vector<std::string> parts_p = splitWithoutTrim(line_p, "\n");
    for (const auto &messageLine : parts_p) {
      this->write(prefix);
      this->write(messageLine);
      *this << std::endl;
    }
  } else *this << std::endl;
}

ProcessLogging &ProcessLogging::operator<<(std::ostream &(*f)(std::ostream &)) {
  // Synchronize our buffer with the main process whenever a newline
  // is written to pout
  if (f == (std::basic_ostream<char> & (*)(std::basic_ostream<char> &)) & std::endl) {
    // Send our buffer to master
    ExchangeBuffer exBuffer(0);
    exBuffer.Send(0) << buffer.str();
    exBuffer.Communicate();

    if (PPM->Master()) {
      mout.startPoutBlock();
      // Collect and log the buffers from all child processes
      for (short p = 0; p < PPM->Size(0); ++p) {
        std::string line_p = "";
        exBuffer.Receive(p) >> line_p;

        if (line_p != "" || p == 0) mout.printFromPLogging(line_p, p);
      }
      mout.endPoutBlock();
    }

    // Clear the buffer
    buffer.str("");
    buffer.clear();
    return *this;
  }
  if (f == (std::basic_ostream<char> & (*)(std::basic_ostream<char> &)) & std::beginD) {
    doubleBlock = true;
    return *this;
  }
  if (f == (std::basic_ostream<char> & (*)(std::basic_ostream<char> &)) & std::endD) {
    doubleBlock = false;
    return *this;
  }
  THROW("Not implemented")
}

ProcessLogging &ProcessLogging::operator<<(const double &s) {
  if (doubleBlock) {
    double absValue = std::abs(s);
    if (std::abs(s - infty) < 1.0) {
      int l = getPrecision() / 2;
      std::string msg(l + 1, ' ');
      msg.append("infty").append(std::string(getPrecision() - l + 1, ' '));
      buffer << msg;
    } else if (std::abs(s + infty) < 1.0) {
      int l = getPrecision() / 2;
      std::string msg(l, ' ');
      msg.append("-infty ").append(std::string(std::max(getPrecision() - l, 0), ' '));
      buffer << msg;
    } else if (s == 0.0 || (0.001 < absValue && absValue < 1000.0)) {
      char buf[32];
      sprintf(buf, doubleOutputFormat.c_str(), s);
      buffer << buf;
    } else {
      if (s >= 0) buffer << " ";
      buffer.precision(getPrecision() - 1);
      buffer << s;
      buffer.precision(getPrecision());
    }
  } else {
    buffer << s;
  }
  return *this;
}

void MasterLogging::CommunicateProcessEnvironment() {
  // Get the current hostname (unix-specific!)
  // + 1 because HOST_NAME_MAX does not include the terminating null byte
  std::string hostname(HOST_NAME_MAX + 1, '\0');
  gethostname(hostname.data(), hostname.size());
  hostname.resize(hostname.find_first_of('\0'));

  ExchangeBuffer exBuffer(0);

  // Tell the master about our hostname
  exBuffer.Send(0) << hostname;
  exBuffer.Communicate();

  if (PPM->Master()) {
    // Create a mapping from <host name> to <number of processes on that host>
    std::unordered_map<std::string, int> environments;
    environments.reserve(8);

    for (short p = 0; p < PPM->Size(0); ++p) {
      std::string processHost;
      exBuffer.Receive(p) >> processHost;
      environments[processHost] += 1;
    }

    // Log the process counts for each host
    std::vector<PrintInfoEntry<int>> entries;
    entries.reserve(environments.size());
    for (const auto &p : environments)
      entries.push_back(PrintInfoEntry<int>("Processes on \"" + p.first + "\"", p.second));
    std::sort(entries.begin(), entries.end());
    this->PrintInfo("Process environment", 1, entries);
  }
}

template<>
std::string vec2str<int>(std::vector<int> vec) {
  std::string str = "[";
  for (unsigned long i = 0; i < vec.size(); i++) {
    if (i != vec.size() - 1) str += std::to_string(vec[i]) + ", ";
    else str += std::to_string(vec[i]);
  }
  return (str + "]");
}

template<>
std::string vec2str<double>(std::vector<double> vec) {
  std::string str = "[";
  for (unsigned long i = 0; i < vec.size(); i++) {
    std::ostringstream streamObj;
    if (i != vec.size() - 1) {
      streamObj << vec[i];
      str += streamObj.str() + ", ";
    } else {
      streamObj << vec[i];
      str += streamObj.str();
    }
  }
  return (str + "]");
}

template<>
std::string vec2str<std::string>(std::vector<std::string> vec) {
  std::string str = "[";
  for (unsigned long i = 0; i < vec.size(); i++) {
    if (i != vec.size() - 1) str += vec[i] + ", ";
    else str += vec[i];
  }
  return (str + "]");
}

template<>
std::string vec2str<int>(std::vector<std::vector<int>> vec) {
  std::string str = "[";
  for (unsigned long i = 0; i < vec.size(); i++) {
    std::ostringstream streamObj;
    if (i != vec.size() - 1) {
      streamObj << vec2str(vec[i]);
      str += streamObj.str() + ", ";
    } else {
      streamObj << vec2str(vec[i]);
      str += streamObj.str();
    }
  }
  return (str + "]");
}

template<>
std::string vec2str<double>(std::vector<std::vector<double>> vec) {
  std::string str = "[";
  for (unsigned long i = 0; i < vec.size(); i++) {
    std::ostringstream streamObj;
    if (i != vec.size() - 1) {
      streamObj << vec2str(vec[i]);
      str += streamObj.str() + ", ";
    } else {
      streamObj << vec2str(vec[i]);
      str += streamObj.str();
    }
  }
  return (str + "]");
}

template<>
std::string vec2str<std::string>(std::vector<std::vector<std::string>> vec) {
  std::string str = "[";
  for (unsigned long i = 0; i < vec.size(); i++) {
    std::ostringstream streamObj;
    if (i != vec.size() - 1) {
      streamObj << vec2str(vec[i]);
      str += streamObj.str() + ", ";
    } else {
      streamObj << vec2str(vec[i]);
      str += streamObj.str();
    }
  }
  return (str + "]");
}