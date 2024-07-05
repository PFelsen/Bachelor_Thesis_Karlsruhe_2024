#include "IAIO.hpp"

#include "StringUtil.hpp"

#include "IACInterval.hpp"

void SetIAOutputPrecision(int width, int digits) { std::cout << cxsc::SetPrecision(width, digits); }

std::string OutputDown(double real) {
  std::string outDown("[");
  outDown << cxsc::SaveOpt << cxsc::RndDown << cxsc::real(real) << cxsc::RestoreOpt;
  outDown += "]";
  return outDown;
}

std::string OutputUp(double real) {
  std::string outUp("[");
  outUp << cxsc::SaveOpt << cxsc::RndUp << cxsc::real(real) << cxsc::RestoreOpt;
  outUp += "]";
  return outUp;
}

std::string OutputDown(double real, int width, int digits) {
  std::string outDown("[");
  outDown << cxsc::SaveOpt << cxsc::RndDown << cxsc::SetPrecision(width, digits) << cxsc::real(real)
          << cxsc::RestoreOpt;
  outDown += "]";
  return outDown;
}

std::string OutputUp(double real, int width, int digits) {
  std::string outUp("[");
  outUp << cxsc::SaveOpt << cxsc::RndUp << cxsc::SetPrecision(width, digits) << cxsc::real(real)
        << cxsc::RestoreOpt;
  outUp += "]";
  return outUp;
}

std::string Output(const IAInterval &value, int width, int digits) {
  std::string out("");
  out << cxsc::SaveOpt << cxsc::SetPrecision(width, digits)
      << cxsc::interval(value.inf(), value.sup()) << cxsc::RestoreOpt;
  return out;
}

std::string Output(const IACInterval &value, int width, int digits) {
  std::string out("");
  out << cxsc::SaveOpt << cxsc::SetPrecision(width, digits)
      << cxsc::cinterval(cxsc::interval(value.real().inf(), value.real().sup()),
                         cxsc::interval(value.imag().inf(), value.imag().sup()))
      << cxsc::RestoreOpt;
  return out;
}

std::string to_string(const IAInterval &value) {
  std::string out("");
  out << cxsc::interval(value.inf(), value.sup());
  return out;
}

std::string to_string(const IACInterval &value) {
  std::string out("");
  out << cxsc::cinterval(cxsc::interval(value.real().inf(), value.real().sup()),
                         cxsc::interval(value.imag().inf(), value.imag().sup()));
  return out;
}

std::string LatexOutputUp(double value, int numOfDigits) {
  std::string value_string = OutputUp(value, numOfDigits + 5, numOfDigits);
  value_string = value_string.substr(1, value_string.size() - 2);
  trim(value_string);
  std::vector<std::string> value_string_split = split(value_string, ".");
  return "$" + value_string_split[0] + "{.}" + value_string_split[1] + "$";
}

std::string LatexOutputDown(double value, int numOfDigits) {
  std::string value_string = OutputDown(value, numOfDigits + 5, numOfDigits);
  value_string = value_string.substr(1, value_string.size() - 2);
  trim(value_string);
  std::vector<std::string> value_string_split = split(value_string, ".");
  return "$" + value_string_split[0] + "{.}" + value_string_split[1] + "$";
}

std::string LatexOutput(const IAInterval &value, int numOfDigits) {
  std::string value_up = LatexOutputUp(sup(value), numOfDigits);
  value_up = value_up.substr(1, value_up.size() - 2);
  std::string value_down = LatexOutputDown(inf(value), numOfDigits);
  value_down = value_down.substr(1, value_down.size() - 2);

  int index;
  std::string output = "$";
  for (index = 0; index < value_up.size(); ++index) {
    if (value_up[index] != value_down[index]) break;
    output += value_up[index];
  }
  if (index < value_up.size()) {
    output += "_{" + value_down.substr(index) + "}^{" + value_up.substr(index) + "}";
  }
  output += "$";
  return output;
}