#ifndef INTERVALARITHMETIC_H
#define INTERVALARITHMETIC_H

#include <string>

class IAInterval;

class IACInterval;

void SetIAOutputPrecision(int width = -1, int digits = -1);

std::string OutputDown(double real);

std::string OutputUp(double real);

std::string OutputDown(double real, int width, int digits);

std::string OutputUp(double real, int width, int digits);

std::string Output(const IAInterval &value, int width, int digits);

std::string Output(const IACInterval &value, int width, int digits);

std::string to_string(const IAInterval &value);

std::string to_string(const IACInterval &value);

std::string LatexOutputUp(double value, int numOfDigits);

std::string LatexOutputDown(double value, int numOfDigits);

std::string LatexOutput(const IAInterval &value, int numOfDigits);

#endif // INTERVALARITHMETIC_H
