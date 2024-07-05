#ifndef _DEBUG_H_
#define _DEBUG_H_

#include <cassert>
#include <exception>
#include <string>
#include <typeinfo>
#include <vector>
#include "mpi.h"

template<class C>
inline const char *Type(const C &c) {
  return typeid(c).name();
}

template<class C>
inline const char *Type(const C *c) {
  return typeid(*c).name();
}

class Assertion {
  static std::vector<std::string> previousWarnings;
public:
  Assertion() = default;

  static void assertion(bool assertion, const std::string &string, const std::string &file,
                        size_t line);

  static void assertion(bool b, const char *s, int n, const char *f, int L);

  static void exit(const char *s, const char *f, int L);

  static void exit(std::string s, const char *f, int L);

  static void warning(const char *s, const char *f, int L);

  static void warning(std::string s, const char *f, int L);

  static void warningOnProc(int proc, const char *s, const char *f, int L);

  static void warningOnProc(int proc, std::string s, const char *f, int L);

  static void warningOnProc(const char *s, const char *f, int L);

  static void warningOnProc(std::string s, const char *f, int L);

  static void synchronizeErrors(bool failed, std::string = "", const char *f = "", int L = -1);

  static void Error(std::string, const char *f, int L);
};

struct MppException : public std::exception {
  std::string error;
  const char *file;
  int line;
  std::string errorMsg{};

  MppException(std::string error, const char *file, int line) :
      error(error), file(file), line(line) {
    errorMsg += error + "\n";
    errorMsg += "           in " + std::string(file) + ":" + std::to_string(line) + "\n";
    errorMsg += "           on line " + std::to_string(line) + "\n\n";
  }

  const char *what() const throw() { return errorMsg.c_str(); }
};

struct NotImplementedException : public std::exception {};

/// Use this for errors definitely occurring on all procs
#define ERROR(s)                                                                                   \
  {                                                                                                \
    Assertion::Error(s, __FILE__, __LINE__);                                                       \
    std::exit(1);                                                                                  \
  }

#define THROW(s) throw MppException(s, __FILE__, __LINE__);

#define TRY try

#define CATCH(s)                                                                                   \
  catch (const MppException &exception) {                                                          \
    Assertion::synchronizeErrors(true, exception.error, exception.file, exception.line);           \
  }                                                                                                \
  catch (...) {                                                                                    \
    Assertion::synchronizeErrors(true, s, __FILE__, __LINE__);                                     \
  }                                                                                                \
  Assertion::synchronizeErrors(false);

#define Exit(s)                                                                                    \
  {                                                                                                \
    Assertion::exit(s, __FILE__, __LINE__);                                                        \
    exit(1);                                                                                       \
  }

#define Warning(s)                                                                                 \
  { Assertion::warning(s, __FILE__, __LINE__); }

#define VerboseWarning(s, i)                                                                       \
  if (verbose >= i) Warning(s)

#define WarningOnProc(s)                                                                           \
  { Assertion::warningOnProc(s, __FILE__, __LINE__); }

#define WarningOnMaster(s)                                                                         \
  { Assertion::warningOnProc(0, s, __FILE__, __LINE__); }

#define LOCATION "file = " << __FILE__ << " line = " << __LINE__

#define Assert(s) Assertion::assertion(s, __STRING(s), __FILE__, __LINE__)

#define CheckMem(a)                                                                                \
  if (!(a)) Assertion::assertion(false, "out of memory", __FILE__, __LINE__);

#endif // of #ifndef _DEBUG_H_
