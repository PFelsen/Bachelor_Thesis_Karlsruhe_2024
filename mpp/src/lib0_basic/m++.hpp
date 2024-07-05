#ifndef _MPLUSPLUS_HPP_
#define _MPLUSPLUS_HPP_

#include "Assertion.hpp"
#include "Config.hpp"

class Mpp {
public:
  static Mpp &initialize(int *argc, char **argv);

  static Mpp &instance();
private:
  static Mpp *mpp;

  Mpp(int *argc, char **argv);

  Mpp(const Mpp &);

  ~Mpp();

  class MppGuard {
  public:
    ~MppGuard();
  };
};

#endif // of #ifndef _MPLUSPLUS_HPP_
