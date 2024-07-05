#ifndef SPACETIME_STSINGLERUNMAIN_HPP
#define SPACETIME_STSINGLERUNMAIN_HPP

#include "STMethod.hpp"


class STSingleRunMethod : public STMethod {
public:

  STSingleRunMethod() : STMethod(STMainBuilder().ReadConfig()) {};

  STSingleRunMethod(STMainBuilder builder);

  void run() override;
};

#endif //SPACETIME_STSINGLERUNMAIN_HPP
