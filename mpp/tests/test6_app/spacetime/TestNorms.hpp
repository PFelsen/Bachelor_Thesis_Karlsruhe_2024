#ifndef SPACETIME_TESTNORMS_HPP
#define SPACETIME_TESTNORMS_HPP

#include "TestEnvironment.hpp"

class TestNorms : public TestWithParam<std::string> {
protected:
  TestNorms() {}
};

INSTANTIATE_TEST_SUITE_P(TestNorms, TestNorms, Values(""));


#endif // SPACETIME_TESTNORMS_HPP
