#include "StatisticalTest.hpp"

#include "ChiSquared.hpp"
#include "KolmogorovSmirnov.hpp"
#include "ShapiroWilk.hpp"

/*
std::unique_ptr<StatisticalTest1D> CreateStatisticalTestUnique(std::string test_name){
  if(test_name == "ShapiroWilk") {
    return std::make_unique<ShapiroWilkTest>();
  }
  if(test_name == "KolmogorovSmirnov") {
    return std::make_unique<KolmogorovSmirnovTest>();
  }
  if(test_name == "ChiSquared") {
    return std::make_unique<ChiSquaredTest>();
  }
}
*/