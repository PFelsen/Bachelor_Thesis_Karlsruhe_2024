#ifndef FUNCTOOLS_HPP
#define FUNCTOOLS_HPP

#include <algorithm>
#include <string>
#include <tuple>
#include <vector>

namespace ft {

template<typename... Ts>
struct Overload : Ts... {
  using Ts::operator()...;
};
template<class... Ts>
Overload(Ts...) -> Overload<Ts...>;

template<typename Cont, typename Pred>
Cont filter(const Cont &container, Pred predicate) {
  Cont result;
  std::copy_if(container.begin(), container.end(), std::back_inserter(result), predicate);
  return result;
}

const auto stod = [](const std::string &s) { return std::stod(s); };

const auto create_unit_vec = [](int size, int i) -> std::vector<double> {
  std::vector<double> vec(size, 0.0);
  vec[i] = 1.0;
  return vec;
};

const auto make_unit_vec = [](std::vector<double> &vec, int i) -> void {
  std::fill(begin(vec), end(vec), 0.0);
  vec[i] = 1.0;
};

struct T4IntHash {
  size_t operator()(std::tuple<int, int, int, int> const &t) const noexcept {
    return 100 * std::get<0>(t) + 10000 * std::get<1>(t) + 1000000 * std::get<2>(t)
           + 100000000 * std::get<3>(t);
  }
};

struct PairIntHash {
  std::size_t operator()(const std::pair<int, int> &k) const { return k.first + 1000 * k.second; }
};

std::vector<double> linspace(double init, double end, unsigned int size);

std::pair<double, double> linearFit(std::vector<double> &x, std::vector<double> &y);
} // namespace ft

#endif // FUNCTOOLS_HPP
