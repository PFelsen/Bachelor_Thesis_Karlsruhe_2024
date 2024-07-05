#pragma once

#ifdef __cpp_concepts
#include <concepts>

namespace mpp::shape {
template<typename T>
concept PolyCalculationType = requires(T a, T b) {
  { a *b } -> std::convertible_to<T>;
  { a / b } -> std::convertible_to<T>;
  { a + b } -> std::convertible_to<T>;
  { a - b } -> std::convertible_to<T>;
  { a * 2.0 } -> std::convertible_to<T>;
  { 2.0 * a } -> std::convertible_to<T>;
  { a / 2.0 } -> std::convertible_to<T>;
  { 2.0 / a } -> std::convertible_to<T>;
  { a - 2.0 } -> std::convertible_to<T>;
  { 2.0 - a } -> std::convertible_to<T>;
  { 2.0 + a } -> std::convertible_to<T>;
  { a + 2.0 } -> std::convertible_to<T>;
};

} // namespace mpp::shape
#else
#define PolyCalculationType typename
#endif