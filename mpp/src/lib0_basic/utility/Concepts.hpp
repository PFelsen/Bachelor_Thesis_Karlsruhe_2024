#ifndef CONCEPTS_HPP
#define CONCEPTS_HPP

#include <concepts>

template<typename T>
concept ReservableContainer = requires(T type) { type.reserve(1); };

template<typename T>
concept SizableIterableLike = requires(T type, typename T::value_type value) {
  typename T::value_type;
  std::begin(type);
  std::end(type);
  type.size();
  type.insert(std::end(type), value);
};

template<typename T>
concept NotIterableLike = !SizableIterableLike<T>;

template<typename T>
concept Ostreamable = requires(const T &type, std::ostream &stream) { stream << type; };

template<typename T>
concept MapIterableLike = requires(T type, typename T::value_type value) {
  typename T::key_type;
  typename T::mapped_type;
} && SizableIterableLike<T>;

template<typename T>
concept TupleLike = requires(T t) {
  std::tuple_size<T>::value;
  std::get<0>(t);
  typename std::tuple_element<0, T>::type;
};

template<typename IterableLike>
struct TypeHolder {};

template<class T>
concept StringLike = std::is_convertible_v<T, std::string_view>;

template<SizableIterableLike IterableLike>
struct TypeHolder<IterableLike> {
  typedef typename IterableLike::value_type value_type;
};

template<MapIterableLike IterableLike>
struct TypeHolder<IterableLike> {
  typedef typename std::pair<typename IterableLike::key_type, typename IterableLike::mapped_type>
      value_type;
};

#endif