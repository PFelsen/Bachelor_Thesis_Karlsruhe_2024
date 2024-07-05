#ifndef BUFFER_HPP
#define BUFFER_HPP

#include <concepts>
#include <cstdio>
#include <cstring>
#include <string>
#include <tuple>
#include <vector>
#include "utility/Concepts.hpp"

class Buffer {
private:
  char *b;

  char *p;

  size_t n;

  size_t max_n;

  size_t BufferSize = 512000;
public:
  Buffer(size_t m = 0);

  ~Buffer();

  void Destruct();

  size_t size() const;

  size_t Size() const;

  void rewind();

  char *operator()();

  void resize(size_t m);

  template<typename T>
  Buffer &fill(const T &type, size_t m) {
    if (size() + m > n) {
      size_t j = size_t(((size() + m) - n) / BufferSize) + 1;
      resize(n + j * BufferSize);
    }
    memcpy(p, &type, m);
    p += m; // Todo double check purpose
    return *this;
  }

  template<typename T>
  Buffer &read(T &type, size_t m) {
    memcpy(&type, p, m);
    p += m; // Todo double check purpose
    return *this;
  }

  template<typename T>
  Buffer &operator<<(const T &type) {
    size_t m = sizeof(type);
    return fill(type, m);
  }

  template<typename T>
  Buffer &operator>>(T &type) {
    size_t m = sizeof(type);
    return read(type, m);
  }

  template<TupleLike TupleType>
  Buffer &operator<<(const TupleType &tuple) {
    std::apply([&buffer =
                    *this](const auto &...tupleArgs) mutable { ((buffer << tupleArgs), ...); },
               tuple);
    return *this;
  }

  template<TupleLike TupleType>
  Buffer &operator>>(TupleType &tuple) {
    std::apply([&buffer = *this](auto &...tupleArgs) mutable { ((buffer >> tupleArgs), ...); },
               tuple);
    return *this;
  }

  template<SizableIterableLike IterableLike>
  Buffer &operator<<(const IterableLike &iterable) {
    Buffer &buffer = *this;
    buffer << iterable.size();
    for (const auto &tuple : iterable) {
      buffer << tuple;
    }
    return buffer;
  }

  template<SizableIterableLike IterableLike>
  Buffer &operator>>(IterableLike &iterable) {
    Buffer &buffer = *this;
    size_t size = 0;
    buffer >> size;
    if constexpr (ReservableContainer<IterableLike>) { iterable.reserve(size); }
    for (size_t i = 0; i < size; i++) {
      typename TypeHolder<IterableLike>::value_type value;
      buffer >> value;
      iterable.insert(std::end(iterable), value);
    }
    return buffer;
  }
};

#endif // BUFFER_HPP
