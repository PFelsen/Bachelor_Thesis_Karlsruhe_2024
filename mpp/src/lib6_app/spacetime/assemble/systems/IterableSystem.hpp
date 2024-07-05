#ifndef __ITERABLE__SYSTEMS
#define __ITERABLE__SYSTEMS

#include <concepts>
#include <cstddef>
#include <functional>
#include <tuple>

#include <stdint.h>

template <bool outer = false>
struct OutputIterator {
  size_t i = 0;
  size_t j = 0;
  Scalar *entries;
  size_t offset;

  OutputIterator(DGRowEntries *entries)
      : entries(entries->a), offset(entries->nf) {}

  template <bool otherIterator>
  OutputIterator(const OutputIterator<otherIterator> &iterator, size_t init)
      : entries(iterator.entries),
        offset(iterator.offset),
        i(iterator.i),
        j(iterator.j) {
    if constexpr (outer) {
      i = init;
    } else {
      j = init;
    }
  }

  inline auto &operator*() const { return entries[i * offset + j]; };
};

struct VectorIterator {
  size_t i = 0;
  Scalar *rowIterator;

  VectorIterator(Vector &vector, const Point &point) : rowIterator(vector(point)) {}

  inline auto &operator*() const { return rowIterator[i]; };
};

template <std::same_as<DGRowEntries>... RowEntries>
inline auto zipOutput(RowEntries *...entries) {
  return std::make_tuple(OutputIterator<false>(entries)...);
}

template <bool outer>
inline bool operator!=(const OutputIterator<outer> &output1,
                       const OutputIterator<outer> &output2) {
  return output1.i != output2.i || output1.j != output2.j;
}

inline auto operator++(VectorIterator &output) {
  output.i++;
  return output;
}

inline bool operator!=(const VectorIterator &output1,
                       const VectorIterator &output2) {
  return output1.i != output2.i || output1.rowIterator != output2.rowIterator;
}

template <bool outer>
inline auto operator++(OutputIterator<outer> &output) {
  if constexpr (outer) {
    output.i++;
  } else {
    output.j++;
  }
  return output;
}

template <bool outer, std::integral Comparable>
inline bool operator<(const OutputIterator<outer> &output,
                      Comparable comparable) {
  if constexpr (outer) {
    return output.i < comparable;
  } else {
    return output.j < comparable;
  }
}

inline bool operator<(const VectorIterator &output,
                      const std::integral auto comparable) {
  return output.i < comparable;
}

template <bool outer>
inline auto &operator+=(OutputIterator<outer> &output, int64_t size) {
  if constexpr (outer) {
    output.i += size;
  } else {
    output.j += size;
  }
  return output;
}

template <typename Tuple>
struct CalculateIterator {
  Tuple types;
  size_t end;

  inline const auto operator*() const {
    return std::apply(
        [](const auto size, const auto &...inputs) {
          return std::make_tuple(size, (*inputs)...);
        },
        types);
  }
};

OutputIterator<false> transposeOutputIterator(
    const OutputIterator<true> &iterator) {
  OutputIterator<false> output(iterator, 0);
  output.i = iterator.i;
  return output;
}

OutputIterator<true> transposeOutputIterator(
    const OutputIterator<false> &iterator) {
  return OutputIterator<true>(iterator, iterator.j);
}

auto transposeOutputIterator(const std::integral auto &ref) {
  return ref;
}

auto transposeOutputIterator(const VectorIterator& ref) {
  return ref;
}

auto transposeOutputIterator(const TupleLike auto &ref) {
  return std::apply(
      [](const auto &...values) {
        return std::make_tuple(transposeOutputIterator(values)...);
      },
      ref);
}

template <typename Tuple>
struct CalculateIteratorStepper : public CalculateIterator<Tuple> {
  size_t stepper = 1;
};

inline auto &operator++(TupleLike auto &tuple) {
  std::apply([&](auto &...inputs) { (++inputs, ...); }, tuple);
  return tuple;
}

inline bool operator<(const TupleLike auto &output,
                      std::integral auto comparable) {
  return std::get<0>(output) < comparable;
}

template <typename Tuple>
inline auto &operator++(CalculateIteratorStepper<Tuple> &iterator) {
  return ++iterator.types;
}

template <typename Tuple>
inline auto &operator++(CalculateIterator<Tuple> &iterator) {
  do {
    ++iterator.types;
  } while (std::apply(
      [end = iterator.end](auto size, auto &...inputs) {
        bool flag = size < end;
        ((flag = flag && (norm(*inputs)) == 0), ...);
        return flag;
      },
      iterator.types));
  return iterator;
}

template <typename Tuple>
inline bool operator!=(const CalculateIterator<Tuple> &iterator,
                       const size_t end) {
  return std::get<0>(iterator.types) < end;
}

template <typename Tuple>
struct CalculateIterable {
  CalculateIterator<Tuple> beginTypes;
  const size_t endSize;

  CalculateIterable(const Tuple &beginTypes, const size_t endSize)
      : beginTypes(CalculateIteratorStepper<Tuple>{beginTypes, endSize}),
        endSize(endSize) {}

  inline const auto &begin() const { return beginTypes; }

  inline const auto end() const { return endSize; }
};

template <typename Tuple>
struct CalculateIterableStepper {
  CalculateIteratorStepper<Tuple> beginTypes;
  const size_t endSize;

  CalculateIterableStepper(const Tuple &beginTypes, const size_t endSize)
      : beginTypes(CalculateIteratorStepper<Tuple>{beginTypes, endSize}),
        endSize(endSize) {}

  inline const auto &begin() const { return beginTypes; }

  inline const auto end() const { return endSize; }
};

template <class... ExternalTypes>
inline auto make_iterable(const size_t size, const ExternalTypes &...types) {
  return make_iterable<ExternalTypes...>((size_t)0, size, types...);
}

template <bool useStepper = false, class... ExternalTypes>
inline auto make_iterable(const size_t startOffset, const size_t endSize,
                          const ExternalTypes &...types) {
  using IteratorType = std::remove_cvref_t<std::remove_const_t<
      std::tuple_element_t<0, std::tuple<ExternalTypes...>>>>;
  const auto valueTuple = std::invoke(
      [&](const auto &size, const auto &...inputs) {
        IteratorType iterator = size;
        return std::make_tuple(iterator, (std::begin(inputs) + startOffset)...);
      },
      types...);
  if constexpr (useStepper) {
    return CalculateIterableStepper(valueTuple, endSize);
  } else {
    return CalculateIterable(valueTuple, endSize);
  }
}

#define DEFAULT_OVERLOAD(name)                                        \
  template <bool useStepper = false, class CalculationInfo>           \
  inline auto name(const size_t startOf, CalculationInfo &info) {     \
    return name<useStepper>(startOf, info,                            \
                            OutputIterator<false>(&info.rowEntries)); \
  }

template <class CalculationInfo>
inline auto quadratureIterable(const CalculationInfo &info) {
  return make_iterable<true>(0, info.quadratureSize, size_t{0},
                             info.weightCache);
}

template <class CalculationInfo>
inline auto quadratureIterableWithPoint(const CalculationInfo &info) {
  return make_iterable<true>(0, info.quadratureSize, size_t{0},
                             info.weightCache, info.element.qPoint);
}

template <bool useStepper = false, class Iterator, class CalculationInfo>
inline auto pressureIterable(const size_t startOf, const CalculationInfo &info,
                             const Iterator &iterator) {
  const auto fullDimension = info.element.GetFullDimensions();
  return make_iterable<useStepper>(startOf * fullDimension, fullDimension,
                                   transposeOutputIterator(iterator),
                                   info.element.pressureST);
}
DEFAULT_OVERLOAD(pressureIterable);

template <bool useStepper = false, class Iterator, class CalculationInfo>
inline auto fluxIterable(const size_t startOf, const CalculationInfo &info,
                         const Iterator &iterator) {
  const auto fullDimension = info.element.GetFullDimensions();
  return make_iterable(startOf * fullDimension, fullDimension,
                       transposeOutputIterator(iterator), info.fluxCache);
}
DEFAULT_OVERLOAD(fluxIterable)

template <bool useStepper = false, class Iterator, class CalculationInfo>
inline auto derivativePressureIterable(const size_t startOf,
                                       const CalculationInfo &info,
                                       const Iterator &iterator) {
  const auto fullDimension = info.element.GetFullDimensions();
  return make_iterable<useStepper>(startOf * fullDimension, fullDimension,
                                   transposeOutputIterator(iterator),
                                   info.element.dtPressureST);
}
DEFAULT_OVERLOAD(derivativePressureIterable)

template <bool useStepper = false, class Iterator, class CalculationInfo>
inline auto derivativeVelocityIterable(const size_t startOf,
                                       const CalculationInfo &info,
                                       const Iterator &iterator) {
  const auto fullDimension = info.element.GetFullDimensions();
  return make_iterable<useStepper>(startOf * fullDimension, fullDimension,
                                   transposeOutputIterator(iterator),
                                   info.element.dtVelocityST);
}
DEFAULT_OVERLOAD(derivativeVelocityIterable)

template <bool useStepper = false, class Iterator, class CalculationInfo>
inline auto gradientPressureIterable(const size_t startOf,
                                     const CalculationInfo &info,
                                     const Iterator &iterator) {
  const auto fullDimension = info.element.GetFullDimensions();
  return make_iterable<useStepper>(startOf * fullDimension, fullDimension,
                                   transposeOutputIterator(iterator),
                                   info.element.gradPressureST);
}
DEFAULT_OVERLOAD(gradientPressureIterable)

template <bool useStepper = false, class Iterator, class CalculationInfo>
inline auto divergenceVelocityIterable(const size_t startOf,
                                       const CalculationInfo &info,
                                       const Iterator &iterator) {
  const auto fullDimension = info.element.GetFullDimensions();
  return make_iterable<useStepper>(startOf * fullDimension, fullDimension,
                                   transposeOutputIterator(iterator),
                                   info.element.divVelocityST);
}
DEFAULT_OVERLOAD(divergenceVelocityIterable)

template <bool useStepper = false, class Iterator, class CalculationInfo>
inline auto velocityIterable(const size_t startOf, const CalculationInfo &info,
                             const Iterator &iterator) {
  const auto fullDimension = info.element.GetFullDimensions();
  return make_iterable<useStepper>(startOf * fullDimension, fullDimension,
                                   transposeOutputIterator(iterator),
                                   info.element.velocityST);
}
DEFAULT_OVERLOAD(velocityIterable)

#endif