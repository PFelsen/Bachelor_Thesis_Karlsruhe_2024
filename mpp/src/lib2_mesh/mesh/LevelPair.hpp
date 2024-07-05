#ifndef LEVELPAIR_HPP
#define LEVELPAIR_HPP

#include <compare>
#include <cstddef>
#include <functional>
#include <ostream>
#include <sstream>
#include <string>

struct LevelPair {
  int space = -1;
  int time = -1;
  int adaptivityLevel = 0;
  int commSplit = 0;

  LevelPair NextInAdaptivity() const {
    return LevelPair{space, time, adaptivityLevel + 1, commSplit};
  };

  LevelPair CoarserInAdaptivity() const {
    return LevelPair{space, time, adaptivityLevel - 1, commSplit};
  };

  LevelPair WithoutAdaptivity() const { return LevelPair{space, time, 0, commSplit}; };

  LevelPair NextInSpace() const { return LevelPair{space + 1, time, adaptivityLevel, commSplit}; }

  LevelPair CoarserInSpace() const {
    return LevelPair{space - 1, time, adaptivityLevel, commSplit};
  }

  LevelPair NextInTime() const { return LevelPair{space, time + 1, adaptivityLevel, commSplit}; }

  LevelPair CoarserInTime() const { return LevelPair{space, time - 1, adaptivityLevel}; }

  static LevelPair Next(LevelPair from, LevelPair to) {
    LevelPair difference = {to.space - from.space, to.time - from.time,
                            to.adaptivityLevel - from.adaptivityLevel};
    return {from.space + (difference.space > 0), from.time + (difference.time > 0),
            from.adaptivityLevel + (difference.adaptivityLevel > 0)};
  }

  LevelPair CoarserInBoth() const { return LevelPair{space - 1, time - 1, commSplit}; }

  LevelPair WithCommSplit(int _commSplit) const {
    return LevelPair{space, time, adaptivityLevel, _commSplit};
  }

  std::string str() const {
    std::stringstream ss;
    ss << "[" << space << ", " << time << ", " << adaptivityLevel << ", " << commSplit << "]";
    return ss.str();
  }

  std::strong_ordering operator<=>(const LevelPair &other) const {
    if (space != other.space) { return space <=> other.space; }
    if (time != other.time) { return time <=> other.time; }
    if (adaptivityLevel != other.adaptivityLevel) {
      return adaptivityLevel <=> other.adaptivityLevel;
    }
    return std::strong_ordering::equivalent;
  }
};

inline bool operator==(const LevelPair &l1, const LevelPair &l2) {
  return (l1.space == l2.space) && (l1.time == l2.time)
         && (l1.adaptivityLevel == l2.adaptivityLevel) && (l1.commSplit == l2.commSplit);
}

inline bool operator!=(const LevelPair &l1, const LevelPair &l2) { return !(l1 == l2); }

namespace std {

template<>
struct hash<LevelPair> {
  size_t operator()(const LevelPair &levels) const {
    return size_t(levels.space * 10000 + levels.time + levels.adaptivityLevel * 10
                  + 100 * levels.commSplit);
  }
};

} // namespace std

inline std::string to_string(const LevelPair &level) {
  if (level.time < 0) return "[ " + std::to_string(level.space) + " ]";
  if (level.adaptivityLevel > 0) {
    return "[ " + std::to_string(level.space) + " | " + std::to_string(level.time) + " | "
           + std::to_string(level.adaptivityLevel) + " ]";
  }
  return "[ " + std::to_string(level.space) + " | " + std::to_string(level.time) + " ]";
}

inline std::ostream &operator<<(std::ostream &os, LevelPair level) {
  return os << to_string(level);
}

#endif // LEVELPAIR_HPP
