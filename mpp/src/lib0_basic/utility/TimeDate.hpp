#ifndef _TIMEDATE_H_
#define _TIMEDATE_H_

#include <string>
#include <vector>
#include <json.hpp>
class Time;

class Date {
  clock_t c;
  time_t t;
public:
  Date() {
    c = clock();
    t = time(0);
  }

  friend Time operator-(const Date &, const Date &);

  friend std::ostream &operator<<(std::ostream &, const Date &);

  friend void to_json(nlohmann::json &j, const Date &d);
};

std::ostream &operator<<(std::ostream &s, const Date &d);

class Time {
public:
  double t;

  Time() { t = 0; }

  Time(double tt) : t(tt) {}

  void Max();

  double Seconds() const { return t; }

  int Minutes() const { return int(t / 60.0); }

  int Hours() const { return int(t / 3600.0); }

  Time &operator=(const Time &T) {
    t = T.t;
    return *this;
  }

  Time &operator+=(const Time &T) {
    t += T.t;
    return *this;
  }

  Time &operator-=(const Time &T) {
    t -= T.t;
    return *this;
  }

  Time &operator/=(double n) {
    t /= n;
    return *this;
  }

  ~Time() = default;

  friend Time operator-(const Date &, const Date &);
};

class MTime : public Time {
public:
  MTime &operator=(const Time &T) {
    t = T.t;
    return *this;
  }
};

class Times {
  std::vector<Time> t{};
public:
  Times() = default;

  int size() const { return t.size(); }

  void AddTime(const Time &time) { t.push_back(time); }

  Time Average() {
    Time avg;
    for (int i = 0; i < t.size(); ++i)
      avg += t[i];
    return avg /= t.size();
  }
};

std::ostream &operator<<(std::ostream &s, const Time &t);

std::ostream &operator<<(std::ostream &s, const MTime &t);

Time operator-(const Date &d2, const Date &d1);

bool operator<(const Time &, const Time &);

void NoDate();

class DateTime {
  Date Start;
  double elapsed;
  std::string name;
  int TimeLvl;
public:
  DateTime() : TimeLvl(1) {}

  DateTime(std::string N) : name(N), TimeLvl(1) {}

  DateTime(const DateTime &DT) : Start(DT.Start), elapsed(DT.elapsed), name(DT.name), TimeLvl(1) {}

  void SetDate() { Start = Date(); }

  void SetTimeLvl(int t = 0) { TimeLvl = t; }

  void AddTime();

  void AddTime(double d) { elapsed += d; }

  double GetTime() const { return elapsed; }

  void SetName(std::string N) { name = N; }

  void ResetTime() { elapsed = 0; }

  std::string GetName() const { return name; }

  void SetMax(short l = -1);

  void SetMin(short l = -1);

  double GetMax(short l = -1);

  double GetMin(short k = -1);

  double SetSum(short l = -1);
};

extern Date startTime;

#endif // of #ifndef _TIMEDATE_H_
