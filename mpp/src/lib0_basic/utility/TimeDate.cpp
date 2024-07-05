#include <climits>
#include <iostream>

#include "Parallel.hpp"
#include "TimeDate.hpp"


using namespace std;

bool nodate = false;

void NoDate() { nodate = true; }

ostream &operator<<(ostream &s, const Date &d) {
  if (nodate) return s << endl;
  char *p = ctime(&d.t);
  p[strcspn(p, "\n")] = '\0';
  return s << p;
}

ostream &operator<<(ostream &s, const Time &t) {
  int M = t.Minutes();
  double S = t.Seconds();
  char c[64];
  sprintf(c, "%5.2f", S + 0.005);
  if (M == 0) return s << c << " seconds";
  int H = t.Hours();
  S -= 60 * M;
  sprintf(c, "%d:%05.2f", M, S);
  if (H == 0) return s << c << " minutes";
  M -= 60 * H;
  sprintf(c, "%d:%02d:%05.2f", H, M, S);
  return s << c << " hours";
}

ostream &operator<<(ostream &s, const MTime &t) {
  int M = t.Minutes();
  double S = t.Seconds() + 0.005;
  char c[64];
  sprintf(c, "%d:%02d", M, int(S - 60 * M));
  return s << c << " minutes";
}

Time operator-(const Date &d2, const Date &d1) {
  const double maxclock = INT_MAX * 0.9999 / CLOCKS_PER_SEC;
  Time t;
  t.t = difftime(d2.t, d1.t);
  if (t.t > maxclock) return t;
  double d = double(d2.c - d1.c) / CLOCKS_PER_SEC;
  if (d < 0.0) d = double(INT_MAX + d2.c - d1.c) / CLOCKS_PER_SEC;
  if (abs(t.t - d) > 2.0) return t;
  t.t = d;
  return t;
}

bool operator<(const Time &t1, const Time &t2) {
  if (t1.Seconds() < t2.Seconds()) return true;
  return false;
}

void Time::Max() { t = PPM->Max(t, 0); }

void DateTime::SetMax(short l) {
  if (TimeLvl >= 1) elapsed = PPM->Max(elapsed, l);
}

void DateTime::SetMin(short l) {
  if (TimeLvl >= 1) elapsed = PPM->Min(elapsed, l);
}

double DateTime::GetMax(short l) {
  if (TimeLvl >= 1) return PPM->Max(elapsed, l);
  return 0;
}

double DateTime::GetMin(short l) {
  if (TimeLvl >= 1) return PPM->Min(elapsed, l);
  return 0;
}

double DateTime::SetSum(short l) {
  if (TimeLvl >= 1) {
    PPM->Sum(&elapsed, 1, l);
    return elapsed;
  }
  return 0;
}

void DateTime::AddTime() {
  if (TimeLvl >= 1) elapsed += (Date() - Start).Seconds();
}

Date startTime = Date();
