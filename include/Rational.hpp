#ifndef JKQZX_INCLUDE_RATIONAL_HPP_
#define JKQZX_INCLUDE_RATIONAL_HPP_

#include <stdint.h>

namespace zx {

/*
 * Representation of fractions as multiples of pi
 * Rationals can only have values in the half-open interval (-1,1],
 * corresponding to the interval (-pi, pi]
*/
class Rational {
  void normalize();

public:
  int32_t num, denom;
  
  Rational() : num(0), denom(1){};
  explicit Rational(int32_t num, int32_t denom) : num(num), denom(denom){normalize();}
  explicit Rational(int32_t num) : num(num), denom(1){normalize();}
  explicit Rational(double val);

  Rational &operator+=(const Rational &rhs);
  Rational &operator+=(const int32_t rhs);

  Rational &operator-=(const Rational &rhs);
  Rational &operator-=(const int32_t rhs);

  Rational &operator*=(const Rational &rhs);
  Rational &operator*=(const int32_t rhs);

  Rational &operator/=(const Rational &rhs);
  Rational &operator/=(const int32_t rhs);

  double to_double() const;
  bool is_integer() const {return denom == 1;}
};

inline Rational operator-(const Rational &rhs) {
  return Rational(-rhs.num, rhs.denom);
}
inline Rational operator+(Rational lhs, const Rational &rhs) {
  lhs += rhs;
  return lhs;
}
inline Rational operator+(Rational lhs, const int32_t rhs) {
  lhs += rhs;
  return lhs;
}
inline Rational operator+(const int32_t lhs, Rational rhs) {
  rhs += lhs;
  return rhs;
}

inline Rational operator-(Rational lhs, const Rational &rhs) {
  lhs -= rhs;
  return lhs;
}
inline Rational operator-(Rational lhs, const int32_t rhs) {
  lhs -= rhs;
  return lhs;
}
inline Rational operator-(const int32_t lhs, Rational rhs) {
  rhs -= lhs;
  return rhs;
}

inline Rational operator*(Rational lhs, const Rational &rhs) {
  lhs *= rhs;
  return lhs;
}
inline Rational operator*(Rational lhs, const int32_t rhs) {
  lhs *= rhs;
  return lhs;
}
inline Rational operator*(const int32_t lhs, Rational rhs) {
  rhs *= lhs;
  return rhs;
}

inline Rational operator/(Rational lhs, const Rational &rhs) {
  lhs /= rhs;
  return lhs;
}
inline Rational operator/(Rational lhs, const int32_t rhs) {
  lhs /= rhs;
  return lhs;
}
inline Rational operator/(const int32_t lhs, Rational rhs) {
  rhs /= lhs;
  return rhs;
}

inline bool operator<(const Rational &lhs, const Rational &rhs) {
  return lhs.num * rhs.denom < rhs.num * lhs.denom;
}

inline bool operator<(const Rational &lhs, int32_t rhs) {
  return lhs.num < rhs * lhs.denom;
}

inline bool operator<(int32_t lhs, const Rational &rhs) {
  return lhs * rhs.denom < rhs.num;
}

inline bool operator<=(const Rational &lhs, const Rational &rhs) {
  return lhs.num * rhs.denom <= rhs.num * lhs.denom;
}

inline bool operator<=(const Rational &lhs, int32_t rhs) {
  return lhs.num <= rhs * lhs.denom;
}

inline bool operator<=(int32_t lhs, const Rational &rhs) {
  return lhs * rhs.denom <= rhs.num;
}

inline bool operator>(const Rational &lhs, const Rational &rhs) {
  return rhs < lhs;
}

inline bool operator>(const Rational &lhs, int32_t rhs) { return rhs < lhs; }

inline bool operator>(int32_t lhs, const Rational &rhs) { return rhs < lhs; }

inline bool operator>=(const Rational &lhs, const Rational &rhs) {
  return rhs <= lhs;
}

inline bool operator>=(const Rational &lhs, int32_t rhs) { return rhs <= lhs; }

inline bool operator>=(int32_t lhs, const Rational &rhs) { return rhs <= lhs; }

  inline bool operator==(const Rational &lhs, const Rational &rhs) {
    return lhs.num == rhs.num && lhs.denom == rhs.denom;
}

inline bool operator==(const Rational &lhs, int32_t rhs) { return lhs.num == rhs && lhs.denom == 1; }

inline bool operator==(int32_t lhs, const Rational &rhs) { return rhs == lhs; }


} // namespace zx
#endif /* JKQZX_INCLUDE_RATIONAL_HPP_ */
