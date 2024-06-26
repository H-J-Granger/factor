#ifndef FACTOR_PRIME_H_
#define FACTOR_PRIME_H_

#include <chrono>
#include <iostream>

#include <gmpxx.h>

#include "utility.h"

namespace factor {

/**
 * @brief Performs a Miller-Rabin test.
 * TODO: Use Montgomery reduction or Barret's reduction to improve speed, but 
 * since this is not the bottleneck, the improvement is not really neccesary.
 * 
 * @param n The number to identify primality.
 * @param base The base.
 * 
 * @return true if the number is a psuedo-prime;
 * @return false if it's a composite number.
 */
bool miller_rabin(mpz_class const& n, mpz_class const& base) {
  int s = mpz_scan1(n.get_mpz_t(), 0);
  std::cout << s << std::endl;
  mpz_class d = n >> s;
  mpz_class now = base;
  if (now > n) [[unlikely]] {
    now %= n;
  }

  for (int k = 1; k <= s; ++k) {
    if ((k == 1 && now == 1) || (k > 1 && now == n - 1)) [[unlikely]] {
      return true;
    }
    now = now * now % n;
  }
  return false;
}

/**
 * @brief Performs multiple Miller-Rabin test, on different random bases in 
 * [1,n]. Each test have a false positive rate of 25%, so multiple tests will
 * reduce it to 4^(-times).
 * 
 * @param n The number to identify primality.
 * @param times The number of times doing Miller-Rabin test.
 * 
 * @return true if the number is a psuedo-prime;
 * @return false if it's a composite number.
 */
bool multiple_miller_rabin(mpz_class const& n, int times) {
  for (int i = 0; i < times; ++i) {
    if (!miller_rabin(n, gmp_random_engine.random_integer(1, n))) {
      return false;
    }
  }
  return true;
}

}  // namespace factor

#endif  // #ifndef FACTOR_PRIME_H_