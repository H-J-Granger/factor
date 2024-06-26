#ifndef FACTOR_UTILITY_H_
#define FACTOR_UTILITY_H_

#include <cassert>
#include <chrono>
#include <random>

#include <gmpxx.h>

namespace factor {

constexpr long ayachi_constant = 0x0d000721;

std::mt19937_64 random_engine =
    std::mt19937_64(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now().time_since_epoch())
                        .count());

inline long long random_integer(long long l, long long r) {
  assert(l <= r);
  return random_engine() % (r - l + 1) + l;
}

class gmp_random_engine_t {
 public:
  gmp_randclass engine = gmp_randclass(gmp_randinit_mt);
  gmp_random_engine_t() {
    engine.seed(std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now().time_since_epoch())
                    .count() ^
                ayachi_constant);
  }
  mpz_class random_integer(mpz_class const& l, mpz_class const& r) {
    assert(l <= r);
    return engine.get_z_range(r - l) + l;
  }
} gmp_random_engine;

/**
 * @brief Return 2 if n is definitely prime, return 1 if n is probably prime 
 * (without being certain), or return 0 if n is definitely non-prime.
 */
inline int is_prime(mpz_class const& n) {
  return mpz_probab_prime_p(n.get_mpz_t(), 20);
}

inline mpz_class gcd(mpz_class const& n, mpz_class const& m) {
  mpz_class res;
  mpz_gcd(res.get_mpz_t(), n.get_mpz_t(), m.get_mpz_t());
  return res;
}

}  // namespace factor

#endif  // #ifndef FACTOR_UTILITY_H_