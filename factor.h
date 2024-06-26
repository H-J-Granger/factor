#ifndef FACTOR_FACTOR_H_
#define FACTOR_FACTOR_H_

#ifdef FACTOR_DEBUG
#include <iostream>
#endif

#include <chrono>
#include <cmath>
#include <vector>

#include <gmpxx.h>

#include "utility.h"

namespace factor {

/**
 * @brief Performs trivial division for divisors under limit.
 * 
 * @return std::vector<int> all divisors found.
 */
std::vector<int> trivial_division(mpz_class n, int limit) {
  std::vector<int> res;
  // Scan 2
  int twos = mpz_scan1(n.get_mpz_t(), 0);
  n >>= twos;
  res.reserve(twos);
  for (int i = 1; i <= twos; ++i) {
    res.push_back(2);
  }
  // Scan 3
  while (n % 3 == 0) {
    n /= 3;
    res.push_back(3);
  }
  // Scan others
  int now = 5;
  int delta[2] = {2, 4};
  int index = 0;
  while (now <= limit && n != 1) {
    while (n % now == 0) {
      n /= now;
      res.push_back(now);
    }
    now += delta[index];
    index ^= 1;
  }
  return res;
}

/**
 * @brief Performs Pollard's rho algorithm.
 * 
 * @param n The number to factor.
 * @param failure_limit The limit of failure counts.
 * @return mpz_class The non-trival divisor found.
 */
mpz_class pollard_rho(mpz_class n, int failure_limit = 50) {
  assert(n >= 10);
  std::vector<int> ans = trivial_division(n, 10);
  if (!ans.empty()) {
    return ans[0];
  }
  for (int times = 1; times <= failure_limit; times++) {
    mpz_class a = gmp_random_engine.random_integer(1, n - 1);
    mpz_class x = gmp_random_engine.random_integer(1, n - 1);
    mpz_class y = (x * x + a) % n;
    mpz_class product = 1;
    int cnt = 1;
    while (true) {
      product = product * (y - x + n) % n;
      if (cnt == 20) {
        cnt = 1;
        mpz_class res = gcd(product, n);
        if (res != 1) {
          if (res == n) {
            break;
          } else {
            return res;
          }
        }
      }
      x = (x * x + a) % n;
      y = (y * y + a) % n;
      y = (y * y + a) % n;
      cnt++;
    }
  }
  throw "factor::pollard_rho: Failure limit reached.";
}

/**
 * @brief Performs Shanks' square form algorithm.
 * If process failed by returning 1, one can try using some multiples of n
 * to factor.
 * Runs extremely fast when n = pq, where prime numbers p, q are very close to 
 * sqrt(n).
 * 
 * @param n The number to factor.
 * @return mpz_class A non-trival factor.
 */
mpz_class shank_square_form(mpz_class n_, long test_limit = 100000,
                            long k = 1) {
  assert(n_ >= 10);
  if (mpz_probab_prime_p(n_.get_mpz_t(), 5)) {
    if (mpz_probab_prime_p(n_.get_mpz_t(), 30)) {
      throw "factor::shank_square_form: The number is prime.";
    }
  }
  mpz_class n = n_ * k;
  mpz_class p, prev_p, prev_q = 1, q;
  mpz_sqrtrem(prev_p.get_mpz_t(), q.get_mpz_t(), n.get_mpz_t());
  mpz_class sqrt_n = prev_p, b;
  if (q == 0) {
    return sqrt_n;
  }
  mpz_class tmp;
  bool success = false;
  for (long i = 1; i <= test_limit; ++i) {
    b = (sqrt_n + prev_p) / q;
    p = b * q - prev_p;
    tmp = q;
    q = prev_q + b * (prev_p - p);
    prev_q = tmp;
    prev_p = p;
    if (mpz_perfect_square_p(q.get_mpz_t())) {
      success = true;
      break;
    }
  }
  if (!success) {
    throw "factor::shank_square_form: The first loop reaches test limit.";
  }
  mpz_class sqrt_q_k;
  mpz_sqrt(sqrt_q_k.get_mpz_t(), q.get_mpz_t());
  b = (sqrt_n - p) / sqrt_q_k;
  p = prev_p = b * sqrt_q_k + p;
  prev_q = sqrt_q_k;
  q = (n - p * p) / prev_q;
  for (long i = 1; i <= test_limit; ++i) {
    b = (sqrt_n + prev_p) / q;
    p = b * q - prev_p;
    tmp = q;
    q = prev_q + b * (prev_p - p);
    prev_q = tmp;
    if (prev_p == p) {
      mpz_class gcd_res;
      mpz_gcd(gcd_res.get_mpz_t(), p.get_mpz_t(), n_.get_mpz_t());
      return gcd_res;
    }
    prev_p = p;
  }
  throw "factor::shank_square_form: The second loop reaches test limit.";
}

struct ec_point {
  mpz_class x, y;
};

ec_point ec_add(ec_point p, ec_point q, mpz_class const& n,
                mpz_class const& a) {
  mpz_class lambda = 0;
  if (p.x == q.x && p.y == q.y) {
    lambda = (3 * p.x * p.x + a) % n;
    mpz_class inverse;
    if (mpz_invert(inverse.get_mpz_t(), mpz_class(2 * p.y).get_mpz_t(),
                   n.get_mpz_t()) == 0) {
      throw p.y;
    }
    lambda = lambda * inverse % n;
  } else {
    lambda = p.y - q.y;
    if (lambda < 0) {
      lambda += n;
    }
    mpz_class inverse, deno = p.x - q.x;
    if (deno < 0) {
      deno += n;
    }
    if (mpz_invert(inverse.get_mpz_t(), deno.get_mpz_t(), n.get_mpz_t()) == 0) {
      throw deno;
    }
  }
  ec_point res;
  res.x = (lambda * lambda - p.x - q.x + n * 2) % n;
  res.y = (lambda * (p.x - res.x + n) - p.y + n) % n;
  return res;
}

std::vector<long> eratosthenes_sieve(long n) {
  std::vector<long> res;
  std::vector<bool> not_prime;
  not_prime.resize(n);
  not_prime[0] = not_prime[1] = true;
  for (int i = 2; i <= n; ++i) {
    if (!not_prime[i]) {
      res.push_back(i);
      if (static_cast<long long>(i) * i > n) {
        continue;
      }
      for (int j = i * i; j <= n; j += i) {
        not_prime[j] = true;
      }
    }
  }
  return res;
}

/**
 * @brief Performs Lenstra's Elliptic-Curve Method.
 * search_limit should be at least exp(sqrt(log(n)*log(log(n))))**(1/sqrt(2)).
 * 
 * @param n The number to factor.
 * @param search_limit 
 * @return mpz_class 
 */
mpz_class lenstra_ecm(mpz_class n, long search_limit = 12000,
                      long test_limit = 100) {
  if (n % 2 == 0) {
    return 2;
  }
  if (n % 3 == 0) {
    return 3;
  }
  for (int times = 1; times <= test_limit; ++times) {
    auto primes = eratosthenes_sieve(search_limit);
    search_limit += 100;
    try {
      mpz_class a = gmp_random_engine.random_integer(1, n - 1);
      ec_point point(0, 1);
      for (auto p : primes) {
        long multiple = 1;
        int power = std::floor(std::log(search_limit) / std::log(p));
        for (int i = 1; i <= power; ++i) {
          multiple *= p;
        }
        ec_point now(point.x, point.y);
        while (multiple) {
          if (multiple & 1) {
            point = ec_add(point, now, n, a);
          }
          now = ec_add(now, now, n, a);
          multiple >>= 1;
        }
      }
    } catch (mpz_class d) {
      mpz_class gcd_res = 0;
      mpz_gcd(gcd_res.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());
      if (gcd_res != n) {
        std::fprintf(stderr, "%d\n", times);
        return gcd_res;
      }
    }
  }
  throw "factor::lenstra_ecm: Test limit reached.";
}

}  // namespace factor

#endif  // #ifndef FACTOR_FACTOR_H_