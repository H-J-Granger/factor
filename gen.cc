#include <gmpxx.h>
#include <cstdio>
#include "utility.h"

char str[100000];

int main(int argc, char** argv) {
  using namespace factor;
  if (argc < 3) {
    return 0;
  }
  int n, m;
  std::sscanf(argv[1], "%d", &n);
  std::sscanf(argv[2], "%d", &m);
  mpz_class tmp;
  mpz_pow_ui(tmp.get_mpz_t(), mpz_class(10).get_mpz_t(), n);
  mpz_class p = gmp_random_engine.random_integer(tmp / 10 * 9, tmp / 10 * 11);
  mpz_class q;
  mpz_nextprime(q.get_mpz_t(), p.get_mpz_t());
  mpz_pow_ui(tmp.get_mpz_t(), mpz_class(10).get_mpz_t(), m);
  p = gmp_random_engine.random_integer(tmp / 10 * 9, tmp / 10 * 11);
  mpz_class q2;
  mpz_nextprime(q2.get_mpz_t(), p.get_mpz_t());
  std::FILE* out = std::fopen("tmp.txt", "w");
  mpz_out_str(out, 10, mpz_class(q * q2).get_mpz_t());
  std::fclose(out);
  std::system("./test < tmp.txt > out.txt");
  return 0;
}