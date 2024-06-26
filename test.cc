#include "factor.h"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
  using namespace factor;
  std::ios::sync_with_stdio(false);
  std::cin.tie(0);

  mpz_class a;
  std::cin >> a;

  double x = std::sqrt(a.get_d());
  long search_limit =
      1.3 * std::pow(std::exp(std::sqrt(std::log(x) * std::log(std::log(x)))),
                     1 / std::sqrt(2));
  std::cout << lenstra_ecm(a, 12000, 10'000'000'000) << std::endl;
  return 0;
  std::cout << pollard_rho(a) << std::endl;

  try {
    int b = 0;
    mpz_class ans = 1;
    while (ans == 1 || ans == b) {
      b++;
      ans = shank_square_form(a, 100'000'000, b);
    }
    std::cout << ans << std::endl;
  } catch (char const* str) {
    std::cout << str << std::endl;
  }
}