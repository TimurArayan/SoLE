#include "header.hpp"

using namespace std;

int main() {
  int n = 50.0;

  std::vector<double> y_i_p(n + 1), y_i_s(n + 1), y_i_r(n + 1), y_i_sp(n + 1),
      b(n);
  vector<vector<double>> A = create_matrix(n);
  b = func_vec(n);
  printf("Progonka\n");
  progonka_method(A, y_i_p, n, b);
  int k_s = Seidel_method(n, y_i_s, A, b);
  int k_r = relax_top(n, y_i_r, A, b);
  int k_dec = spusk(n, y_i_sp, A, b);
  printf("Zeidel\n");
  m_print(k_s, y_i_p, y_i_s);
  printf("High Relaxation\n");
  m_print(k_r, y_i_p, y_i_r);
  printf("Spusk\n");
  m_print(k_dec, y_i_p, y_i_sp);
  system("pause");
}