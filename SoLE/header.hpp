#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double f_x(double x) {
  return 3 * pow(x, 2) - 2 * x + 3;
}

double u_x(double x) { return x * (1 - x); }

double p_x(double x) { return 1 + pow(x, 2); }

double q_x(double x) { return 1 + x; }

vector<double> func_vec(int n) {
  double h = 1. / n;
  vector<double> b(n + 1, 0);
  for (int i = 1; i <= n; ++i) {
    b[i - 1] = f_x(i * h) * pow(h, 2);
  }
  return b;
}

vector<double> calculate_mtx_vec_mult(vector<vector<double>> &A,
                                      vector<double> &x) {
  int n = A.size();
  int m = x.size();
  vector<double> result(n, 0.0);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      result[i] += A[i][j] * x[j];
    }
  }

  return result;
}

vector<double> calculate_r(vector<vector<double>> &A, vector<double> &b,
                           vector<double> &x) {
  vector<double> Ax = calculate_mtx_vec_mult(A, x);
  vector<double> r(Ax.size());

  for (size_t i = 0; i < r.size(); ++i) {
    r[i] = Ax[i] - b[i];
  }

  return r;
}

double calculate_error_dec(vector<double> &x, vector<vector<double>> &A,
                           vector<double> &b) {
  vector<double> r = calculate_r(A, b, x);
  double max_err = 0.0;
  for (int i = 0; i < r.size(); ++i) {
    if (abs(r[i]) > max_err) max_err = abs(r[i]);
  }

  return max_err;
}
double calculate_error(vector<double> &new_x, vector<double> &old_x) {
  double max_err = 0.0;
  for (int i = 0; i < old_x.size(); ++i) {
    double err = abs(abs(new_x[i]) - abs(old_x[i]));
    if (err > max_err) max_err = err;
  }

  return max_err;
}

vector<vector<double>> create_matrix(int n) {
  double h = 1. / n;

  vector<vector<double>> matrix_res(n + 1, vector<double>(n + 1, 0.0));

  for (int i = 0; i < n; ++i) {
    if (i == 0) {
      double b = (p_x((i + 1) * h) + p_x((i + 2) * h) +
                  (pow(h, 2) * q_x((i + 1) * h)));
      double c = -p_x((i + 2) * h);
      matrix_res[i][i] = b;
      matrix_res[i][i + 1] = c;
      continue;
    }

    if (i == (n - 1)) {
      double b = (p_x((i + 1) * h) + p_x((i + 2) * h) +
                  (pow(h, 2) * q_x((i + 1) * h)));
      double a = -p_x((i + 1) * h);
      matrix_res[i][i] = b;
      matrix_res[i][i - 1] = a;
      continue;
    }

    double a = -p_x((i + 1) * h);
    double b =
        (p_x((i + 1) * h) + p_x((i + 2) * h) + (pow(h, 2) * q_x((i + 1) * h)));
    double c = -p_x((i + 2) * h);

    matrix_res[i][i] = b;
    matrix_res[i][i + 1] = c;
    matrix_res[i][i - 1] = a;
  }

  return matrix_res;
}

//Алгорит Томаса
void progonka_method(vector<vector<double>> &A, vector<double> &x, int n,
                     vector<double> &b) {
  vector<double> alpha(n + 1), betta(n + 1);

  double h = 1.0 / n;

  // прямой ход
  alpha[0] = A[0][1] / A[0][0];
  betta[0] = (b[0]) / A[0][0];

  for (int i = 1; i < n; ++i) {
    double del = 1.0 / (A[i][i] - alpha[i - 1] * A[i][i - 1]);
    alpha[i] = A[i][i + 1] * del;
    betta[i] = (-A[i][i - 1] * betta[i - 1] + b[i]) * del;
  }

  // обратный ход
  x[n - 1] = betta[n - 1];
  for (int i = n - 2; i >= 0; --i) {
    x[i] = -alpha[i] * x[i + 1] + betta[i];
  }

  for (int i = 1; i <= n; ++i) {
    printf(
        "ih = %4.2lf | y_i = %9.6lf | u(ih) = %8.6lf | |y_i - u(ih)| = "
        "%8.6lf\n",
        i * h, x[i - 1], u_x(i * h), abs(x[i - 1] - u_x(i * h)));
  }
}

double calculate_new_x(int i, vector<double> &x, int n,
                       vector<vector<double>> &A, vector<double> &b) {
  double h = 1.0 / n;
  double sum = 0.0;

  if (i > 0) {
    for (int j = 0; j <= i - 1; ++j) {
      sum += A[i][j] * x[j];
    }
  }

  if (i < n - 1) {
    for (int j = i + 1; j < n + 1; ++j) {
      sum += A[i][j] * x[j];
    }
  }

  return (b[i] - sum) * (1.0 / A[i][i]);
}

int Seidel_method(int n, vector<double> &x, vector<vector<double>> &A,
                  vector<double> &b) {
  double h = 1.0 / n;
  int k = 0;
  vector<double> new_x(n, 0.);
  double error = 1.0;

  while (error > 1.0 / pow(n, 3)) {
    for (int i = 0; i < n; ++i) {
      new_x[i] = calculate_new_x(i, new_x, n, A, b);
    }

    error = calculate_error_dec(new_x, A, b);

    x = new_x;

    k++;
  }

  return k;
}
double calculate_new_x_relax(int i, vector<double> &x, int n,
                             vector<vector<double>> &A, double omega,
                             vector<double> &b) {
  double sum = 0.0;

  if (i > 0) {
    for (int j = 0; j <= i - 1; ++j) {
      sum += A[i][j] * x[j];
    }
  }

  if (i < n - 1) {
    for (int j = i + 1; j < n + 1; ++j) {
      sum += A[i][j] * x[j];
    }
  }

  double new_x = (b[i] - sum) * (1.0 / A[i][i]);
  return x[i] + omega * (new_x - x[i]);
}

int relax_top(int n, vector<double> &x, vector<vector<double>> &A,
                 vector<double> &b) {
  double h = 1.0 / n;
  double omega = 1.9;
    double error = 1.0;
    int k = 0;
    vector<double> new_x(n, 0.);

    while (error > 1.0 / pow(n, 3)) {
      for (int i = 0; i < n; ++i) {
        new_x[i] = calculate_new_x_relax(i, new_x, n, A, omega, b);
      }

      error = calculate_error_dec(new_x, A, b);

      x = new_x;

      k++;
    }

    return k;
}
double calculate_tau(vector<double> &r, vector<double> &Ar) {
  double a = 0.;
  double b = 0.;
  for (size_t i = 0; i < r.size(); ++i) {
    a += r[i] * r[i];
    b += Ar[i] * r[i];
  }
  if (abs(b) < 1e-10) {
    return 0.0;
  }
  return a * (1.0 / b);
}

double calculate_new_x_spusk(int i, vector<double> &x, int n,
                             vector<vector<double>> &A, vector<double> &b) {
  vector<double> r = calculate_r(A, b, x);
  vector<double> Ar = calculate_mtx_vec_mult(A, r);

  return x[i] - calculate_tau(r, Ar) * r[i];
}

int spusk(int n, vector<double> &x, vector<vector<double>> &A,
          vector<double> &b) {
  double h = 1.0 / n;
  int k = 1;
  vector<double> new_x(n, 0.);
  double error = 1.0;

  while (error >= 1.0 / pow(n, 3)) {
    for (int i = 0; i < n; ++i) {
      new_x[i] = calculate_new_x_spusk(i, new_x, n, A, b);
    }

    x = new_x;

    error = calculate_error_dec(new_x, A, b);

    k++;
  }

  return k;
}

void m_print(int k, vector<double> &y_i, vector<double> &y_ik) {
  for (int i = 0; i < y_i.size() - 1; ++i) {
    printf(
        "ih = %3d | y_i = %9.6lf | y_ik = %9.6lf | |y_i - y_ik| = %9.6lf | k = "
        "%d\n",
        i, y_i[i], y_ik[i], abs(y_i[i] - y_ik[i]), k);
  }
}