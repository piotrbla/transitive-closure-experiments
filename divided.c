#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define ceild(n, d) ceil(((double)(n)) / ((double)(d)))
#define floord(n, d) floor(((double)(n)) / ((double)(d)))
#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x, y) ((x) < (y) ? (x) : (y))

#define S1(i, j, k) reach[i][j] = reach[i][j] || (reach[i][k] && reach[k][j])
void printMatrix(int **, int, int);
int **allocateMatrix(int);

void computeTC0(int **matrix, int n) {
  int **reach = allocateMatrix(n + 1);
  int i, j, k;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      reach[i][j] = matrix[i][j];
  double start = omp_get_wtime();
  register int lbp = 0;
  register int ubp = floord(n + 1, 2) - 1;
  for (int c0 = 0; c0 < ubp; c0 += 1)
    for (int c1 = 0; c1 <= min(c0 + 1, n / 2 - 1); c1 += 1)
      for (int c2 = 0; c2 <= min(c0 - c1 + 1, n / 2 - 1); c2 += 1)
        for (int c4 = max(2 * c0 + 1, 2 * c1 + 2 * c2); c4 <= 2 * c0 + 2;
             c4 += 1)
          for (int c5 = 2 * c1 + 1; c5 <= min(2 * c1 + 2, -2 * c2 + c4 + 1);
               c5 += 1)
            for (int c6 = 2 * c2 + 1; c6 <= min(2 * c2 + 2, c4 - c5 + 2);
                 c6 += 1)
              reach[c5][c6] = reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                                reach[c4 - c5 - c6 + 2][c6]);

  for (int c0 = ubp; c0 < 2 * n - 2; c0 += 1)
    for (int c1 = 0; c1 < n / 2; c1 += 1)
      for (int c2 = max(0, -2 * n + c0 - c1 + (n + 1) / 2);
           c2 <= min(c0 - c1 + 1, n / 2 - 1); c2 += 1) {
        if (3 * n + 2 * c1 == 2 * c0 + 2 && 2 * c2 + 2 == n)
          for (int c5 = -3 * n + 2 * c0 + 3; c5 <= -3 * n + 2 * c0 + 4; c5 += 1)
            reach[c5][n - 1] =
                reach[c5][n - 1] || (reach[c5][-3 * n + 2 * c0 - c5 + 4] &&
                                     reach[-3 * n + 2 * c0 - c5 + 4][n - 1]);
        for (int c4 = max(2 * c0 + 1, 2 * c1 + 2 * c2);
             c4 <= min(2 * c0 + 2, n + 2 * c1 + 2 * c2 - 1); c4 += 1)
          for (int c5 = 2 * c1 + 1;
               c5 <= min(min(n - 1, 2 * c1 + 2), -2 * c2 + c4 + 1); c5 += 1)
            for (int c6 = 2 * c2 + 1;
                 c6 <= min(min(n - 1, 2 * c2 + 2), c4 - c5 + 2); c6 += 1)
              reach[c5][c6] = reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                                reach[c4 - c5 - c6 + 2][c6]);
        if (2 * n + c2 >= c0 + c1 + 3) {
          for (int c4 = max(2 * c0 + 1, n + 2 * c1 + 2 * c2);
               c4 <= min(min(min(4 * n - 5, 2 * c0 + 2), 3 * n + 2 * c1 - 2),
                         n + c0 + c1 + c2 + 1);
               c4 += 1) {
            if (n + c0 + c1 + c2 >= c4) {
              for (int c5 = 2 * c1 + 1; c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
                for (int c6 = 2 * c2 + 1; c6 <= min(n - 1, 2 * c2 + 2);
                     c6 += 1) {
                  if (n + c5 + c6 >= c4 + 3) {
                    reach[c5][c6] =
                        reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                          reach[c4 - c5 - c6 + 2][c6]);
                  } else
                    reach[c5][c6] =
                        reach[c5][c6] || (reach[c5][-n + c4 - c5 - c6 + 2] &&
                                          reach[-n + c4 - c5 - c6 + 2][c6]);
                }
            } else
              for (int c5 = 2 * c1 + 1; c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
                for (int c6 = 2 * c2 + 1; c6 <= min(n - 1, 2 * c2 + 2);
                     c6 += 1) {
                  if (n + c5 + c6 >= c0 + c1 + c2 + 4) {
                    reach[c5][c6] = reach[c5][c6] ||
                                    (reach[c5][c0 + c1 + c2 - c5 - c6 + 3] &&
                                     reach[c0 + c1 + c2 - c5 - c6 + 3][c6]);
                  } else
                    reach[c5][c6] =
                        reach[c5][c6] ||
                        (reach[c5][-n + c0 + c1 + c2 - c5 - c6 + 3] &&
                         reach[-n + c0 + c1 + c2 - c5 - c6 + 3][c6]);
                }
          }
          if (c0 + 3 == 2 * n && 2 * c1 + 2 == n && 2 * c2 + 2 == n)
            reach[n - 1][n - 1] =
                reach[n - 1][n - 1] || (reach[n - 1][0] && reach[0][n - 1]);
        } else if (2 * c1 + 2 == n && 3 * n + 2 * c2 == 2 * c0 + 2)
          for (int c6 = -3 * n + 2 * c0 + 3; c6 <= -3 * n + 2 * c0 + 4; c6 += 1)
            reach[n - 1][c6] =
                reach[n - 1][c6] || (reach[n - 1][-3 * n + 2 * c0 - c6 + 4] &&
                                     reach[-3 * n + 2 * c0 - c6 + 4][c6]);
        for (int c4 = max(2 * c0 + 1, 2 * n + 2 * c1 + 2 * c2 + 2);
             c4 <= min(2 * c0 + 2, 3 * n + 2 * c1 + 2 * c2 + 1); c4 += 1)
          for (int c5 = max(2 * c1 + 1, -3 * n - 2 * c2 + c4 + 1);
               c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
            for (int c6 = max(2 * c2 + 1, -3 * n + c4 - c5 + 3);
                 c6 <= min(n - 1, 2 * c2 + 2); c6 += 1)
              reach[c5][c6] =
                  reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                    reach[-2 * n + c4 - c5 - c6 + 2][c6]);
      }

  for (int c0 = 2 * n - 2; c0 < 2 * n + floord(n, 2) - 2; c0 += 1)
    for (int c1 = -2 * n + c0 + 2; c1 < n / 2; c1 += 1)
      for (int c2 = max(-2 * n + c0 + 2, -2 * n + c0 - c1 + (n + 1) / 2);
           c2 < n / 2; c2 += 1)
        for (int c4 = 2 * c0 + 1;
             c4 <= min(min(5 * n - 5, 2 * c0 + 2), 3 * n + 2 * c1 + 2 * c2 + 1);
             c4 += 1)
          for (int c5 = max(max(2 * c1 + 1, -3 * n - 2 * c2 + c4 + 1),
                            -4 * n + c4 + 4);
               c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
            for (int c6 = max(2 * c2 + 1, -3 * n + c4 - c5 + 3);
                 c6 <= min(n - 1, 2 * c2 + 2); c6 += 1)
              reach[c5][c6] =
                  reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                    reach[-2 * n + c4 - c5 - c6 + 2][c6]);

  printf("Sequential: %lf\n", omp_get_wtime() - start);
  printMatrix(reach, n, 0);
}

void computeTC1(int **matrix, int n) {
  int **reach = allocateMatrix(n + 1);
  int i, j, k;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      reach[i][j] = matrix[i][j];
  double start = omp_get_wtime();
  register int lbp = 0;
  register int ubp = floord(n + 1, 2) - 1;
  //#pragma omp parallel for private(c1, c2, c3, c4, c5, c6)
  for (int c0 = 0; c0 < ubp; c0 += 1)
    for (int c1 = 0; c1 <= min(c0 + 1, n / 2 - 1); c1 += 1)
      for (int c2 = 0; c2 <= min(c0 - c1 + 1, n / 2 - 1); c2 += 1)
        for (int c4 = max(2 * c0 + 1, 2 * c1 + 2 * c2); c4 <= 2 * c0 + 2;
             c4 += 1)
#pragma omp parallel for private(c5, c6)
          for (int c5 = 2 * c1 + 1; c5 <= min(2 * c1 + 2, -2 * c2 + c4 + 1);
               c5 += 1)
            for (int c6 = 2 * c2 + 1; c6 <= min(2 * c2 + 2, c4 - c5 + 2);
                 c6 += 1)
              reach[c5][c6] = reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                                reach[c4 - c5 - c6 + 2][c6]);

  //#pragma omp parallel for private(c1, c2, c3, c4, c5, c6)
  for (int c0 = ubp; c0 < 2 * n - 2; c0 += 1)
    for (int c1 = 0; c1 < n / 2; c1 += 1)
      for (int c2 = max(0, -2 * n + c0 - c1 + (n + 1) / 2);
           c2 <= min(c0 - c1 + 1, n / 2 - 1); c2 += 1) {
        if (3 * n + 2 * c1 == 2 * c0 + 2 && 2 * c2 + 2 == n)
          for (int c5 = -3 * n + 2 * c0 + 3; c5 <= -3 * n + 2 * c0 + 4; c5 += 1)
            reach[c5][n - 1] =
                reach[c5][n - 1] || (reach[c5][-3 * n + 2 * c0 - c5 + 4] &&
                                     reach[-3 * n + 2 * c0 - c5 + 4][n - 1]);
        for (int c4 = max(2 * c0 + 1, 2 * c1 + 2 * c2);
             c4 <= min(2 * c0 + 2, n + 2 * c1 + 2 * c2 - 1); c4 += 1)
#pragma omp parallel for private(c5, c6)
          for (int c5 = 2 * c1 + 1;
               c5 <= min(min(n - 1, 2 * c1 + 2), -2 * c2 + c4 + 1); c5 += 1)
            for (int c6 = 2 * c2 + 1;
                 c6 <= min(min(n - 1, 2 * c2 + 2), c4 - c5 + 2); c6 += 1)
              reach[c5][c6] = reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                                reach[c4 - c5 - c6 + 2][c6]);
        if (2 * n + c2 >= c0 + c1 + 3) {
          for (int c4 = max(2 * c0 + 1, n + 2 * c1 + 2 * c2);
               c4 <= min(min(min(4 * n - 5, 2 * c0 + 2), 3 * n + 2 * c1 - 2),
                         n + c0 + c1 + c2 + 1);
               c4 += 1) {
            if (n + c0 + c1 + c2 >= c4) {
#pragma omp parallel for private(c5, c6)
              for (int c5 = 2 * c1 + 1; c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
                for (int c6 = 2 * c2 + 1; c6 <= min(n - 1, 2 * c2 + 2);
                     c6 += 1) {
                  if (n + c5 + c6 >= c4 + 3) {
                    reach[c5][c6] =
                        reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                          reach[c4 - c5 - c6 + 2][c6]);
                  } else
                    reach[c5][c6] =
                        reach[c5][c6] || (reach[c5][-n + c4 - c5 - c6 + 2] &&
                                          reach[-n + c4 - c5 - c6 + 2][c6]);
                }
            } else
#pragma omp parallel for private(c5, c6)
              for (int c5 = 2 * c1 + 1; c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
                for (int c6 = 2 * c2 + 1; c6 <= min(n - 1, 2 * c2 + 2);
                     c6 += 1) {
                  if (n + c5 + c6 >= c0 + c1 + c2 + 4) {
                    reach[c5][c6] = reach[c5][c6] ||
                                    (reach[c5][c0 + c1 + c2 - c5 - c6 + 3] &&
                                     reach[c0 + c1 + c2 - c5 - c6 + 3][c6]);
                  } else
                    reach[c5][c6] =
                        reach[c5][c6] ||
                        (reach[c5][-n + c0 + c1 + c2 - c5 - c6 + 3] &&
                         reach[-n + c0 + c1 + c2 - c5 - c6 + 3][c6]);
                }
          }
          if (c0 + 3 == 2 * n && 2 * c1 + 2 == n && 2 * c2 + 2 == n)
            reach[n - 1][n - 1] =
                reach[n - 1][n - 1] || (reach[n - 1][0] && reach[0][n - 1]);
        } else if (2 * c1 + 2 == n && 3 * n + 2 * c2 == 2 * c0 + 2)
          for (int c6 = -3 * n + 2 * c0 + 3; c6 <= -3 * n + 2 * c0 + 4; c6 += 1)
            reach[n - 1][c6] =
                reach[n - 1][c6] || (reach[n - 1][-3 * n + 2 * c0 - c6 + 4] &&
                                     reach[-3 * n + 2 * c0 - c6 + 4][c6]);
        for (int c4 = max(2 * c0 + 1, 2 * n + 2 * c1 + 2 * c2 + 2);
             c4 <= min(2 * c0 + 2, 3 * n + 2 * c1 + 2 * c2 + 1); c4 += 1)
#pragma omp parallel for private(c5, c6)
          for (int c5 = max(2 * c1 + 1, -3 * n - 2 * c2 + c4 + 1);
               c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
            for (int c6 = max(2 * c2 + 1, -3 * n + c4 - c5 + 3);
                 c6 <= min(n - 1, 2 * c2 + 2); c6 += 1)
              reach[c5][c6] =
                  reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                    reach[-2 * n + c4 - c5 - c6 + 2][c6]);
      }

  //#pragma omp parallel for private(c1, c2, c3, c4, c5, c6)
  for (int c0 = 2 * n - 2; c0 < 2 * n + floord(n, 2) - 2; c0 += 1)
    for (int c1 = -2 * n + c0 + 2; c1 < n / 2; c1 += 1)
      for (int c2 = max(-2 * n + c0 + 2, -2 * n + c0 - c1 + (n + 1) / 2);
           c2 < n / 2; c2 += 1)
        for (int c4 = 2 * c0 + 1;
             c4 <= min(min(5 * n - 5, 2 * c0 + 2), 3 * n + 2 * c1 + 2 * c2 + 1);
             c4 += 1)
#pragma omp parallel for private(c5, c6)
          for (int c5 = max(max(2 * c1 + 1, -3 * n - 2 * c2 + c4 + 1),
                            -4 * n + c4 + 4);
               c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
            for (int c6 = max(2 * c2 + 1, -3 * n + c4 - c5 + 3);
                 c6 <= min(n - 1, 2 * c2 + 2); c6 += 1)
              reach[c5][c6] =
                  reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                    reach[-2 * n + c4 - c5 - c6 + 2][c6]);

  printf("Parallel: %lf\n", omp_get_wtime() - start);
  printMatrix(reach, n, 1);
}
#define TILESIZE 32
void computeTC2(int **matrix, int n) {

  int **reach = allocateMatrix(n + 1);
  int i, j, k;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      reach[i][j] = matrix[i][j];
  double start = omp_get_wtime();
  register int lbp = 0;
  register int ubp = floord(n + 1, 32) - 1;
  if (n >= 16 && n <= 31) {
    for (int c4 = 1; c4 <= 32; c4 += 1)
      for (int c5 = 1; c5 <= min(n, c4); c5 += 1)
        for (int c6 = 1; c6 <= min(n, c4 - c5 + 1); c6 += 1) {
          if (c4 + 1 >= n + c5 + c6) {
            reach[c5][c6] =
                reach[c5][c6] || (reach[c5][-n + c4 - c5 - c6 + 2] &&
                                  reach[-n + c4 - c5 - c6 + 2][c6]);
          } else
            reach[c5][c6] = reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                              reach[c4 - c5 - c6 + 2][c6]);
        }
  } else {
    if (n >= 32)
      for (int c0 = 0; c0 < n / 16 + n / 32; c0 += 1)
        for (int c1 = 0; c1 <= min(c0, (n - 1) / 32); c1 += 1) {
          for (int c2 = 0; c2 <= min(c0 - c1, n / 32 - 1); c2 += 1)
            for (int c4 = 32 * c0 + 1; c4 <= 32 * c0 + 32; c4 += 1)
              for (int c5 = 32 * c1 + 1;
                   c5 <= min(min(n, 32 * c1 + 32), -32 * c2 + c4); c5 += 1)
                for (int c6 = 32 * c2 + 1; c6 <= min(32 * c2 + 32, c4 - c5 + 1);
                     c6 += 1) {
                  if (n + c5 + c6 >= c4 + 2) {
                    reach[c5][c6] =
                        reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                          reach[c4 - c5 - c6 + 2][c6]);
                  } else if (c4 + 1 >= 2 * n + c5 + c6) {
                    reach[c5][c6] = reach[c5][c6] ||
                                    (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                     reach[-2 * n + c4 - c5 - c6 + 2][c6]);
                  } else
                    reach[c5][c6] =
                        reach[c5][c6] || (reach[c5][-n + c4 - c5 - c6 + 2] &&
                                          reach[-n + c4 - c5 - c6 + 2][c6]);
                }
          if (32 * c0 + 31 >= n + 32 * c1 && n % 32 >= 1) {
            for (int c4 = 32 * c0 + 1;
                 c4 <= min(32 * c0 + 32, -((n - 1) % 32) + 2 * n + 32 * c1 - 1);
                 c4 += 1)
              for (int c5 = 32 * c1 + 1;
                   c5 <= min(min(n, 32 * c1 + 32), ((n - 1) % 32) - n + c4 + 1);
                   c5 += 1)
                for (int c6 = -((n - 1) % 32) + n; c6 <= min(n, c4 - c5 + 1);
                     c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                        reach[c4 - c5 - c6 + 2][c6]);
            for (int c4 = max(32 * c0 + 1, -((n - 1) % 32) + 2 * n + 32 * c1);
                 c4 <= 32 * c0 + 32; c4 += 1) {
              for (int c5 = 32 * c1 + 1;
                   c5 <= min(32 * c1 + 32, ((n - 1) % 32) - 2 * n + c4 + 1);
                   c5 += 1)
                for (int c6 = -((n - 1) % 32) + n; c6 <= n; c6 += 1) {
                  if (n + c5 + c6 >= c4 + 2) {
                    reach[c5][c6] =
                        reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                          reach[c4 - c5 - c6 + 2][c6]);
                  } else
                    reach[c5][c6] =
                        reach[c5][c6] || (reach[c5][-n + c4 - c5 - c6 + 2] &&
                                          reach[-n + c4 - c5 - c6 + 2][c6]);
                }
              for (int c5 = ((n - 1) % 32) - 2 * n + c4 + 2;
                   c5 <= min(n, 32 * c1 + 32); c5 += 1)
                for (int c6 = -((n - 1) % 32) + n; c6 <= n; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                        reach[c4 - c5 - c6 + 2][c6]);
            }
          }
        }
    if (n % 32 == 0)
      for (int c0 = 3 * n / 32; c0 < 5 * n / 32; c0 += 1)
        for (int c1 = max(0, (-n / 8) + c0); c1 < n / 32; c1 += 1)
          for (int c2 = max(0, (-3 * n / 32) + c0 - c1 - 1); c2 < n / 32;
               c2 += 1) {
            for (int c4 = 32 * c0 + 1;
                 c4 <= min(32 * c0 + 32, 2 * n + 32 * c1 + 32 * c2 + 62);
                 c4 += 1) {
              for (int c5 = 32 * c1 + 1; c5 < -2 * n - 32 * c2 + c4 - 30;
                   c5 += 1)
                for (int c6 = max(32 * c2 + 1, -3 * n + c4 - c5 + 2);
                     c6 <= 32 * c2 + 32; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                        reach[-2 * n + c4 - c5 - c6 + 2][c6]);
              for (int c5 = max(32 * c1 + 1, -2 * n - 32 * c2 + c4 - 30);
                   c5 <= 32 * c1 + 32; c5 += 1) {
                for (int c6 = 32 * c2 + 1; c6 <= -2 * n + c4 - c5 + 1; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                        reach[-2 * n + c4 - c5 - c6 + 2][c6]);
                for (int c6 = max(32 * c2 + 1, -2 * n + c4 - c5 + 2);
                     c6 <= 32 * c2 + 32; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][-n + c4 - c5 - c6 + 2] &&
                                        reach[-n + c4 - c5 - c6 + 2][c6]);
              }
            }
            for (int c4 = max(32 * c0 + 1, 2 * n + 32 * c1 + 32 * c2 + 63);
                 c4 <= min(32 * c0 + 32, 3 * n + 32 * c1 + 32 * c2 + 62);
                 c4 += 1)
              for (int c5 = max(32 * c1 + 1, -3 * n - 32 * c2 + c4 - 30);
                   c5 <= 32 * c1 + 32; c5 += 1)
                for (int c6 = max(32 * c2 + 1, -3 * n + c4 - c5 + 2);
                     c6 <= 32 * c2 + 32; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                        reach[-2 * n + c4 - c5 - c6 + 2][c6]);
          }
  }
  if (n >= 2 && n % 32 >= 1)
    for (int c0 = n / 16 + n / 32; c0 <= (5 * n - 3) / 32; c0 += 1) {
      if (n >= 8 * c0 + 1 && 32 * c0 + 33 >= 5 * n) {
        for (int c4 = 32 * c0 + 1; c4 < 5 * n - 1; c4 += 1)
          for (int c5 = max(1, -4 * n + c4 + 2); c5 <= min(n, c4); c5 += 1)
            for (int c6 = max(1, -3 * n + c4 - c5 + 2);
                 c6 <= min(n, c4 - c5 + 1); c6 += 1) {
              if (n + c5 + c6 >= c4 + 2) {
                reach[c5][c6] = reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                                  reach[c4 - c5 - c6 + 2][c6]);
              } else if (c4 + 1 >= 2 * n + c5 + c6) {
                reach[c5][c6] =
                    reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                      reach[-2 * n + c4 - c5 - c6 + 2][c6]);
              } else
                reach[c5][c6] =
                    reach[c5][c6] || (reach[c5][-n + c4 - c5 - c6 + 2] &&
                                      reach[-n + c4 - c5 - c6 + 2][c6]);
            }
      } else if (8 * c0 >= n) {
        for (int c1 = c0 - (n + 7) / 8; c1 <= (n - 1) / 32; c1 += 1) {
          for (int c2 = max(c0 - (n + 7) / 8, c0 - c1 - (3 * n + 29) / 32 - 1);
               c2 < n / 32; c2 += 1)
            for (int c4 = 32 * c0 + 1;
                 c4 <= min(min(32 * c0 + 32, 4 * n + 32 * c2 + 30),
                           3 * n + 32 * c1 + 32 * c2 + 62);
                 c4 += 1)
              for (int c5 = max(32 * c1 + 1, -3 * n - 32 * c2 + c4 - 30);
                   c5 <= min(n, 32 * c1 + 32); c5 += 1)
                for (int c6 = max(32 * c2 + 1, -3 * n + c4 - c5 + 2);
                     c6 <= 32 * c2 + 32; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                        reach[-2 * n + c4 - c5 - c6 + 2][c6]);
          for (int c4 = 32 * c0 + 1;
               c4 <= min(min(5 * n - 2, 32 * c0 + 32), 4 * n + 32 * c1 + 30);
               c4 += 1)
            for (int c5 = max(32 * c1 + 1, -4 * n + c4 + 2);
                 c5 <= min(n, 32 * c1 + 32); c5 += 1)
              for (int c6 = max(-3 * n + c4 - c5 + 2, -((n - 1) % 32) + n);
                   c6 <= n; c6 += 1)
                reach[c5][c6] =
                    reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                      reach[-2 * n + c4 - c5 - c6 + 2][c6]);
        }
      } else
        for (int c1 = 0; c1 <= floord(n - 1, 32); c1 += 1) {
          for (int c2 = max(0, c0 - c1 - (3 * n + 29) / 32 - 1); c2 < n / 32;
               c2 += 1)
            for (int c4 = 32 * c0 + 1;
                 c4 <= min(32 * c0 + 32, 3 * n + 32 * c1 + 32 * c2 + 62);
                 c4 += 1)
              for (int c5 = max(32 * c1 + 1, -3 * n - 32 * c2 + c4 - 30);
                   c5 <= min(n, 32 * c1 + 32); c5 += 1) {
                for (int c6 = max(32 * c2 + 1, -3 * n + c4 - c5 + 2);
                     c6 <= min(32 * c2 + 32, -n + c4 - c5 + 1); c6 += 1) {
                  if (2 * n + c5 + c6 >= c4 + 2) {
                    reach[c5][c6] =
                        reach[c5][c6] || (reach[c5][-n + c4 - c5 - c6 + 2] &&
                                          reach[-n + c4 - c5 - c6 + 2][c6]);
                  } else if (c0 % 2 == 0) {
                    reach[c5][c6] = reach[c5][c6] ||
                                    (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                     reach[-2 * n + c4 - c5 - c6 + 2][c6]);
                  } else if (c6 >= 33) {
                    reach[c5][c6] = reach[c5][c6] ||
                                    (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                     reach[-2 * n + c4 - c5 - c6 + 2][c6]);
                  } else if (n >= 64 && c4 >= 2 * n + 16 * c0 + 17) {
                    reach[c5][c6] = reach[c5][c6] ||
                                    (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                     reach[-2 * n + c4 - c5 - c6 + 2][c6]);
                  } else if (n <= 63 && c4 >= 2 * n + 16 * c0 + 17) {
                    reach[c5][c6] = reach[c5][c6] ||
                                    (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                     reach[-2 * n + c4 - c5 - c6 + 2][c6]);
                  } else
                    reach[c5][c6] = reach[c5][c6] ||
                                    (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                     reach[-2 * n + c4 - c5 - c6 + 2][c6]);
                }
                for (int c6 = -n + c4 - c5 + 2; c6 <= 32 * c2 + 32; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                        reach[c4 - c5 - c6 + 2][c6]);
              }
          if (3 * c1 == c0)
            for (int c4 = 32 * c0 + 1; c4 <= (64 * c0 / 3) + n; c4 += 1)
              for (int c5 = (32 * c0 / 3) + 1;
                   c5 <= min(n, (-32 * c0 / 3) + c4); c5 += 1)
                for (int c6 = (32 * c0 / 3) + 1; c6 <= min(n, c4 - c5 + 1);
                     c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                        reach[c4 - c5 - c6 + 2][c6]);
          for (int c4 = max(32 * c0 + 1, -((n - 1) % 32) + 2 * n + 32 * c1);
               c4 <= min(min(4 * n - 2, 32 * c0 + 32), 3 * n + 32 * c1 + 30);
               c4 += 1) {
            for (int c5 = 32 * c1 + 1; c5 <= -3 * n + c4 + 1; c5 += 1) {
              if (n >= 32) {
                for (int c6 = -((n - 1) % 32) + n; c6 <= n; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                        reach[-2 * n + c4 - c5 - c6 + 2][c6]);
              } else if (n >= 16 && 2 * n + 32 >= c4) {
                for (int c6 = -3 * n + c4 - c5 + 2; c6 <= n; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                        reach[-2 * n + c4 - c5 - c6 + 2][c6]);
              } else if (c4 >= 3 * n + 1 &&
                         16 * ((-2 * n + c4 - 1) / 32) + 15 >= n) {
                for (int c6 = -3 * n + c4 - c5 + 2; c6 <= n; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                        reach[-2 * n + c4 - c5 - c6 + 2][c6]);
              } else
                for (int c6 = 1; c6 <= n; c6 += 1)
                  reach[1][c6] = reach[1][c6] || (reach[1][n - c6 + 1] &&
                                                  reach[n - c6 + 1][c6]);
            }
            for (int c5 = max(32 * c1 + 1, -3 * n + c4 + 2);
                 c5 <=
                 min(min(n, 32 * c1 + 32), ((n - 1) % 32) - 2 * n + c4 + 1);
                 c5 += 1)
              for (int c6 = -((n - 1) % 32) + n; c6 <= n; c6 += 1) {
                if (n + c5 + c6 >= c4 + 2) {
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                        reach[c4 - c5 - c6 + 2][c6]);
                } else if (c4 + 1 >= 2 * n + c5 + c6) {
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                        reach[-2 * n + c4 - c5 - c6 + 2][c6]);
                } else
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][-n + c4 - c5 - c6 + 2] &&
                                        reach[-n + c4 - c5 - c6 + 2][c6]);
              }
            for (int c5 = ((n - 1) % 32) - 2 * n + c4 + 2;
                 c5 <= min(n, 32 * c1 + 32); c5 += 1) {
              if (c4 + 32 >= ((2 * n + c4 - 1) % 32) + 3 * n) {
                for (int c6 = -((n - 1) % 32) + n; c6 <= min(n, c4 - c5 + 1);
                     c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                        reach[c4 - c5 - c6 + 2][c6]);
              } else if (n >= 32) {
                for (int c6 = -((n - 1) % 32) + n; c6 <= n; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                        reach[c4 - c5 - c6 + 2][c6]);
              } else if (2 * n >= c4 + 2) {
                for (int c6 = 1; c6 <= c4 - c5 + 1; c6 += 1)
                  reach[c5][c6] =
                      reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] &&
                                        reach[c4 - c5 - c6 + 2][c6]);
              } else
                for (int c6 = 1; c6 <= n; c6 += 1)
                  reach[n][c6] = reach[n][c6] || (reach[n][n - c6 + 1] &&
                                                  reach[n - c6 + 1][c6]);
            }
          }
          for (int c4 = max(32 * c0 + 1, 3 * n + 32 * c1 + 31);
               c4 <= min(4 * n - 2, 32 * c0 + 32); c4 += 1)
            for (int c5 = 32 * c1 + 1; c5 <= 32 * c1 + 32; c5 += 1)
              for (int c6 = max(-3 * n + c4 - c5 + 2, -((n - 1) % 32) + n);
                   c6 <= n; c6 += 1)
                reach[c5][c6] =
                    reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                      reach[-2 * n + c4 - c5 - c6 + 2][c6]);
          for (int c4 = 4 * n - 1; c4 <= 32 * c0 + 32; c4 += 1)
            for (int c5 = max(32 * c1 + 1, -4 * n + c4 + 2);
                 c5 <= min(n, 32 * c1 + 32); c5 += 1)
              for (int c6 = max(-3 * n + c4 - c5 + 2, -((n - 1) % 32) + n);
                   c6 <= n; c6 += 1)
                reach[c5][c6] =
                    reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] &&
                                      reach[-2 * n + c4 - c5 - c6 + 2][c6]);
        }
    }

printf("Parallel2: %lf\n", omp_get_wtime() - start);
printMatrix(reach, n, 2);
}

void computeTC3(int **matrix, int n) {
  int **reach = allocateMatrix(n + 1);
  int i, j, k;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      reach[i][j] = matrix[i][j];
  double start = omp_get_wtime();
  for (k = 0; k < n; k++) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        reach[i][j] = reach[i][j] || (reach[i][k] && reach[k][j]);
      }
    }
  }
  printf("Seq Normal: %lf\n", omp_get_wtime() - start);
  printMatrix(reach, n, 2);
}

void printMatrix(int **matrix, int N, int fileno) {
  /*
      for (int i = 0; i < N; i++)
      {
          for (int j = 0; j < N; j++)
              printf ("%d ", matrix[i][j]);
          printf("\n");
      }
  */
  char filename[10];
  sprintf(filename, "result%d", fileno);
  FILE *f = fopen(filename, "wt");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++)
      fprintf(f, "%d ", matrix[i][j]);
    fprintf(f, "\n");
  }
  fclose(f);
}
int **allocateMatrix(int N) {
  int **t = (int **)malloc(sizeof(int *) * N);
  for (int i = 0; i < N; i++) {
    t[i] = (int *)malloc(sizeof(int) * N);
  }
  return t;
}
int main(void) {
  const int N = 400;
  int **graph = allocateMatrix(N);
  int g[4][4] = {{1, 1, 0, 1}, {0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      graph[i][j] = g[i][j];
  for (int i = 0; i < N; i++)
    graph[i][i] = 1;

  printMatrix(graph, 6, 9);
  computeTC0(graph, N);
  computeTC1(graph, N);
  computeTC2(graph, N);
  computeTC3(graph, N);
  return 0;
}
