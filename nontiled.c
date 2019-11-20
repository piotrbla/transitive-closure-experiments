#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define ceild(n, d) ceil(((double)(n)) / ((double)(d)))
#define floord(n, d) floor(((double)(n)) / ((double)(d)))
#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x, y) ((x) < (y) ? (x) : (y))

#define S1(aaa, i, j, k) reach[i][j] = reach[i][j] || (reach[i][k] && reach[k][j])
void printMatrix(int **, int, int);
int **allocateMatrix(int);

void computeTC0(int **matrix, int n) {
  int **reach = allocateMatrix(n + 1);
  int i, j, k;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      reach[i][j] = matrix[i][j];
  double start = omp_get_wtime();
  for (int c0 = 1; c0 < 3 * n - 1; c0 += 1) {
    for (int c1 = 1; c1 <= -2 * n + c0 + 1; c1 += 1)
      for (int c2 = 1; c2 <= n; c2 += 1) {
        if (2 * n + c1 + c2 >= c0 + 2) {
          reach[c1][c2] = reach[c1][c2] || (reach[c1][-n + c0 - c1 - c2 + 2] && reach[-n + c0 - c1 - c2 + 2][c2]);
        } else
          reach[c1][c2] = reach[c1][c2] || (reach[c1][-2 * n + c0 - c1 - c2 + 2] && reach[-2 * n + c0 - c1 - c2 + 2][c2]);
      }
    for (int c1 = max(1, -2 * n + c0 + 2); c1 <= min(n, c0); c1 += 1)
      for (int c2 = 1; c2 <= min(n, c0 - c1 + 1); c2 += 1) {
        if (n + c1 + c2 >= c0 + 2) {
          reach[c1][c2] = reach[c1][c2] || (reach[c1][c0 - c1 - c2 + 2] && reach[c0 - c1 - c2 + 2][c2]);
        } else
          reach[c1][c2] = reach[c1][c2] || (reach[c1][-n + c0 - c1 - c2 + 2] && reach[-n + c0 - c1 - c2 + 2][c2]);
      }
  }
  for (int c0 = 3 * n - 1; c0 < 4 * n - 1; c0 += 1) {
    for (int c1 = 1; c1 <= -3 * n + c0 + 1; c1 += 1)
      for (int c2 = -3 * n + c0 - c1 + 2; c2 <= n; c2 += 1)
        reach[c1][c2] = reach[c1][c2] || (reach[c1][-2 * n + c0 - c1 - c2 + 2] && reach[-2 * n + c0 - c1 - c2 + 2][c2]);
    for (int c1 = -3 * n + c0 + 2; c1 <= n; c1 += 1)
      for (int c2 = 1; c2 <= n; c2 += 1) {
        if (2 * n + c1 + c2 >= c0 + 2) {
          reach[c1][c2] = reach[c1][c2] || (reach[c1][-n + c0 - c1 - c2 + 2] && reach[-n + c0 - c1 - c2 + 2][c2]);
        } else
          reach[c1][c2] = reach[c1][c2] || (reach[c1][-2 * n + c0 - c1 - c2 + 2] && reach[-2 * n + c0 - c1 - c2 + 2][c2]);
      }
  }
  for (int c0 = 4 * n - 1; c0 < 5 * n - 1; c0 += 1)
    for (int c1 = -4 * n + c0 + 2; c1 <= n; c1 += 1)
      for (int c2 = -3 * n + c0 - c1 + 2; c2 <= n; c2 += 1)
        reach[c1][c2] = reach[c1][c2] || (reach[c1][-2 * n + c0 - c1 - c2 + 2] && reach[-2 * n + c0 - c1 - c2 + 2][c2]);

  printf("Sequential: %lf\n", omp_get_wtime() - start);
  printMatrix(reach, n, 0);
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
  return 0;
}
