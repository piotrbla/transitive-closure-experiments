#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S0(a, i, j, k) c[i][j] = c[i][k] + c[k][j]

void printMatrix(int **, int, int);
int **allocateMatrix(int);

void computeDYN0(int **matrix, int n) {
  int **c = allocateMatrix(n + 1);
  int i, j, k, a;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      c[i][j] = matrix[i][j];
  double start = omp_get_wtime();
  for (int c0 = 2; c0 < n; c0 += 1)
    for (int c1 = 1; c1 <= n - c0; c1 += 1)
      for (int c2 = c0 + c1; c2 <= min(n, 2 * c0 + c1 - 2); c2 += 1) {
	if (2 * c0 + c1 >= c2 + 3)
	  S0(c0, c1, c2, -c0 + c2 + 1);
	S0(c0, c1, c2, c0 + c1 - 1);
  }
  printf("nontiled: %lf\n", omp_get_wtime() - start);
  printMatrix(c, n, 0);
}

void computeDYN1(int **matrix, int n) {
  int **c = allocateMatrix(n + 1);
  int i, j, k, a;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      c[i][j] = matrix[i][j];
  double start = omp_get_wtime();
  for (int c0 = 2; c0 < n; c0 += 1)
    #pragma omp parallel for 
    for (int c1 = 1; c1 <= n - c0; c1 += 1)
      for (int c2 = c0 + c1; c2 <= min(n, 2 * c0 + c1 - 2); c2 += 1) {
        if (2 * c0 + c1 >= c2 + 3)
          S0(c0, c1, c2, -c0 + c2 + 1);
        S0(c0, c1, c2, c0 + c1 - 1);
      }
  printf("parallel: %lf\n", omp_get_wtime() - start);
  printMatrix(c, n, 1);
}

void printMatrix(int **matrix, int N, int fileno) {
  char filename[10];
  sprintf(filename, "nontiled%d", fileno);
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
  const int N = 2000;
  int **graph = allocateMatrix(N);
  int g[4][4] = {{1, 1, 0, 1}, {0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      graph[i][j] = g[i][j];
  for (int i = 0; i < N; i++)
    graph[i][i] = 1;

  printMatrix(graph, 6, 9);
  computeDYN0(graph, N);
  computeDYN1(graph, N);
  return 0;
}

