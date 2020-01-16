#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S0(a, i, j, k) c[i][j] = c[i][k] + c[k][j]

void printMatrix(int**, int, int);
int** allocateMatrix(int);
void deallocateMatrix(int**, int);

void write_results(int , double , char );
void write_results(int , double );

void computeDYN0(int** matrix, int n) {
  int** c = allocateMatrix(n + 1);
  int** d = allocateMatrix(n + 1);
  int i, j;
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
  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(c, n, 0);
  deallocateMatrix(c, n + 1);
  deallocateMatrix(d, n + 1);
}

void computeDYN1(int** matrix, int n) {
  int** c = allocateMatrix(n + 1);
  int** d = allocateMatrix(n + 1);
  int i, j;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      c[i][j] = matrix[i][j];
  double start = omp_get_wtime();
  for (int c0 = 2; c0 < n; c0 += 1)
#pragma omp parallel for private(c1, c2, c0)
    for (int c1 = 1; c1 <= n - c0; c1 += 1)
      for (int c2 = c0 + c1; c2 <= min(n, 2 * c0 + c1 - 2); c2 += 1) {
        if (2 * c0 + c1 >= c2 + 3)
          S0(c0, c1, c2, -c0 + c2 + 1);
        S0(c0, c1, c2, c0 + c1 - 1);
      }
  double execution_time = omp_get_wtime() - start;
  printf("parallel: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(c, n, 1);
  deallocateMatrix(c, n + 1);
  deallocateMatrix(d, n + 1);
}

void computeDYN2(int** matrix, int n) {
  int** c = allocateMatrix(n + 1);
  int** d = allocateMatrix(n + 1);
  int i, j;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      c[i][j] = matrix[i][j];
  double start = omp_get_wtime();

for (int c0 = 2; c0 <= floord(9 * n - 14, 80) + 2; c0 += 1) {
  #pragma omp parallel for num_threads(6)
  for (int c1 = max(max(-3 * c0 + 4, -((n + 19) / 20)), -c0 - (n + 28) / 30 + 2); c1 <= min(-c0 + (n - 2) / 16 + 1, -((4 * c0 + 5) / 9)); c1 += 1) {
    for (int c2 = max(-c0 + (8 * c0 - 5 * c1 - 8) / 23 + 1, -c0 - c1 - (n - 16 * c0 - 16 * c1 + 61) / 46 + 1); c2 <= min(min(-1, -c0 - c1), -((15 * c0 + 5 * c1 + 38) / 31) + 1); c2 += 1) {
      for (int c4 = max(max(2, 15 * c0 + 5 * c1 + 15 * c2 - 23), -16 * c2 - 15); c4 <= min(min(min(min(n - 1, -20 * c1 - 1), -16 * c2), 30 * c0 + 10 * c1 + 30 * c2), n + 30 * c0 + 30 * c1 + 30 * c2); c4 += 1) {
        for (int c5 = max(max(1, -30 * c0 - 30 * c1 - 30 * c2), -20 * c1 - 2 * c4 - 17); c5 <= min(min(-30 * c0 - 30 * c1 - 30 * c2 + 29, -20 * c1 - c4), n - c4); c5 += 1) {
          for (int c6 = max(-20 * c1 - 19, c4 + c5); c6 <= min(min(n, -20 * c1), 2 * c4 + c5 - 2); c6 += 1) {
            if (2 * c4 + c5 >= c6 + 3) {
              c[c5][c6] = (c[c5][-c4 + c6 + 1] + c[-c4 + c6 + 1][c6]);
            }
            c[c5][c6] = (c[c5][c4 + c5 - 1] + c[c4 + c5 - 1][c6]);
          }
        }
      }
    }
  }
}



  double execution_time = omp_get_wtime() - start;
  printf("tiles: %lf\n", execution_time);
  write_results(n, execution_time, '\n');
  printMatrix(c, n, 2);
  deallocateMatrix(c, n + 1);
  deallocateMatrix(d, n + 1);
}

void printMatrix(int** matrix, int N, int fileno) {
  char filename[10];
  sprintf(filename, "nontiled%d", fileno);
  FILE* f= fopen(filename, "wt");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++)
      fprintf(f, "%d ", matrix[i][j]);
    fprintf(f, "\n");
  }
  fclose(f);
}

int** allocateMatrix(int N) {
  int** t = (int**)malloc(sizeof(int*) * N);
  for (int i = 0; i < N; i++) {
    t[i] = (int*)malloc(sizeof(int) * N);
  }
  return t;
}

void deallocateMatrix(int **t, int N) {
  for (int i = 0; i < N; i++) {
    free(t[i]);
  }
  free(t);
}

void write_results(int n, double execution_time, char end_char)
{
  FILE* f = fopen("results.txt", "at");
  fprintf(f, "%d:%lf%c", n, execution_time, end_char);
  fclose(f);
}
void write_results(int n, double execution_time)
{
  write_results(n, execution_time, ';');
}

int main(void) {
  const int ZMAX = 2010;
  int** graph = allocateMatrix(ZMAX);
  int g[4][4] = { {1, 1, 0, 1}, {0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1} };
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      graph[i][j] = g[i][j];
  for (int i = 0; i < ZMAX; i++)
    graph[i][i] = 1;
  int N = 2000;
  while (N < ZMAX)
  {

  //printMatrix(graph, 6, 9);
    computeDYN0(graph, N);
    computeDYN1(graph, N);
    computeDYN2(graph, N);
    N += 10;
  }
  deallocateMatrix(graph, ZMAX);
  return 0;
}
