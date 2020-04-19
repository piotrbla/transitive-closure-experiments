#pragma warning(disable : 4996)
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define _CRT_SECURE_NO_WARNINGS
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S0(a, i, j, k) c[i][j] = c[i][k] + c[k][j]

void printMatrix(int**, int, int);
int** allocateMatrix(int);
void deallocateMatrix(int**, int);

void write_results_full(int , double , char );
void write_results(int , double );

void computeSEQ0(int* p, int n)
{
  int** w = allocateMatrix(n + 1);
  int i, j, k, m, h, optimal_w;
  for (i = 0; i <= n; i++)
  {
    for (j = 0; j <= n; j++)
        w[i][j] = 0;
  }

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
        w[i][j] = 0;
    if (i<n-1)
    {
      w[i][i+1] = p[i];
    }
  }
  for (i = 0; i < n; i++)
     for (j = i+1; j < n; j++)
        w[i][j] = 99999;
  double start = omp_get_wtime();
#pragma scop
  for (j = 1; j < n; j++)
    for (i=1; i < n-j; i++){
        for (k = i; k <= i+j; k++){
           w[i][i+j] = min(w[i][k-1] + w[k+1][i+j], w[i][i+j]);
        }
        for (k = i; k <= i+j; k++)
                w[i][i+j] += p[k];
     }
#pragma endscop
  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(w, n, 0);
  deallocateMatrix(w, n + 1);
}

void computeOBSTN(int* p, int n)
{
  int** w = allocateMatrix(n + 1);
  int** c = allocateMatrix(n + 1);
  int i, j, k, m, h, optimal_w;
  for (i = 0; i <= n; i++)
  {
    for (j = 0; j <= n; j++)
    {
        w[i][j] = 0;
        c[i][j] = 0;
    }
  }
  for (i = 0; i < n; i++)
     for (j = i+1; j < n; j++)
        w[i][j] = 99999;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      c[i][j] = p[i];
      w[i][j] = p[i];
      if (i<n-1)
      {
        c[i][i+1] = p[i];
        w[i][i+1] = p[i];
      }
    }
  }
  double start = omp_get_wtime();

  #pragma scop
  for (i=n-1 ; i>=1 ; i--)
  for (j = i+1 ; j<=n ; j+= 1)
    for (k = i+1 ; k<j; k += 1) {
      c[i][j] = min(c[i][j], w[i][j]+c[i][k]+c[k][j]);
    }
#pragma endscop
  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(c, n, 1);
  deallocateMatrix(w, n + 1);
  deallocateMatrix(c, n + 1);
}

void computeOBSTK(int* p, int n)
{
  int** w = allocateMatrix(n + 1);
  int** c = allocateMatrix(n + 1);
  int i, j, k, m, h, optimal_w;
  for (i = 0; i <= n; i++)
  {
    for (j = 0; j <= n; j++)
    {
        w[i][j] = 0;
        c[i][j] = 0;
    }
  }
  for (i = 0; i < n; i++)
     for (j = i+1; j < n; j++)
        w[i][j] = 99999;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      c[i][j] = p[i];
      w[i][j] = p[i];
      if (i<n-1)
      {
        c[i][i+1] = p[i];
        w[i][i+1] = p[i];
      }
    }
  }
  double start = omp_get_wtime();

  #pragma scop
  for (i=n-1; i >=1 ; i--)
  for (j = i+1; j < n ; j += 1)
    for (k = i; k <=j; k += 1) {
      c[i][j]=min(c[i][j], w[i][j]+c[i][k-1]+c[k+1][j]);
    }
#pragma endscop
  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
//  write_results(n, execution_time);
  write_results(n, execution_time);
  printMatrix(c, n, 2);
  deallocateMatrix(w, n + 1);
  deallocateMatrix(c, n + 1);
}

void computeOBSTI(int* p, int n)
{
  int** w = allocateMatrix(n + 1);
  int** c = allocateMatrix(n + 1);
  int i, j, k, m, h, optimal_w;
  for (i = 0; i <= n; i++)
  {
    for (j = 0; j <= n; j++)
    {
        w[i][j] = 0;
        c[i][j] = 0;
    }
  }
  for (i = 0; i < n; i++)
     for (j = i+1; j < n; j++)
        w[i][j] = 99999;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      c[i][j] = p[i];
      w[i][j] = p[i];
      if (i<n-1)
      {
        c[i][i+1] = p[i];
        w[i][i+1] = p[i];
      }
    }
  }
  double start = omp_get_wtime();

  #pragma scop
  #define SI0(a, i, j, k)       c[i][j] = min(c[i][j], w[i][j]+c[i][k]+c[k][j]);

 for (int c0 = 2; c0 < n; c0 += 1)
    for (int c1 = 1; c1 <= n - c0; c1 += 1)
      for (int c2 = c0 + c1; c2 <= min(n, 2 * c0 + c1 - 2); c2 += 1) {
        if (2 * c0 + c1 >= c2 + 3)
          SI0(c0, c1, c2, -c0 + c2 + 1);
        SI0(c0, c1, c2, c0 + c1 - 1);
      }
#pragma endscop
  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
//  write_results(n, execution_time);
  write_results_full(n, execution_time, '\n');
  printMatrix(c, n, 3);
  deallocateMatrix(w, n + 1);
  deallocateMatrix(c, n + 1);
}

void computeDYN2(int** matrix, int n) {
  int** c = allocateMatrix(n + 1);
  int** d = allocateMatrix(n + 1);
  int i, j;
  for (i = 0; i <= n; i++)
    for (j = 0; j <= n; j++)
      c[i][j] = 0;

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
  printf("tiled if: %lf\n", execution_time);
  write_results_full(n, execution_time, '\n');
  printMatrix(c, n, 4);
  deallocateMatrix(c, n + 1);
  deallocateMatrix(d, n + 1);
}


void printMatrix(int** matrix, int N, int fileno) {
  char filename[10];
  sprintf(filename, "nontiled%d", fileno);
  FILE* f = fopen(filename, "wt");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++)
    {
      fprintf(f, "%d ", matrix[i][j]);
      //if (fileno==2)
        //printf("%4d ", matrix[i][j]);
    }
    fprintf(f, "\n");
    //if (fileno==2)
      //printf("\n");
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

int* allocateVector(int N) {
  int* t = (int*)malloc(sizeof(int) * N);
  return t;
}

void deallocateMatrix(int **t, int N) {
  for (int i = 0; i < N; i++) {
    free(t[i]);
  }
  free(t);
}

void write_results_full(int n, double execution_time, char end_char)
{
  FILE* f = fopen("results.txt", "at");
  fprintf(f, "%d:%lf%c", n, execution_time, end_char);
  fclose(f);
}
void write_results(int n, double execution_time)
{
  write_results_full(n, execution_time, ';');
}


int main(void) {//vector
  const int ZMAX = 2500; //ZMAXSIZE;
  int* seq = allocateVector(ZMAX);
  for (int i = 0; i < ZMAX; i++)
    seq[i] = i;
  int** graph = allocateMatrix(ZMAX);
  int g[4][4] = { {1, 1, 0, 1}, {0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1} };
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      graph[i][j] = g[i][j];
  for (int i = 0; i < ZMAX; i++)
    graph[i][i] = 1;
  int N = 500;//ZMAX;
  while (N <= ZMAX)
  {
    //computeSEQ0(seq, N);
    computeOBSTN(seq, N);
    computeOBSTK(seq, N);
    computeOBSTI(seq, N);
    //computeDYN2(graph, N);
    N += 100;
  }
  free(seq);
  deallocateMatrix(graph, ZMAX);
  return 0;
}
