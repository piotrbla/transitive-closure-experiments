//kompilacja np.: gcc -DMYTHREADS=32 -o bstexe_32 -O2 bst.cpp -fopenmp -lm
//kompilacja np.: gcc -DMYTHREADS=16 -o bstexe_16 -O2 bst.cpp -fopenmp -lm

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

void computeOBSTNSeq(int* p, int n)
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

  for (i=n-1 ; i>=1 ; i--)
  for (j = i+1 ; j<=n ; j+= 1)
    for (k = i+1 ; k<j; k += 1) {
      c[i][j] = min(c[i][j], w[i][j]+c[i][k]+c[k][j]);
    }

  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(c, n, 1);
  deallocateMatrix(w, n + 1);
  deallocateMatrix(c, n + 1);
}

void computeOBSTISeq(int* p, int n)
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

#define SI0(a, i, j, k)       c[i][j] = min(c[i][j], w[i][j]+c[i][k]+c[k][j]);

 for (int c0 = 2; c0 < n; c0 += 1)
    for (int c1 = 1; c1 <= n - c0; c1 += 1)
      for (int c2 = c0 + c1; c2 <= min(n, 2 * c0 + c1 - 2); c2 += 1) {
        if (2 * c0 + c1 >= c2 + 3)
          SI0(c0, c1, c2, -c0 + c2 + 1);
        SI0(c0, c1, c2, c0 + c1 - 1);
      }

  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
//  write_results(n, execution_time);
  write_results(n, execution_time);
  printMatrix(c, n, 3);
  deallocateMatrix(w, n + 1);
  deallocateMatrix(c, n + 1);
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

for (int c0 = 1; c0 <= floord(n - 1, 12) + 1; c0 += 1) {
    #pragma omp parallel for num_threads(MYTHREADS)
    for (int c1 = max(max(-((n + 11) / 12), -c0 - (n + 62) / 64 + 1), -c0 - (3 * c0 + 12) / 13 + 1); c1 <= -c0; c1 += 1) {
      for (int c4 = max(max(-n + 2, 64 * c0 + 64 * c1 - 63), 12 * c1 + 2); c4 <= min(-1, 64 * c0 + 64 * c1); c4 += 1) {
        for (int c5 = max(-12 * c1 - 11, -c4 + 2); c5 <= min(n, -12 * c1); c5 += 1) {
          for (int c6 = -c4 + 1; c6 < c5; c6 += 1) {
            c[-c4][c5] = ((c[-c4][c5] < ((w[-c4][c5] + c[-c4][c6]) + c[c6][c5])) ? c[-c4][c5] : ((w[-c4][c5] + c[-c4][c6]) + c[c6][c5]));
          }
        }
      }
    }
  }

  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(c, n, 1);
  deallocateMatrix(w, n + 1);
  deallocateMatrix(c, n + 1);
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

for (int c0 = 1; c0 <= floord(n - 1, 12) + 1; c0 += 1) {
    #pragma omp parallel for num_threads(MYTHREADS)
    for (int c1 = max(max(-((n + 11) / 12), -c0 - (n + 62) / 64 + 1), -c0 - (3 * c0 + 12) / 13 + 1); c1 <= -c0; c1 += 1) {
      for (int c4 = max(max(-n + 2, 64 * c0 + 64 * c1 - 63), 12 * c1 + 2); c4 <= min(-1, 64 * c0 + 64 * c1); c4 += 1) {
        for (int c5 = max(-12 * c1 - 11, -c4 + 2); c5 <= min(n, -12 * c1); c5 += 1) {
          for (int c6 = -c4 + 1; c6 < c5; c6 += 1) {
            c[-c4][c5] = ((c[-c4][c5] < ((w[-c4][c5] + c[-c4][c6]) + c[c6][c5])) ? c[-c4][c5] : ((w[-c4][c5] + c[-c4][c6]) + c[c6][c5]));
          }
        }
      }
    }
  }

  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(c, n, 1);
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

  #define SI0(a, i, j, k)       c[i][j] = min(c[i][j], w[i][j]+c[i][k]+c[k][j]);

  for (int c0 = 2; c0 <= floord(19 * n - 22, 192) + 2; c0 += 1) {
    #pragma omp parallel for num_threads(MYTHREADS)
    for (int c1 = max(max(-4 * c0 + 5, -((n + 11) / 12)), -c0 - (n + 14) / 16 + 2); c1 <= min(-c0 + (n - 2) / 64 + 1, -c0 + (3 * c0 - 4) / 19 + 1); c1 += 1) {
      for (int c2 = max(-((4 * c0 + c1 + 15) / 20), -((n + 16 * c0 + 16 * c1 + 63) / 80)); c2 <= min(min(-1, -c0 - c1), -((4 * c0 + c1 + 29) / 36)); c2 += 1) {
        for (int c5 = max(max(2, -64 * c2 - 63), 8 * c0 + 2 * c1 + 8 * c2 - 12); c5 <= min(min(min(min(n - 1, -12 * c1 - 1), -64 * c2), 16 * c0 + 4 * c1 + 16 * c2), n + 16 * c0 + 16 * c1 + 16 * c2); c5 += 1) {
          for (int c6 = max(max(1, -16 * c0 - 16 * c1 - 16 * c2), -12 * c1 - 2 * c5 - 9); c6 <= min(min(-16 * c0 - 16 * c1 - 16 * c2 + 15, -12 * c1 - c5), n - c5); c6 += 1) {
            for (int c7 = max(-12 * c1 - 11, c5 + c6); c7 <= min(min(n, -12 * c1), 2 * c5 + c6 - 2); c7 += 1) {
              if (2 * c5 + c6 >= c7 + 3) {
                c[c6][c7] = ((c[c6][c7] < ((w[c6][c7] + c[c6][-c5 + c7 + 1]) + c[-c5 + c7 + 1][c7])) ? c[c6][c7] : ((w[c6][c7] + c[c6][-c5 + c7 + 1]) + c[-c5 + c7 + 1][c7]));
              }
              c[c6][c7] = ((c[c6][c7] < ((w[c6][c7] + c[c6][c5 + c6 - 1]) + c[c5 + c6 - 1][c7])) ? c[c6][c7] : ((w[c6][c7] + c[c6][c5 + c6 - 1]) + c[c5 + c6 - 1][c7]));
            }
          }
        }
      }
    }
  }


  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
//  write_results(n, execution_time);
  write_results_full(MYTHREADS, execution_time, '\n');
  printMatrix(c, n, 3);
  deallocateMatrix(w, n + 1);
  deallocateMatrix(c, n + 1);
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
  const int ZMAX = 2300; //ZMAXSIZE;
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
  int N = 2100;//ZMAX;
  while (N <= ZMAX)
  {
    computeOBSTNSeq(seq, N);
    computeOBSTISeq(seq, N);
    computeOBSTN(seq, N);
    computeOBSTI(seq, N);
    N += 100;
  }
  free(seq);
  deallocateMatrix(graph, ZMAX);
  return 0;
}
                              