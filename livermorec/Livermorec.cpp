#pragma warning(disable : 4996)
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define _CRT_SECURE_NO_WARNINGS

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

int** allocateMatrix(int N) {
  int** t = (int**)malloc(sizeof(int*) * N);
  for (int i = 0; i < N; i++) {
    t[i] = (int*)malloc(sizeof(int) * N);
  }
  return t;
}

int **getFullCopy(int ** table, int N)
{
  int **S = allocateMatrix(N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      S[i][j] = table[i][j];
  return S;
}

int* allocateVector(int N) {
  int* t = (int*)malloc(sizeof(int) * N);
  return t;
}

int *getVectorCopy(int *table, int N)
{
  int *S = allocateVector(N);
  for (int i = 0; i < N; i++)
      S[i] = table[i];
  return S;
}


void deallocateMatrix(int **t, int N) {
  for (int i = 0; i < N; i++) {
    free(t[i]);
  }
  free(t);
}

void printVector(int* vector, int N) {
  static int fileno=1;
  char filename[12];
  sprintf(filename, "resVec_%d", fileno);
  FILE* f = fopen(filename, "wt");
  for (int i = 0; i < N; i++) {
    fprintf(f, "%d ", vector[i]);
    fprintf(f, "\n");
  }
  fclose(f);
  fileno++;
}

void printMatrix(int** matrix, int N) {
  static int fileno=1;
  char filename[10];
  sprintf(filename, "resMat_%d", fileno);
  FILE* f = fopen(filename, "wt");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++)
      fprintf(f, "%d ", matrix[i][j]);
    fprintf(f, "\n");
  }
  fclose(f);
  fileno++;
}


/*
  *******************************************************************
  *   Kernel 6 -- general linear recurrence equations
  *******************************************************************
  *    DO  6  L= 1,Loop
  *    DO  6  i= 2,n
  *    DO  6  k= 1,i-1
  *        W(i)= W(i)  + B(i,k) * W(i-k)
  *  6 CONTINUE
*/

void LMKernel6_01_Base(int loop, int n, int *input_w, int **input_b)
{
  int** b = getFullCopy(input_b, n);
  int* w = getVectorCopy(input_w, n);

  double start = omp_get_wtime();

  int l, i, k;

  for (l = 1; l <= loop; l++) {
    for (i = 1; i < n; i++) {
      for (k = 0; k < i; k++) {
        w[i] += b[k][i] * w[(i - k) - 1];
      }
    }
  }
  double execution_time = omp_get_wtime() - start;

  printf("BAS_01: %lf\n", execution_time);
  write_results(n, execution_time);
  printVector(w, n);
  printMatrix(b, n);
  deallocateMatrix(b, n);
  free(w);
  return;
}

void LMKernel6_02_Modification(int loop, int n, int *input_w, int **input_b)
{
  int** b = getFullCopy(input_b, n);
  int* w = getVectorCopy(input_w, n);

  double start = omp_get_wtime();

  int l, i, k;

  for (l = 1; l <= loop; l++) {
    for (i = 1; i < n; i++) {
      for (k = 0; k < i; k++) {
        w[i] += b[k][i] * w[(i - k) - 1];
      }
    }
  }
  double execution_time = omp_get_wtime() - start;

  printf("MOD_02: %lf\n", execution_time);
  write_results(n, execution_time);
  printVector(w, n);
  deallocateMatrix(b, n);
  free(w);
  return;
}

#define PERFORMANCE_TEST 0

int main()
{
  //TODO: change int to double
#if PERFORMANCE_TEST==1
    const int ZMAX = 2600;
    const int LMAX = 10;
#else
  const int ZMAX = 16;
  const int LMAX = 10;
#endif

  int** input_b = allocateMatrix(ZMAX);
  int* result_w = allocateVector(ZMAX);
  for (int i = 0; i < ZMAX; i++)
    for (int j = 0; j < ZMAX; j++)
      input_b[i][j] = 0;
  for (int i = 0; i < ZMAX; i++)
    result_w[i] = 1;
  const char* seqTest = "1234432432123412";
#if PERFORMANCE_TEST==1
  for (int i = 0; i < ZMAX; i++)
    for (int j = 0; j < ZMAX; j++)
      input_b[i][j] = rand()%4+1;
#else
  for (int i = 0; i < ZMAX; i++)
    for (int j = 0; j < ZMAX; j++)
      input_b[i][j] = seqTest[i]-'0';
      
#endif


  LMKernel6_01_Base(LMAX, ZMAX, result_w, input_b);
  LMKernel6_02_Modification(LMAX, ZMAX, result_w, input_b);
  deallocateMatrix(input_b, ZMAX);
  free(result_w);

}

