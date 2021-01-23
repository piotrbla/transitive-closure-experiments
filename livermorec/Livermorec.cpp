#pragma warning(disable : 4996)
#include <fstream>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#define _CRT_SECURE_NO_WARNINGS

bool is_files_equal(const char * , const char * ) ;
void write_results_full(int, double, char);
void write_results(int, double );
int** allocate_matrix(int) ;

int **get_full_copy(int ** , int);
int* allocate_vector(int);
int *get_vector_copy(int *, int);

void deallocate_matrix(int **, int ) ;
void print_vector(int* , int );
void print_matrix(int**, int);

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
  int** b = get_full_copy(input_b, n);
  int* w = get_vector_copy(input_w, n);

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
  print_vector(w, n);
  print_matrix(b, n);
  deallocate_matrix(b, n);
  free(w);
  return;
}

void LMKernel6_02_Modification(int loop, int n, int *input_w, int **input_b)
{
  int** b = get_full_copy(input_b, n);
  int* w = get_vector_copy(input_w, n);

  double start = omp_get_wtime();

for (int c0 = 1; c0 < n; c0 += 1)
  for (int c1 = 1; c1 <= loop; c1 += 1)
    for (int c2 = c0; c2 < n; c2 += 1)
        w[c2] += b[-c0 + c2][c2] * w[(c2 - (-c0 + c2)) - 1];
  double execution_time = omp_get_wtime() - start;

  printf("MOD_02: %lf\n", execution_time);
  write_results(n, execution_time);
  print_vector(w, n);
  deallocate_matrix(b, n);
  free(w);
  return;
}

void LMKernel6_03_Modification(int loop, int n, int *input_w, int **input_b)
{
  int** b = get_full_copy(input_b, n);
  int* w = get_vector_copy(input_w, n);

  double start = omp_get_wtime();
for (int c0 = 1; c0 <= loop; c0 += 1)
  for (int c1 = 1; c1 < n; c1 += 1)
    for (int c2 = c1; c2 < n; c2 += 1)
       w[c2] += b[-c1 + c2][c2] * w[(c2 - (-c1 + c2)) - 1];

  double execution_time = omp_get_wtime() - start;

  printf("MOD_03: %lf\n", execution_time);
  write_results(n, execution_time);
  print_vector(w, n);
  deallocate_matrix(b, n);
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

  int** input_b = allocate_matrix(ZMAX);
  int* result_w = allocate_vector(ZMAX);
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
  LMKernel6_03_Modification(LMAX, ZMAX, result_w, input_b);
  deallocate_matrix(input_b, ZMAX);
  free(result_w);

  if (is_files_equal("resVec_1", "resVec_2"))
    std::cout << "OK\n";
  else
    std::cout << "ERROR\n";
  if (is_files_equal("resVec_1", "resVec_3"))
    std::cout << "OK\n";
  else
    std::cout << "ERROR\n";


}

bool is_files_equal(const char * filename1, const char * filename2) {
  std::ifstream f1(filename1, std::ifstream::binary | std::ifstream::ate);
  std::ifstream f2(filename2, std::ifstream::binary | std::ifstream::ate);
  if (f1.fail() || f2.fail()) {
    return false;
  }
  if (f1.tellg() != f2.tellg()) {
    return false;
  }
  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator <char >(f1.rdbuf()),
    std::istreambuf_iterator <char >(),
    std::istreambuf_iterator <char >(f2.rdbuf()));
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

int** allocate_matrix(int N) {
  int** t = (int**)malloc(sizeof(int*) * N);
  for (int i = 0; i < N; i++) {
    t[i] = (int*)malloc(sizeof(int) * N);
  }
  return t;
}

int **get_full_copy(int ** table, int N)
{
  int **S = allocate_matrix(N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      S[i][j] = table[i][j];
  return S;
}

int* allocate_vector(int N) {
  int* t = (int*)malloc(sizeof(int) * N);
  return t;
}

int *get_vector_copy(int *table, int N)
{
  int *S = allocate_vector(N);
  for (int i = 0; i < N; i++)
      S[i] = table[i];
  return S;
}


void deallocate_matrix(int **t, int N) {
  for (int i = 0; i < N; i++) {
    free(t[i]);
  }
  free(t);
}

void print_vector(int* vector, int N) {
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

void print_matrix(int** matrix, int N) {
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
