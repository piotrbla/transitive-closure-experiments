#include <math.h>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))

bool is_files_equal(const char * , const char * ) ;
void print_files_equal(const char * , const char * );
void write_results_full(int, double, char);
void write_results(int, double );
void write_results_last(int, double );
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

void LMKernel6_02_Modification_Bad(int loop, int n, int *input_w, int **input_b)
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

void LMKernel6_04_Modification_DAPT(int loop, int n, int *input_w, int **input_b)
{
  int** b = get_full_copy(input_b, n);
  int* w = get_vector_copy(input_w, n);

  double start = omp_get_wtime();


for (int i0 = 1; i0 <= loop; i0 += 1) {
  for (int w0 = 2; w0 <= floord(n - 2, 8) + 2; w0 += 1) {
#pragma omp parallel for
    for (int h0 = max(-w0 + 1, -((n + 14) / 16)); h0 <= -w0 + w0 / 2; h0 += 1) {
      for (int i1 = 16 * w0 + 16 * h0 - 15; i1 <= min(n - 1, 16 * w0 + 16 * h0); i1 += 1) {
        for (int i2 = max(-16 * h0 - 15, i1); i2 <= min(n - 1, -16 * h0); i2 += 1) {
          w[i2] += (b[-i1 + i2][i2] * w[i1 - 1]);
        }
      }
    }
  }
}

  double execution_time = omp_get_wtime() - start;

  printf("MOD_05: %lf\n", execution_time);
  write_results(n, execution_time);
  print_vector(w, n);
  deallocate_matrix(b, n);
  free(w);
  return;
}


void LMKernel6_05_Modification_DAPT2(int loop, int n, int *input_w, int **input_b)
{
  int** b = get_full_copy(input_b, n);
  int* w = get_vector_copy(input_w, n);

  double start = omp_get_wtime();

  for (int i0 = 1; i0 <= loop; i0 += 1){
    for (int w0 = 0; w0 <= floord(n - 1, 8); w0 += 1){
#pragma omp parallel for
      for (int h0 = max(0, w0 - (n + 15) / 16 + 1); h0 <= w0 / 2; h0 += 1){
        for (int t0 = max(1, 8 * w0); t0 <= min(min(n - 1, 8 * w0 + 15), 8 * h0 + n / 2 + 7); t0 += 1){
          for (int i1 = max(max(max(1, 16 * h0), -16 * w0 + 16 * h0 + 2 * t0 - 15), -n + 2 * t0 + 1); i1 <= min(min(16 * h0 + 15, t0), -16 * w0 + 16 * h0 + 2 * t0 + 1); i1 += 1){
            for (int i2 = max(16 * w0 - 16 * h0, 2 * t0 - i1); i2 <= min(min(n - 1, 16 * w0 - 16 * h0 + 15), 2 * t0 - i1 + 1); i2 += 1){
              w[i2] += (b[-i1 + i2][i2] * w[i1 - 1]);
            }
          }
        }
      }
    }
  }

  double execution_time = omp_get_wtime() - start;

  printf("MOD_05: %lf\n", execution_time);
  write_results(n, execution_time);
  print_vector(w, n);
  deallocate_matrix(b, n);
  free(w);
  return;
}

void LMKernel6_06_Modification_Pluto(int loop, int n, int *input_w, int **input_b)
{
  int** b = get_full_copy(input_b, n);
  int* w = get_vector_copy(input_w, n);

  double start = omp_get_wtime();

    int t1, t2, t3, t4, t5;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if ((loop >= 1) && (n >= 2)) {
  for (t1=1;t1<=loop;t1++) {
    for (t2=0;t2<=floord(n-1,16);t2++) {
      lbp=ceild(t2,2);
      ubp=min(floord(n-1,32),t2);
#pragma omp parallel for private(lbv,ubv,t4,t5)
      for (t3=lbp;t3<=ubp;t3++) {
        for (t4=max(1,32*t3);t4<=min(n-1,32*t3+31);t4++) {
          for (t5=max(1,32*t2-32*t3);t5<=min(t4,32*t2-32*t3+31);t5++) {
            w[t4] += b[-t5 + t4][t4] * w[t5 - 1];;
          }
        }
      }
    }
  }
}

  double execution_time = omp_get_wtime() - start;

  printf("MOD_06: %lf\n", execution_time);
  write_results_last(n, execution_time);
  print_vector(w, n);
  deallocate_matrix(b, n);
  free(w);
  return;
}


#define PERFORMANCE_TEST 1

void make_work_one_size(const int ZMAX, const int LMAX)
{
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

  print_files_equal("resVec_1", "resVec_2");
  LMKernel6_01_Base(LMAX, ZMAX, result_w, input_b);
  LMKernel6_02_Modification_Bad(LMAX, ZMAX, result_w, input_b);
  LMKernel6_03_Modification(LMAX, ZMAX, result_w, input_b);
  LMKernel6_04_Modification_DAPT(LMAX, ZMAX, result_w, input_b);
  LMKernel6_05_Modification_DAPT2(LMAX, ZMAX, result_w, input_b);
  LMKernel6_06_Modification_Pluto(LMAX, ZMAX, result_w, input_b);
  deallocate_matrix(input_b, ZMAX);
  free(result_w);
}


int main()
{
  //TODO: change int to double
#if PERFORMANCE_TEST==1
    const int ZMAX = 5500;
    const int LMAX = 10;
    for (int z = 1000 ; z <ZMAX ; z+=500)
      make_work_one_size(z, LMAX);
#else
  const int ZMAX = 16;
  const int LMAX = 10;
  make_work_one_size(ZMAX, LMAX);
#endif 

  print_files_equal("resVec_1", "resVec_2");
  print_files_equal("resVec_1", "resVec_3");
  print_files_equal("resVec_1", "resVec_4");
  print_files_equal("resVec_1", "resVec_5");
  print_files_equal("resVec_1", "resVec_6");
  
}

bool is_files_equal(const char * filename_template, const char * filename_compared) {
  std::ifstream f1(filename_template, std::ifstream::binary | std::ifstream::ate);
  std::ifstream f2(filename_compared, std::ifstream::binary | std::ifstream::ate);
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

void print_files_equal(const char* filename_template, const char* filename_compared)
{
  if (is_files_equal(filename_template, filename_compared))
    std::cout << filename_compared << ": OK\n";
  else
    std::cout << filename_compared << ": ERROR\n";
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

void write_results_last(int n, double execution_time)
{
  write_results_full(n, execution_time, '\n');
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
