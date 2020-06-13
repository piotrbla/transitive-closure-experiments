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
#define match(b1, b2) (((b1)+(b2)) == 3 ? 1 : 0)
#define sigma(i, j) (match(seq[i], seq[j]))
#define max_score(s1, s2) ((s1 >= s2) ? s1 : s2)
#define max_sc(s1, s2, s3) ((s1 >= s2) ? ((s1 >= s3) ? s1 : s3) : ((s2 >= s3) ? s2 : s3))



void printMatrix(int**, int, int);
int ** getFullCopy(int ** table, int N);
int** allocateMatrix(int);
void deallocateMatrix(int**, int);

void write_results_full(int , double , char );
void write_results(int , double );

void computeDYN0Imperfect(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();
  for (int i = n - 1; i >= 1; i--) {
    for (int j = i + 1; j < n; j++) {
      for (int k = i; k < j - 1; k++) {
        S[i][j] = max_score(S[i][k - 1] + S[k + 1][j - 1], S[i][j]); // s1
      }
      S[i][j] = max_score(S[i][j - 1], S[i][j] + sigma(i, j)); // s2
    }
  }
  double execution_time = omp_get_wtime() - start;

  printf("IMP: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(S, n, 0);
  deallocateMatrix(S, n);
}

void computeDYN0Perfect(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();
  //Listing 1.2: Perfectly nested Nussinov loops
  for (int i = n - 1; i >= 1; i--) {
    for (int j = i + 1; j < n; j++) {
      for (int k = i; k < j - 1; k++) {
        S[i][j] = max_sc(S[i][j], S[i][j-1], S[i][k-1] + S[k + 1][j - 1] + sigma(i, j));
      }
    }
  }
  double execution_time = omp_get_wtime() - start;

  printf("PER: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(S, n, 1);
  deallocateMatrix(S, n);
}

void computeDYN0TestbenchImplementation(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();
  for (int i = n - 1; i >= 0; i--) {
    for (int j = i + 1; j < n; j++) {

      if (j - 1 >= 0)
        S[i][j] = max_score(S[i][j], S[i][j - 1]);
      if (i + 1 < n)
        S[i][j] = max_score(S[i][j], S[i + 1][j]);

      if (j - 1 >= 0 && i + 1 < n) {
        if (i < j - 1)
          S[i][j] = max_score(S[i][j], S[i + 1][j - 1] + sigma(i, j));
        else
          S[i][j] = max_score(S[i][j], S[i + 1][j - 1]);
      }

      for (int k = i + 1; k < j; k++) {
        S[i][j] = max_score(S[i][j], S[i][k] + S[k + 1][j]);
      }
    }
  }

  double execution_time = omp_get_wtime() - start;
  printf("TSB: %lf\n", execution_time);
  write_results_full(n, execution_time, '\n');
  printMatrix(S, n, 2);
  deallocateMatrix(S, n);
}


void printMatrix(int** matrix, int N, int fileno) {
  char filename[10];
  sprintf(filename, "nontiled%d", fileno);
  FILE* f = fopen(filename, "wt");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++)
      fprintf(f, "%d ", matrix[i][j]);
    fprintf(f, "\n");
  }
  fclose(f);
}

int **getFullCopy(int ** table, int N)
{
  int **S = allocateMatrix(N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      S[i][j] = table[i][j];
  return S;
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

int getValue(const char c)
{
  if(c=='A')
    return 0;
  if(c=='G')
    return 1;
  if(c=='U')
    return 2;
  if(c=='C')
    return 3;
  return 4;
}

int main(void) {
  const int ZMAX = 16;
  int** graph = allocateMatrix(ZMAX);
  int* seq = allocateVector(ZMAX);
  for (int i = 0; i < ZMAX; i++)
    for (int j = 0; j < ZMAX; j++)
      graph[i][j] = 0;
  for (int i = 0; i < ZMAX; i++)
    graph[i][i] = 0;
  //
  const char *seqTest = "GCGUCCACGGCUAGCU";
  //for (int i=0 ; i<ZMAX ; i++)
  //  seq[i] = rand()%4;
  for (int i=0 ; i<ZMAX ; i++)
    seq[i] = getValue(seqTest[i]);
  int N = ZMAX - 10;
  //while (N < ZMAX)
  //{
    N += 10;
    computeDYN0Imperfect(graph, N, seq);
    computeDYN0Perfect(graph, N, seq);
    computeDYN0TestbenchImplementation(graph, N, seq);
    //N += 10;
  //}
  deallocateMatrix(graph, ZMAX);
  free(seq);
  return 0;
}
