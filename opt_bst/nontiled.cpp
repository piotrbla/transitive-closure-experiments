#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S0(a, i, j, k) c[i][j] = c[i][k] + c[k][j]
#define match(b1, b2) (((b1)+(b2)) == 3 ? 1 : 0)
#define max_score(s1, s2) ((s1 >= s2) ? s1 : s2)

void printMatrix(int**, int, int);
int** allocateMatrix(int);
void deallocateMatrix(int**, int);

void write_results(int , double , char );
void write_results(int , double );

void computeSEQ0(int* p, int n)
{
  int** w = allocateMatrix(n + 1);
  int i, j, k, m, h, optimal_w;
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



  for (j = 1; j < n; j++)
    for (i=0; i < n-j; i++){    
        for (k = i; k <= i+j; k++){
           optimal_w = w[i][k-1] + w[k+1][i+j]; 
           if (optimal_w < w[i][i+j]){   
              w[i][i+j] = optimal_w;	  
           }
        }
        for (k = i; k <= i+j; w[i][i+j] += p[k++])
          ;
     }
  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(w, n, 0);
  deallocateMatrix(w, n + 1);
}


void computeDYN0(int** p, int n) {
  int** w = allocateMatrix(n + 1);
  int** d = allocateMatrix(n + 1);
  int i, j, k, optimal_w;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      w[i][j] = p[i][j];
  double start = omp_get_wtime();

  for (i = 0; i < n; i++) {//, i < k < j
    for (j = 0; j < i; j++) {
      optimal_w = 9999;
      for (k = 0; k < j; k++) {
        optimal_w = min(w[i][k], optimal_w);
      }
      w[i][j] = p[i][j] + optimal_w;
    }
  }

  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(w, n, 0);
  deallocateMatrix(w, n + 1);
  deallocateMatrix(d, n + 1);
}

void computeDYN1(int** table, int _PB_N, int* seq) {
  int** c = allocateMatrix(_PB_N + 1);
  int** d = allocateMatrix(_PB_N + 1);
  int i, j, k;
  for (i = 0; i < _PB_N; i++)
    for (j = 0; j < _PB_N; j++)
      c[i][j] = table[i][j];
  double start = omp_get_wtime();

#pragma omp parallel for
  for (int i = _PB_N - 1; i >= 0; i--) {
    for (int j = i + 1; j < _PB_N; j++) {

      if (j - 1 >= 0)
        c[i][j] = max_score(c[i][j], c[i][j - 1]);
      if (i + 1 < _PB_N)
        c[i][j] = max_score(c[i][j], c[i + 1][j]);

      if (j - 1 >= 0 && i + 1 < _PB_N) {
        if (i < j - 1)
          c[i][j] = max_score(c[i][j], c[i + 1][j - 1] + match(seq[i], seq[j]));
        else
          c[i][j] = max_score(c[i][j], c[i + 1][j - 1]);
      }

      for (k = i + 1; k < j; k++) {
        c[i][j] = max_score(c[i][j], c[i][k] + c[k + 1][j]);
      }
    }
  }

  double execution_time = omp_get_wtime() - start;
  printf("parallel: %lf\n", execution_time);
  write_results(_PB_N, execution_time);
  printMatrix(c, _PB_N, 1);
  deallocateMatrix(c, _PB_N + 1);
  deallocateMatrix(d, _PB_N + 1);
}

void computeDYN2(int** table, int _PB_N, int* seq) {
  int** c = allocateMatrix(_PB_N + 1);
  int** d = allocateMatrix(_PB_N + 1);
  int i, j, k;
  for (i = 0; i < _PB_N; i++)
    for (j = 0; j < _PB_N; j++)
      c[i][j] = table[i][j];
  double start = omp_get_wtime();

  for (int i = _PB_N - 1; i >= 0; i--) {
    for (int j = i + 1; j < _PB_N; j++) {

      if (j - 1 >= 0)
        c[i][j] = max_score(c[i][j], c[i][j - 1]);
      if (i + 1 < _PB_N)
        c[i][j] = max_score(c[i][j], c[i + 1][j]);

      if (j - 1 >= 0 && i + 1 < _PB_N) {
        if (i < j - 1)
          c[i][j] = max_score(c[i][j], c[i + 1][j - 1] + match(seq[i], seq[j]));
        else
          c[i][j] = max_score(c[i][j], c[i + 1][j - 1]);
      }
      
      for (k = i + 1; k < j; k++) {
        c[i][j] = max_score(c[i][j], c[i][k] + c[k + 1][j]);
      }
    }
  }

  double execution_time = omp_get_wtime() - start;
  printf("tiles: %lf\n", execution_time);
  //write_results(_PB_N, execution_time, '\n');
  printMatrix(c, _PB_N, 2);
  deallocateMatrix(c, _PB_N + 1);
  deallocateMatrix(d, _PB_N + 1);
}

void printMatrix(int** matrix, int N, int fileno) {
  char filename[10];
  sprintf_s(filename, "nontiled%d", fileno);
  FILE* f;
  fopen_s(&f, filename, "wt");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++)
    {
      fprintf(f, "%d ", matrix[i][j]);
      printf("%d ", matrix[i][j]);
    }
    fprintf(f, "\n");
    printf("\n");
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

void write_results(int n, double execution_time, char end_char)
{
  FILE* f;
  fopen_s(&f, "results.txt", "at");
  fprintf(f, "%d:%lf%c", n, execution_time, end_char);
  fclose(f);
}
void write_results(int n, double execution_time)
{
  write_results(n, execution_time, ';');
}


int main(void) {//vector
  const int ZMAX = 7;
  int* seq = allocateVector(ZMAX);
  for (int i = 0; i < ZMAX; i++)
    seq[i] = i;
  int N = ZMAX;
  while (N <= ZMAX)
  {

    computeSEQ0(seq, N);
    N += 10;
  }
  free(seq);
  return 0;
}

//int mainMatrix(void) {
//  const int ZMAX = 510;
//  int** graph = allocateMatrix(ZMAX);
//  int* seq = allocateVector(ZMAX);
//  int g[4][4] = { {1, 1, 0, 1}, {0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1} };
//  for (int i = 0; i < 4; i++)
//    for (int j = 0; j < 4; j++)
//      graph[i][j] = g[i][j];
//  for (int i = 0; i < ZMAX; i++)
//    graph[i][i] = 1;
//  int N = ZMAX - 10;
//  while (N < ZMAX)
//  {
//
//    //printMatrix(graph, 6, 9);
//    computeDYN0(graph, N);
//    //computeDYN0(graph, N, seq);
//    //computeDYN1(graph, N, seq);
//    //computeDYN2(graph, N, seq);
//    N += 10;
//  }
//  deallocateMatrix(graph, ZMAX);
//  free(seq);
//  return 0;
//}
