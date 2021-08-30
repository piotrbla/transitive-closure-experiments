#include <omp.h>
//kompilacja np.: gcc -DMYTHREADS=32 -o bstexe_32 -O2 bst.cpp -fopenmp -lm
//kompilacja np.: gcc -DMYTHREADS=16 -o bstexe_16 -O2 bst.cpp -fopenmp -lm

#pragma warning(disable : 4996)
#include <math.h>
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
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (n >= 3) {
  for (t2=-n+2;t2<=-1;t2++) {
    lbp=-t2+2;
    ubp=n;
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10)
    for (t4=lbp;t4<=ubp;t4++) {
      for (t8=-t2+1;t8<=t4-1-7;t8+=8) {
        c[-t2][t4] = min(c[-t2][t4], w[-t2][t4]+c[-t2][t8]+c[t8][t4]);;
        c[-t2][t4] = min(c[-t2][t4], w[-t2][t4]+c[-t2][(t8+1)]+c[(t8+1)][t4]);;
        c[-t2][t4] = min(c[-t2][t4], w[-t2][t4]+c[-t2][(t8+2)]+c[(t8+2)][t4]);;
        c[-t2][t4] = min(c[-t2][t4], w[-t2][t4]+c[-t2][(t8+3)]+c[(t8+3)][t4]);;
        c[-t2][t4] = min(c[-t2][t4], w[-t2][t4]+c[-t2][(t8+4)]+c[(t8+4)][t4]);;
        c[-t2][t4] = min(c[-t2][t4], w[-t2][t4]+c[-t2][(t8+5)]+c[(t8+5)][t4]);;
        c[-t2][t4] = min(c[-t2][t4], w[-t2][t4]+c[-t2][(t8+6)]+c[(t8+6)][t4]);;
        c[-t2][t4] = min(c[-t2][t4], w[-t2][t4]+c[-t2][(t8+7)]+c[(t8+7)][t4]);;
      }
      for (;t8<=t4-1;t8++) {
        c[-t2][t4] = min(c[-t2][t4], w[-t2][t4]+c[-t2][t8]+c[t8][t4]);;
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
    N += 100;
  }
  free(seq);
  deallocateMatrix(graph, ZMAX);
  return 0;
}
                              
