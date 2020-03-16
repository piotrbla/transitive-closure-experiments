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

void computeSEQ0Pluto(int* p, int n)
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
   int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (n >= 3) {
  for (t1=0;t1<=floord(44*n-69,475);t1++) {
    lbp=max(ceild(19*t1-n+2,19),ceild(19*t1-23,44));
    ubp=min(floord(n-1,25),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(1,19*t1-19*t2);t3<=min(min(n-2,25*t2+23),19*t1-19*t2+18);t3++) {
        for (t4=max(25*t2,t3+1);t4<=min(n-1,25*t2+24);t4++) {
          for (t6=-t3+t4;t6<=t4;t6++) {
            w[(-t3+t4)][(-t3+t4)+t3] = min(w[(-t3+t4)][t6-1] + w[t6+1][(-t3+t4)+t3], w[(-t3+t4)][(-t3+t4)+t3]);;
          }
          for (t6=-t3+t4;t6<=t4;t6++) {
            w[(-t3+t4)][(-t3+t4)+t3] += p[t6];;
          }
        }
      }
    }
  }
}

  double execution_time = omp_get_wtime() - start;
  printf("pluto: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(w, n, 0);
  deallocateMatrix(w, n + 1);
}

void computeSEQ1(int* p, int n)
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
    #pragma omp parallel for
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
  printf("paralell: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(w, n, 0);
  deallocateMatrix(w, n + 1);
}

void computeSEQ2(int* p, int n)
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
for (int c0 = floord(-n + 1, 16) + 2; c0 <= floord(n - 3, 16) + 1; c0 += 1) {
    #pragma omp parallel for
    for (int c1 = -c0 - (n - 16 * c0 + 46) / 32 + 1; c1 <= min(-1, -c0); c1 += 1) {
      for (int c4 = -16 * c1 - 15; c4 <= min(min(n - 2, n + 16 * c0 + 16 * c1 - 1), -16 * c1); c4 += 1) {
        for (int c5 = max(1, -16 * c0 - 16 * c1); c5 <= min(-16 * c0 - 16 * c1 + 15, n - c4 - 1); c5 += 1) {
          for (int c7 = c5; c7 <= c4 + c5; c7 += 1) {
            w[c5][c4 + c5] = (((w[c5][c7 - 1] + w[c7 + 1][c4 + c5]) < w[c5][c4 + c5]) ? (w[c5][c7 - 1] + w[c7 + 1][c4 + c5]) : w[c5][c4 + c5]);
          }
          for (int c7 = c5; c7 <= c4 + c5; c7 += 1) {
            w[c5][c4 + c5] += p[c7];
          }
        }
      }
    }
    if (c0 == 0) {
    }
  }
  double execution_time = omp_get_wtime() - start;
  printf("stencil: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(w, n, 0);
  deallocateMatrix(w, n + 1);
}




void printMatrix(int** matrix, int N, int fileno) {
  char filename[10];
  sprintf(filename, "nontiled%d", fileno);
  FILE* f = fopen(filename, "wt");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++)
    {
      fprintf(f, "%d ", matrix[i][j]);
      //printf("%4d ", matrix[i][j]);
    }
    fprintf(f, "\n");
//    printf("\n");
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
  FILE* f = fopen("results.txt", "at");
  fprintf(f, "%d:%lf%c", n, execution_time, end_char);
  fclose(f);
}
void write_results(int n, double execution_time)
{
  write_results(n, execution_time, ';');
}


int main(void) {//vector
  const int ZMAX = 900;
  int* seq = allocateVector(ZMAX);
  for (int i = 0; i < ZMAX; i++)
    seq[i] = i;
  int N = ZMAX;
  while (N <= ZMAX)
  {
    computeSEQ0(seq, N);
    computeSEQ0Pluto(seq, N);
    computeSEQ1(seq, N);
    computeSEQ2(seq, N);
    N += 10;
  }
  free(seq);
  return 0;
}


