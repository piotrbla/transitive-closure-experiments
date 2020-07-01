#include <omp.h>
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
//#define match(b1, b2) (((b1)+(b2)) == 3 ? 1 : 0)
#define sigma(i, j) (match(seq[i], seq[j]))
int max_score(int s1, int s2)
{
  if (s1 >= s2)
    return s1;
  return s2;
}

int max_sc(int s1, int s2, int s3) {
  if (s1>=s2 && s1>=s3)
    return s1;
  if (s2>=s3)
    return s2;
  return s3;
}

int match(const int e1, const int e2)
{
  /*
   *  'A' => 0
   *  'G' => 1
   *  'C' => 2
   *  'U' => 3
  */
  const int match =
    (e1 == 0 && e2 == 3) || (e1 == 3 && e2 == 0) ||
    (e1 == 1 && e2 == 2) || (e1 == 2 && e2 == 1) ||
    (e1 == 1 && e2 == 3) || (e1 == 3 && e2 == 1);
  return match;
  //(e1 == "A" && e2 == "U") ||
  //(e1 == "U" && e2 == "A") ||
  //(e1 == "G" && e2 == "C") ||
  //(e1 == "C" && e2 == "G") ||
  //(e1 == "G" && e2 == "U") ||
  //(e1 == "U" && e2 == "G");

}



void printMatrix(int**, int, int);
int ** getFullCopy(int ** table, int N);
int** allocateMatrix(int);
void deallocateMatrix(int**, int);

void write_results_full(int , double , char );
void write_results(int , double );

void computeDYN1Imperfect(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (n >= 2) {
  for (t2=ceild(-n-16,19);t2<=0;t2++) {
    lbp=max(0,ceild(-19*t2-41,25));
    ubp=floord(n-1,25);
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10)
    for (t4=lbp;t4<=ubp;t4++) {
      for (t5=max(max(19*t2,-n+2),-25*t4-23);t5<=min(0,19*t2+18);t5++) {
        for (t7=max(25*t4,-t5+1);t7<=min(n-1,25*t4+24);t7++) {
          for (t9=-t5;t9<=t7-1;t9++) {
            S[-t5][t7] = max_score(S[-t5][t9] + S[t9+1][t7], S[-t5][t7]);;
          }
          S[-t5][t7] = max_score(S[-t5][t7], S[-t5+1][t7-1] + sigma(-t5, t7));;
        }
      }
    }
  }
}
  double execution_time = omp_get_wtime() - start;

  printf("IMP: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(S, n, 1);
  deallocateMatrix(S, n);
}

void computeDYN2Perfect(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();
  //Listing 1.2: Perfectly nested Nussinov loops
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (n >= 2) {
  for (t2=ceild(-n-16,19);t2<=0;t2++) {
    lbp=max(0,ceild(-19*t2-41,25));
    ubp=floord(n-1,25);
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10)
    for (t4=lbp;t4<=ubp;t4++) {
      for (t5=max(max(19*t2,-n+2),-25*t4-23);t5<=min(0,19*t2+18);t5++) {
        for (t7=max(25*t4,-t5+1);t7<=min(n-1,25*t4+24);t7++) {
          for (t9=-t5;t9<=t7-1;t9++) {
            S[-t5][t7] = max_sc(S[-t5][t9] + S[t9+1][t7], S[-t5][t7], max_score(S[-t5][t7], S[-t5+1][t7-1] + sigma(-t5, t7)));;
          }
        }
      }
    }
  }
}
  double execution_time = omp_get_wtime() - start;

  printf("PER: %lf\n", execution_time);
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
  if(c=='A')    return 0;
  if(c=='G')    return 1;
  if(c=='C')    return 2;
  if(c=='U')    return 3;
  return 4;
}

int main(void) {
  const int ZMAX = 1600;
  int** graph = allocateMatrix(ZMAX);
  int* seq = allocateVector(ZMAX);
  for (int i = 0; i < ZMAX; i++)
    for (int j = 0; j < ZMAX; j++)
      graph[i][j] = 0;
  for (int i = 0; i < ZMAX; i++)
    graph[i][i] = 0;
  //
  const char* seqTest = "GCGUCCACGGCUAGCU";
  ///////////////////////GCGUCCACGGCUAGCU
  for (int i=0 ; i<ZMAX ; i++)
    seq[i] = rand()%4;
  //for (int i = 0; i < ZMAX; i++)
  //  seq[i] = getValue(seqTest[i]);
  
  int N = ZMAX - 10;
  //while (N < ZMAX)
  //{
  N += 10;
  computeDYN1Imperfect(graph, N, seq);
  computeDYN2Perfect(graph, N, seq);
  //N += 10;
//}
  deallocateMatrix(graph, ZMAX);
  free(seq);
  return 0;
}
