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
   *  'A' => 0 -> bitowo 0001 -> 1
   *  'G' => 1 -> bitowo 0010 -> 2
   *  'C' => 2 -> bitowo 0100 -> 4
   *  'U' => 3 -> bitowo 1000 -> 8
  */
  //const bool match =
  //  (e1 == 0 && e2 == 3) || (e1 == 3 && e2 == 0) ||
  //  (e1 == 1 && e2 == 2) || (e1 == 2 && e2 == 1) ||
  //  (e1 == 1 && e2 == 3) || (e1 == 3 && e2 == 1);
  //return match;
  const int match =
    (e1 + e2 == 9) ||
    (e1 + e2 == 6) || 
    (e1 + e2 == 10) ;
  return match;
}



void printMatrix(int**, int, int);
int ** getFullCopy(int ** table, int N);
int** allocateMatrix(int);
void deallocateMatrix(int**, int);

void write_results_full(int , double , char );
void write_results(int , double );

void computeDYN2Perfect(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();
  //Listing 1.2: Perfectly nested Nussinov loops
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (n >= 3) {
  for (t1=1;t1<=n-2;t1++) {
    lbp=ceild(t1-25,26);
    ubp=floord(-t1+2*n-2,26);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(max(0,ceild(-t1+13*t2-28,30)),ceild(26*t2-n-28,30));t3<=min(floord(-t1+n-1,30),floord(-t1+26*t2+25,60));t3++) {
        if (t1 >= 2) {
          for (t4=max(30*t3,-t1+13*t2+1);t4<=min(min(-2*t1+n,30*t3+29),-t1+13*t2+13);t4++) {
            lbv=max(26*t2,t1+2*t4);
            ubv=2*t1+2*t4-2;
#pragma ivdep
#pragma vector always
            for (t5=lbv;t5<=ubv;t5++) {
              S[t4][(-t4+t5)] = max_sc(S[t4][-t1 + (-t4+t5)] + S[-t1 + (-t4+t5) + 1][(-t4+t5)], S[t4][(-t4+t5)], S[t4 + 1][(-t4+t5) - 1] + sigma(t4, (-t4+t5)));;
              S[t4][(-t4+t5)] = max_sc(S[t4][t1 + t4 - 1] + S[t1 + t4 - 1 + 1][(-t4+t5)], S[t4][(-t4+t5)], S[t4 + 1][(-t4+t5) - 1] + sigma(t4, (-t4+t5)));;
            }
            S[t4][(2*t1+t4-1)] = max_sc(S[t4][t1 + t4 - 1] + S[t1 + t4 - 1 + 1][(2*t1+t4-1)], S[t4][(2*t1+t4-1)], S[t4 + 1][(2*t1+t4-1) - 1] + sigma(t4, (2*t1+t4-1)));;
          }
        }
        for (t4=max(max(30*t3,-2*t1+n+1),26*t2-n+1);t4<=min(min(30*t3+29,-t1+n-1),-t1+13*t2+13);t4++) {
          lbv=max(26*t2,t1+2*t4);
          ubv=t4+n-1;
#pragma ivdep
#pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            S[t4][(-t4+t5)] = max_sc(S[t4][-t1 + (-t4+t5)] + S[-t1 + (-t4+t5) + 1][(-t4+t5)], S[t4][(-t4+t5)], S[t4 + 1][(-t4+t5) - 1] + sigma(t4, (-t4+t5)));;
            S[t4][(-t4+t5)] = max_sc(S[t4][t1 + t4 - 1] + S[t1 + t4 - 1 + 1][(-t4+t5)], S[t4][(-t4+t5)], S[t4 + 1][(-t4+t5) - 1] + sigma(t4, (-t4+t5)));;
          }
        }
        for (t4=max(max(30*t3,-t1+13*t2+14),26*t2-n+1);t4<=min(min(floord(-t1+26*t2+25,2),30*t3+29),-t1+n-1);t4++) {
          lbv=max(26*t2,t1+2*t4);
          ubv=min(26*t2+25,t4+n-1);
#pragma ivdep
#pragma vector always
          for (t5=lbv;t5<=ubv;t5++) {
            S[t4][(-t4+t5)] = max_sc(S[t4][-t1 + (-t4+t5)] + S[-t1 + (-t4+t5) + 1][(-t4+t5)], S[t4][(-t4+t5)], S[t4 + 1][(-t4+t5) - 1] + sigma(t4, (-t4+t5)));;
            S[t4][(-t4+t5)] = max_sc(S[t4][t1 + t4 - 1] + S[t1 + t4 - 1 + 1][(-t4+t5)], S[t4][(-t4+t5)], S[t4 + 1][(-t4+t5) - 1] + sigma(t4, (-t4+t5)));;
          }
        }
        if (t1 == 1) {
          for (t4=max(13*t2,30*t3);t4<=min(min(n-2,13*t2+12),30*t3+29);t4++) {
            S[t4][(t4+1)] = max_sc(S[t4][1 + t4 - 1] + S[1 + t4 - 1 + 1][(t4+1)], S[t4][(t4+1)], S[t4 + 1][(t4+1) - 1] + sigma(t4, (t4+1)));;
          }
        }
      }
    }
  }
}
  double execution_time = omp_get_wtime() - start;

  printf("PERF: %lf\n", execution_time);
  write_results(n, execution_time);
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
  /*
   *  'A' => 0 -> bitowo 0001 -> 1
   *  'G' => 1 -> bitowo 0010 -> 2
   *  'C' => 2 -> bitowo 0100 -> 4
   *  'U' => 3 -> bitowo 1000 -> 8
  */

  if(c=='A')    return 1;
  if(c=='G')    return 2;
  if(c=='C')    return 4;
  if(c=='U')    return 8;
  return 16;
}

#define PERFORMANCE_TEST 0

int main(void) {
#if PERFORMANCE_TEST==1
    const int ZMAX = 1600;
#else
    const int ZMAX = 16;
#endif
  int** graph = allocateMatrix(ZMAX);
  int* seq = allocateVector(ZMAX);
  for (int i = 0; i < ZMAX; i++)
    for (int j = 0; j < ZMAX; j++)
      graph[i][j] = 0;
  for (int i = 0; i < ZMAX; i++)
    graph[i][i] = 0;
  //
  const char* seqTest = "GCGUCCACGGCUAGCU";
#if PERFORMANCE_TEST==1
  for (int i=0 ; i<ZMAX ; i++)
  {
    seq[i] = 1 << (rand()%4+1);
  }
#else
  for (int i = 0; i < ZMAX; i++)
    seq[i] = getValue(seqTest[i]);
#endif
  
  int N = ZMAX - 10;
  //while (N < ZMAX)
  //{
  N += 10;
  computeDYN2Perfect(graph, N, seq);
  //N += 10;
//}
  deallocateMatrix(graph, ZMAX);
  free(seq);
  return 0;
}
