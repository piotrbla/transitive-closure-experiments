#pragma warning(disable : 4996)
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
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

  for (int c0 = floord(-n + 113, 1680) + 1; c0 <= floord(n - 2, 120) + 1; c0 += 1) {
    #pragma omp parallel for
    for (int c1 = max(-c0 - (n + 110) / 112 + 1, -((n + 118) / 120)); c1 <= min(14 * c0 - 1, -c0); c1 += 1) {
      for (int c4 = max(max(-n + 2, 112 * c0 + 112 * c1 - 111), 120 * c1 + 1); c4 <= 112 * c0 + 112 * c1; c4 += 1) {
        for (int c5 = max(-120 * c1 - 119, -c4 + 1); c5 <= min(n - 1, -120 * c1); c5 += 1) {
          for (int c7 = -c4; c7 < c5; c7 += 1) {
            S[-c4][c5] = max_score(S[-c4][c7] + S[c7 + 1][c5], S[-c4][c5]);
          }
          S[-c4][c5] = max_score(S[-c4][c5], S[-c4 + 1][c5 - 1] + match(seq[-c4], seq[c5]));
        }
      }
    }
  }
  double execution_time = omp_get_wtime() - start;

  printf("IMPE: %lf\n", execution_time);
  write_results_full(n, execution_time, '\n');
  printMatrix(S, n, 1);
  deallocateMatrix(S, n);
}

void computeDYN2Perfect(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();
  //Listing 1.2: Perfectly nested Nussinov loops
#pragma scop
  for (int c0 = 1; c0 < n - 1; c0 += 1)
  {
    for (int c1 = 0; c1 < n - c0; c1 += 1)
    {
      for (int c2 = c0 + c1; c2 < min(n, 2 * c0 + c1); c2 += 1) {
        if (2 * c0 + c1 >= c2 + 2)
        {
          // (c0, c1, c2, -c0 + c2);
          S[c1][c2] = max_sc(S[c1][-c0 + c2] + S[-c0 + c2 + 1][c2], S[c1][c2], S[c1 + 1][c2 - 1] + sigma(c1, c2));
        }
        //(c0, c1, c2, c0 + c1 - 1);
        S[c1][c2] = max_sc(S[c1][c0 + c1 - 1] + S[c0 + c1 - 1 + 1][c2], S[c1][c2], S[c1 + 1][c2 - 1] + sigma(c1, c2));
      }
    }
  }
#pragma endscop
  double execution_time = omp_get_wtime() - start;

  printf("PERF: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(S, n, 2);
  deallocateMatrix(S, n);
}

void computeDYN3ImperfA(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();
  //Listing 1.2: Perfectly nested Nussinov loops
//c1-> "i" (drugi element pseudo instrukcji)
//trzeci element każdej  pseudo instrukcji jest to "j"
//czwarty element jest to "k"
//piąty element jeśli jest "0" oznacza "s1", 
//natomiast jeśli "1" oznacza, ze maja byc  2 instrukcje wstawione - s1 i s2  
#pragma scop
  for (int c0 = 1; c0 < n; c0 += 1)
  {
      for (int c1 = 0; c1 < n - c0; c1 += 1)
      {
          for (int c2 = c0 + c1; c2 < min(n, 2 * c0 + c1); c2 += 1)
          {
              if (2 * c0 + c1 >= c2 + 2)
              {
                  S[c1][c2] = max_score(S[c1][-c0 + c2] + S[-c0 + c2+1][c2], S[c1][c2]); // s1
                  //(c0, c1, c2, -c0 + c2, 0);
                  if (c2 == c0 + c1)
                  {
                      S[c1][c2] = max_score(S[c1][c2], S[c1+1][c2-1] + sigma(c1, c2)); // s2
                      //(c0, c1, c0 + c1, c1, 1);
                  }
              }
              S[c1][c2] = max_score(S[c1][c0 + c1 - 1] + S[c0 + c1 - 1 + 1][c2], S[c1][c2]); // s1
              //(c0, c1, c2, c0 + c1 - 1, 0);
              if (c2 == c0 + c1)
              {
                  S[c1][c2] = max_score(S[c1][c2], S[c1+1][c2-1] + sigma(c1, c2)); // s2
                  //(c0, c1, c0 + c1, c0 + c1 - 1, 1);
              }
          }
      }
  }
#pragma endscop
  double execution_time = omp_get_wtime() - start;

  printf("IMPA: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(S, n, 3);
  deallocateMatrix(S, n);
}

void computeDYN4ImperfB(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();
  //Listing 1.2: Perfectly nested Nussinov loops
#pragma scop
  for (int c0 = 1; c0 < n; c0 += 1)
  {
      for (int c1 = 0; c1 < n - c0; c1 += 1)
      {
          for (int c2 = c0 + c1; c2 < min(n, 2 * c0 + c1); c2 += 1)
          {
              if (2 * c0 + c1 >= c2 + 2)
              {
                  S[c1][c2] = max_score(S[c1][-c0 + c2] + S[-c0 + c2 + 1][c2], S[c1][c2]); // s1
                  //(c0, c1, c2, -c0 + c2, 0);
              }
              S[c1][c2] = max_score(S[c1][c0 + c1 - 1] + S[c0 + c1 - 1 + 1][c2], S[c1][c2]); // s1
              //(c0, c1, c2, c0 + c1 - 1, 0);
              if (c2 == c0 + c1)
              {
                  S[c1][c0 + c1] = max_score(S[c1][c0 + c1], S[c1+1][c0 + c1 - 1] + sigma(c1, c0 + c1)); // s2
                  //(c0, c1, c0 + c1, c0 + c1 - 1, 1);
              }
          }
      }
  }

#pragma endscop
  double execution_time = omp_get_wtime() - start;

  printf("IMPB: %lf\n", execution_time);
  write_results_full(n, execution_time, '\n');
  printMatrix(S, n, 4);
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
  setlocale(LC_NUMERIC, "Polish");
  fprintf(f, "%d;%lf%c", n, execution_time, end_char);
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

#define PERFORMANCE_TEST 1

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
  computeDYN1Imperfect(graph, N, seq);
  //computeDYN2Perfect(graph, N, seq);
  //computeDYN3ImperfA(graph, N, seq);
  //computeDYN4ImperfB(graph, N, seq);
  //N += 10;
//}
  deallocateMatrix(graph, ZMAX);
  free(seq);
  return 0;
}
