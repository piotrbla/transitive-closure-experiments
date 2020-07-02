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
  write_results(n, execution_time);
  printMatrix(S, n, 1);
  deallocateMatrix(S, n);
}

void computeDYN2Perfect(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();
  //Listing 1.2: Perfectly nested Nussinov loops
  for (int c0 = 1; c0 <= floord(n - 3, 8) + 1; c0 += 1) {
    for (int c3 = 8 * c0 - 7; c3 <= min(n - 2, 8 * c0); c3 += 1) {
      for (int c4 = 0; c4 < n - c3; c4 += 1) {
        for (int c5 = c3 + c4; c5 < min(n, 2 * c3 + c4); c5 += 1) {
          if (2 * c3 + c4 >= c5 + 2) {
            S[c4][c5] = max_sc(S[c4][-c3 + c5] + S[-c3 + c5 + 1][c5], S[c4][c5], S[c4 + 1][c5 - 1] + match(seq[c4], seq[c5]));
          }
          S[c4][c5] = max_sc(S[c4][c3 + c4 - 1] + S[c3 + c4][c5], S[c4][c5], S[c4 + 1][c5 - 1] + match(seq[c4], seq[c5]));
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

void computeDYN3ImperfA(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();

  for (int c0 = floord(-31 * n + 115, 3132) + 2; c0 <= floord(79 * n - 158, 2436) + 2; c0 += 1) {
    #pragma omp parallel for
    for (int c1 = max(-c0 - (n + 52) / 54 + 2, -((n + 114) / 116)); c1 <= min(min(-c0 + (n - 2) / 42 + 1, c0 + floord(-4 * c0 + 3, 31) - 1), floord(-21 * c0 + 20, 79)); c1 += 1) {
      for (int c2 = max(-c0 + c1 + floord(21 * c0 - 17 * c1 - 21, 48) + 1, -c0 - c1 - (n - 42 * c0 - 42 * c1 + 136) / 96 + 1); c2 <= min(min(-1, -c0 - c1), -((27 * c0 - 31 * c1 + 54) / 69) + 1); c2 += 1) {
        for (int c5 = max(27 * c0 - 31 * c1 + 27 * c2 - 83, -42 * c2 - 41); c5 <= min(min(n + 54 * c0 + 54 * c1 + 54 * c2 - 1, -42 * c2), 54 * c0 - 62 * c1 + 54 * c2); c5 += 1) {
          for (int c6 = max(-54 * c0 - 54 * c1 - 54 * c2, -116 * c1 - 2 * c5 - 114); c6 <= min(min(-54 * c0 - 54 * c1 - 54 * c2 + 53, n - c5 - 1), -116 * c1 - c5); c6 += 1) {
            for (int c7 = max(-116 * c1 - 115, c5 + c6); c7 <= min(min(n - 1, -116 * c1), 2 * c5 + c6 - 1); c7 += 1) {
              if (2 * c5 + c6 >= c7 + 2) {
                S[c6][c7] = max_score(S[c6][-c5 + c7] + S[-c5 + c7 + 1][c7], S[c6][c7]);
                if (c7 == c5 + c6) {
                  S[c6][c5 + c6] = max_score(S[c6][c5 + c6], S[c6 + 1][c5 + c6 - 1] + match(seq[c6], seq[c5 + c6]));
                }
              }
              S[c6][c7] = max_score(S[c6][c5 + c6 - 1] + S[c5 + c6][c7], S[c6][c7]);
              if (c7 == c5 + c6) {
                S[c6][c5 + c6] = max_score(S[c6][c5 + c6], S[c6 + 1][c5 + c6 - 1] + match(seq[c6], seq[c5 + c6]));
              }
            }
          }
        }
      }
    }
  }



  double execution_time = omp_get_wtime() - start;

  printf("IMPA: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(S, n, 3);
  deallocateMatrix(S, n);
}

void computeDYN4ImperfB(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

  double start = omp_get_wtime();
  for (int c0 = max(floord(-n - 13, 28) + 3, floord(-37 * n + 115, 2436) + 2); c0 <= floord(3 * n - 6, 116) + 2; c0 += 1) {
    #pragma omp parallel for
    for (int c1 = max(-c0 - (n + 40) / 42 + 2, -((n + 114) / 116)); c1 <= min(min(-c0 + (n - 2) / 58 + 1, floord(-c0, 3)), c0 + floord(-16 * c0 + 15, 37) - 1); c1 += 1) {
      for (int c2 = max(c1 + floord(-21 * c0 - 13 * c1 + 21, 50), -((n + 42 * c0 + 42 * c1 + 56) / 100)); c2 <= min(min(-1, -c0 - c1), -((21 * c0 - 37 * c1 + 80) / 79) + 1); c2 += 1) {
        for (int c5 = max(21 * c0 - 37 * c1 + 21 * c2 - 77, -58 * c2 - 57); c5 <= min(min(n + 42 * c0 + 42 * c1 + 42 * c2 - 1, -58 * c2), 42 * c0 - 74 * c1 + 42 * c2); c5 += 1) {
          for (int c6 = max(-42 * c0 - 42 * c1 - 42 * c2, -116 * c1 - 2 * c5 - 114); c6 <= min(min(-42 * c0 - 42 * c1 - 42 * c2 + 41, n - c5 - 1), -116 * c1 - c5); c6 += 1) {
            for (int c7 = max(-116 * c1 - 115, c5 + c6); c7 <= min(min(n - 1, -116 * c1), 2 * c5 + c6 - 1); c7 += 1) {
              if (2 * c5 + c6 >= c7 + 2) {
                S[c6][c7] = max_score(S[c6][-c5 + c7] + S[-c5 + c7 + 1][c7], S[c6][c7]);
              }
              S[c6][c7] = max_score(S[c6][c5 + c6 - 1] + S[c5 + c6][c7], S[c6][c7]);
              if (c7 == c5 + c6) {
                S[c6][c5 + c6] = max_score(S[c6][c5 + c6], S[c6 + 1][c5 + c6 - 1] + match(seq[c6], seq[c5 + c6]));
              }
            }
          }
        }
      }
    }
  }

  double execution_time = omp_get_wtime() - start;

  printf("IMPB: %lf\n", execution_time);
  write_results_full(n, execution_time, '\n');
  printMatrix(S, n, 4);
  deallocateMatrix(S, n);
}

void computeNusMP(int** table, int n, int *seq) {
  int** S = getFullCopy(table, n);

        //int c3, c1, c4, c5, c9, c10, c11;
  double start = omp_get_wtime();
for( int c1 = 1; c1 < n + floord(n - 2, 128); c1 += 1)
  #pragma omp parallel for
  for(int c3 = max(0, -n + c1 + 1); c3 <= (c1 - 1) / 129; c3 += 1)
    for(int c4 = 0; c4 <= 1; c4 += 1) {
      if (c4 == 1) {
        for(int c9 = n - c1 + 129 * c3; c9 <= min(n - 1, n - c1 + 129 * c3 + 127); c9 += 1)
          for(int c10 = max(0, n - c1 + 129 * c3 - c9 + 1); c10 <= 1; c10 += 1) {
            if (c10 == 1) {
              S[(n-c1+c3-1)][c9] = max(S[(n-c1+c3-1)][c9], S[(n-c1+c3-1)+1][c9-1] + sigma((n-c1+c3-1), c9));
            } else {
              for(int c11 = n - c1 + 129 * c3; c11 < c9; c11 += 1)
                S[(n-c1+c3-1)][c9] = max(S[(n-c1+c3-1)][c11] + S[c11+1][c9], S[(n-c1+c3-1)][c9]);
            }
          }
      } else {
        for(int c5 = 0; c5 <= 8 * c3; c5 += 1)
          for(int c9 = n - c1 + 129 * c3; c9 <= min(n - 1, n - c1 + 129 * c3 + 127); c9 += 1)
            for(int c11 = n - c1 + c3 + 16 * c5 - 1; c11 <= min(n - c1 + 129 * c3 - 1, n - c1 + c3 + 16 * c5 + 14); c11 += 1)
              S[(n-c1+c3-1)][c9] = max(S[(n-c1+c3-1)][c11] + S[c11+1][c9], S[(n-c1+c3-1)][c9]);
      }
    }
  double execution_time = omp_get_wtime() - start;

  printf("COMP: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(S, n, 0);
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
  computeNusMP(graph, N, seq);
  computeDYN1Imperfect(graph, N, seq);
  computeDYN2Perfect(graph, N, seq);
  computeDYN3ImperfA(graph, N, seq);
  computeDYN4ImperfB(graph, N, seq);
  //N += 10;
//}
  deallocateMatrix(graph, ZMAX);
  free(seq);
  return 0;
}
