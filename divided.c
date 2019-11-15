#include <omp.h>
#include <math.h>
#include<stdio.h> 
#include<stdlib.h> 
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S1(i, j, k) reach[i][j] = reach[i][j] || (reach[i][k] && reach[k][j])
void printMatrix(int **, int); 
int **allocateMatrix(int);

void computeTC(int **matrix, int n)
{
    int **reach = allocateMatrix(n+1);
    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            reach[i][j] = matrix[i][j];
//#pragma omp parallel for
  for (int c0 = 0; c0 < floord(n + 1, 2) - 1; c0 += 1)
    for (int c1 = 0; c1 <= min(c0 + 1, n / 2 - 1); c1 += 1)
      for (int c2 = 0; c2 <= min(c0 - c1 + 1, n / 2 - 1); c2 += 1)
        for (int c4 = max(2 * c0 + 1, 2 * c1 + 2 * c2); c4 <= 2 * c0 + 2; c4 += 1)
          for (int c5 = 2 * c1 + 1; c5 <= min(2 * c1 + 2, -2 * c2 + c4 + 1); c5 += 1)
            for (int c6 = 2 * c2 + 1; c6 <= min(2 * c2 + 2, c4 - c5 + 2); c6 += 1)
              reach[c5][c6] = reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] && reach[c4 - c5 - c6 + 2][c6]);
//#pragma omp parallel for
  for (int c0 = floord(n + 1, 2) - 1; c0 < 2 * n - 2; c0 += 1)
    for (int c1 = 0; c1 < n / 2; c1 += 1)
      for (int c2 = max(0, -2 * n + c0 - c1 + (n + 1) / 2); c2 <= min(c0 - c1 + 1, n / 2 - 1); c2 += 1) {
        if (3 * n + 2 * c1 == 2 * c0 + 2 && 2 * c2 + 2 == n)
          for (int c5 = -3 * n + 2 * c0 + 3; c5 <= -3 * n + 2 * c0 + 4; c5 += 1)
            reach[c5][n - 1] = reach[c5][n - 1] || (reach[c5][-3 * n + 2 * c0 - c5 + 4] && reach[-3 * n + 2 * c0 - c5 + 4][n - 1]);
        for (int c4 = max(2 * c0 + 1, 2 * c1 + 2 * c2); c4 <= min(2 * c0 + 2, n + 2 * c1 + 2 * c2 - 1); c4 += 1)
          for (int c5 = 2 * c1 + 1; c5 <= min(min(n - 1, 2 * c1 + 2), -2 * c2 + c4 + 1); c5 += 1)
            for (int c6 = 2 * c2 + 1; c6 <= min(min(n - 1, 2 * c2 + 2), c4 - c5 + 2); c6 += 1)
              reach[c5][c6] = reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] && reach[c4 - c5 - c6 + 2][c6]);
        if (2 * n + c2 >= c0 + c1 + 3) {
          for (int c4 = max(2 * c0 + 1, n + 2 * c1 + 2 * c2); c4 <= min(min(min(4 * n - 5, 2 * c0 + 2), 3 * n + 2 * c1 - 2), n + c0 + c1 + c2 + 1); c4 += 1) {
            if (n + c0 + c1 + c2 >= c4) {
              for (int c5 = 2 * c1 + 1; c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
                for (int c6 = 2 * c2 + 1; c6 <= min(n - 1, 2 * c2 + 2); c6 += 1) {
                  if (n + c5 + c6 >= c4 + 3) {
                    reach[c5][c6] = reach[c5][c6] || (reach[c5][c4 - c5 - c6 + 2] && reach[c4 - c5 - c6 + 2][c6]);
                  } else
                    reach[c5][c6] = reach[c5][c6] || (reach[c5][-n + c4 - c5 - c6 + 2] && reach[-n + c4 - c5 - c6 + 2][c6]);
                }
            } else
              for (int c5 = 2 * c1 + 1; c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
                for (int c6 = 2 * c2 + 1; c6 <= min(n - 1, 2 * c2 + 2); c6 += 1) {
                  if (n + c5 + c6 >= c0 + c1 + c2 + 4) {
                    reach[c5][c6] = reach[c5][c6] || (reach[c5][c0 + c1 + c2 - c5 - c6 + 3] && reach[c0 + c1 + c2 - c5 - c6 + 3][c6]);
                  } else
                    reach[c5][c6] = reach[c5][c6] || (reach[c5][-n + c0 + c1 + c2 - c5 - c6 + 3] && reach[-n + c0 + c1 + c2 - c5 - c6 + 3][c6]);
                }
          }
          if (c0 + 3 == 2 * n && 2 * c1 + 2 == n && 2 * c2 + 2 == n)
            reach[n - 1][n - 1] = reach[n - 1][n - 1] || (reach[n - 1][0] && reach[0][n - 1]);
        } else if (2 * c1 + 2 == n && 3 * n + 2 * c2 == 2 * c0 + 2)
          for (int c6 = -3 * n + 2 * c0 + 3; c6 <= -3 * n + 2 * c0 + 4; c6 += 1)
            reach[n - 1][c6] = reach[n - 1][c6] || (reach[n - 1][-3 * n + 2 * c0 - c6 + 4] && reach[-3 * n + 2 * c0 - c6 + 4][c6]);
        for (int c4 = max(2 * c0 + 1, 2 * n + 2 * c1 + 2 * c2 + 2); c4 <= min(2 * c0 + 2, 3 * n + 2 * c1 + 2 * c2 + 1); c4 += 1)
          for (int c5 = max(2 * c1 + 1, -3 * n - 2 * c2 + c4 + 1); c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
            for (int c6 = max(2 * c2 + 1, -3 * n + c4 - c5 + 3); c6 <= min(n - 1, 2 * c2 + 2); c6 += 1)
              reach[c5][c6] = reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] && reach[-2 * n + c4 - c5 - c6 + 2][c6]);
      }
  for (int c0 = 2 * n - 2; c0 < 2 * n + floord(n, 2) - 2; c0 += 1)
    for (int c1 = -2 * n + c0 + 2; c1 < n / 2; c1 += 1)
      for (int c2 = max(-2 * n + c0 + 2, -2 * n + c0 - c1 + (n + 1) / 2); c2 < n / 2; c2 += 1)
        for (int c4 = 2 * c0 + 1; c4 <= min(min(5 * n - 5, 2 * c0 + 2), 3 * n + 2 * c1 + 2 * c2 + 1); c4 += 1)
          for (int c5 = max(max(2 * c1 + 1, -3 * n - 2 * c2 + c4 + 1), -4 * n + c4 + 4); c5 <= min(n - 1, 2 * c1 + 2); c5 += 1)
            for (int c6 = max(2 * c2 + 1, -3 * n + c4 - c5 + 3); c6 <= min(n - 1, 2 * c2 + 2); c6 += 1)
              reach[c5][c6] = reach[c5][c6] || (reach[c5][-2 * n + c4 - c5 - c6 + 2] && reach[-2 * n + c4 - c5 - c6 + 2][c6]);
    printMatrix(reach, 6);
}

void printMatrix(int **matrix, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            printf ("%d ", matrix[i][j]);
        printf("\n");
    }
}
int **allocateMatrix(int N)
{
	int **t = (int **) malloc(sizeof(int *) * N);
	for (int i = 0 ; i < N ; i++)
	{
		t[i] = (int *) malloc(sizeof(int) * N);
	}
	return t;
}
int main(void)
{
	const int N = 400;
	int **graph = allocateMatrix(N);
	int g[4][4]=	{ {1, 1, 0, 1},
                        {0, 1, 1, 0},
                        {0, 0, 1, 1},
                        {0, 0, 0, 1}
                      };
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            graph[i][j] = g [i][j];
    for (int i = 0; i < N; i++)
        graph[i][i] = 1;
	
    printMatrix(graph, 6);
    computeTC(graph, N);
    return 0;
}
