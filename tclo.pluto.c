#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#include<stdio.h> 
#include<stdlib.h> 
  
void printMatrix(int **, int); 
int **allocateMatrix(int);
void computeTC(int **matrix, int N) 
{ 
    int **reach = allocateMatrix(N);
    int i, j, k; 
    for (i = 0; i < N; i++) 
        for (j = 0; j < N; j++) 
            reach[i][j] = matrix[i][j]; 
  int t1, t2, t3, t4, t5;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (N >= 1) {
  for (t1=0;t1<=N-1;t1++) {
    for (t2=0;t2<=floord(N-1,16);t2++) {
      lbp=max(0,ceild(32*t2-N+1,32));
      ubp=min(floord(N-1,32),t2);
#pragma omp parallel for private(lbv,ubv,t4,t5)
      for (t3=lbp;t3<=ubp;t3++) {
        for (t4=32*t2-32*t3;t4<=min(N-1,32*t2-32*t3+31);t4++) {
          for (t5=32*t3;t5<=min(N-1,32*t3+31);t5++) {
            reach[t4][t5] = reach[t4][t5] || (reach[t4][t1] && reach[t1][t5]);;
          }
        }
      }
    }
  }
}
    printMatrix(reach, N); 
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
	const int N = 4;
	int **graph = allocateMatrix(N);
	int g[4][4]=	{ {1, 1, 0, 1},
                        {0, 1, 1, 0},
                        {0, 0, 1, 1},
                        {0, 0, 0, 1}
                      };
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            graph[i][j] = g [i][j];
	
    printMatrix(graph, 4);
    computeTC(graph, 4);
    return 0;
}
