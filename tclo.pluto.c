#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#include<stdio.h> 
  
#define V 4 
void printMatrix(int matrix[][V]); 

void computeTC(int graph[][V]) 
{ 
    int reach[V][V], i, j, k; 
    for (i = 0; i < V; i++) 
        for (j = 0; j < V; j++) 
            reach[i][j] = graph[i][j]; 
  int t1, t2, t3, t4, t5;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (V >= 1) {
  for (t1=0;t1<=V-1;t1++) {
    for (t2=0;t2<=floord(V-1,16);t2++) {
      lbp=max(0,ceild(32*t2-V+1,32));
      ubp=min(floord(V-1,32),t2);
#pragma omp parallel for private(lbv,ubv,t4,t5)
      for (t3=lbp;t3<=ubp;t3++) {
        for (t4=32*t2-32*t3;t4<=min(V-1,32*t2-32*t3+31);t4++) {
          for (t5=32*t3;t5<=min(V-1,32*t3+31);t5++) {
            reach[t4][t5] = reach[t4][t5] || (reach[t4][t1] && reach[t1][t5]);;
          }
        }
      }
    }
  }
}
    printMatrix(reach); 
} 

void printMatrix(int matrix[][V])
{
    for (int i = 0; i < V; i++)
    {
        for (int j = 0; j < V; j++)
            printf ("%d ", matrix[i][j]);
        printf("\n");
    }
}

int main(void)
{
    int graph[V][V] = { {1, 1, 0, 1},
                        {0, 1, 1, 0},
                        {0, 0, 1, 1},
                        {0, 0, 0, 1}
                      };
    printMatrix(graph);
    computeTC(graph);
    return 0;
}
