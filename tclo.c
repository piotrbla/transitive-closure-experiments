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
#pragma scop  
    for (k = 0; k < N; k++) 
    { 
        for (i = 0; i < N; i++) 
        { 
            for (j = 0; j < N; j++) 
            { 
                reach[i][j] = reach[i][j] || (reach[i][k] && reach[k][j]); 
            } 
        } 
    } 
#pragma endscop
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
