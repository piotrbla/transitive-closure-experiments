#include<stdio.h> 
  
#define V 4 
void printMatrix(int matrix[][V]); 

void computeTC(int graph[][V]) 
{ 
    int reach[V][V], i, j, k; 
    for (i = 0; i < V; i++) 
        for (j = 0; j < V; j++) 
            reach[i][j] = graph[i][j]; 
#pragma scop  
    for (k = 0; k < V; k++) 
    { 
        for (i = 0; i < V; i++) 
        { 
            for (j = 0; j < V; j++) 
            { 
                reach[i][j] = reach[i][j] || (reach[i][k] && reach[k][j]); 
            } 
        } 
    } 
#pragma endscop
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
