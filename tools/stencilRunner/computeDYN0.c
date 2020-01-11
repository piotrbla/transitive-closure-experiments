#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S0(a, i, j, k) c[i][j] = c[i][k] + c[k][j]

void printMatrix(int **, int, int);
int **allocateMatrix(int);

void computeDYN0(int **matrix, int n) {
  int **c = allocateMatrix(n + 1);
  int i, j, k, a;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      c[i][j] = matrix[i][j];
  double start = omp_get_wtime();
#pragma scop
  for (int c0 = 2; c0 < n; c0 += 1)
    for (int c1 = 1; c1 <= n - c0; c1 += 1)
      for (int c2 = c0 + c1; c2 <= min(n, 2 * c0 + c1 - 2); c2 += 1) {
	if (2 * c0 + c1 >= c2 + 3)
	  S0(c0, c1, c2, -c0 + c2 + 1);
	S0(c0, c1, c2, c0 + c1 - 1);
  }
#pragma endscop
  printf("nontiled: %lf\n", omp_get_wtime() - start);
  printMatrix(c, n, 0);
}
