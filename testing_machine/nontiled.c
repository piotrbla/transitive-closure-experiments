//
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define S0(a, i, j, k) c[i][j] = c[i][k] + c[k][j]
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
  const int match =
    (e1 + e2 == 9) ||
    (e1 + e2 == 6) || 
    (e1 + e2 == 10) ;
  return match;
}
void LMKernel6_03_Modification(int loop, int n, int *input_w, int **input_b)
{
  int** b = get_full_copy(input_b, n);
  int* w = get_vector_copy(input_w, n);

  double start = omp_get_wtime();
#pragma scop
for (int c0 = 1; c0 <= loop; c0 += 1)
  for (int c1 = 1; c1 < n; c1 += 1)
    for (int c2 = c1; c2 < n; c2 += 1)
       w[c2] += b[-c1 + c2][c2] * w[(c2 - (-c1 + c2)) - 1];
#pragma endscop

  double execution_time = omp_get_wtime() - start;

  printf("MOD_03: %lf\n", execution_time);
  write_results(n, execution_time);
  print_vector(w, n);
  deallocate_matrix(b, n);
  free(w);
  return;
}


