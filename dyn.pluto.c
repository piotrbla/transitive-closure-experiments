#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define S0(a, i, j, k) d[i][j] = c[i][k] + c[k][j]

void printMatrix(int**, int, int);
int** allocateMatrix(int);
void deallocateMatrix(int**, int);

void write_results(int , double , char );
void write_results(int , double );

void computeDYN0(int** matrix, int n) {
  int** c = allocateMatrix(n + 1);
  int** d = allocateMatrix(n + 1);
  int i, j;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      c[i][j] = matrix[i][j];
  double start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if (n >= 3) {
  lbp=0;
  ubp=floord(n-1,19);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6)
  for (t1=lbp;t1<=ubp;t1++) {
    for (t2=0;t2<=min(floord(n-2,25),floord(-19*t1+n,25));t2++) {
      for (t3=max(max(ceild(19*t1-27,29),ceild(25*t2-26,29)),ceild(19*t1+25*t2-28,29));t3<=min(floord(n,29),floord(38*t1+25*t2+58,29));t3++) {
        if ((t1 <= floord(-25*t2+29*t3-22,38)) && (t2 <= floord(29*t3-28,25))) {
          if ((t2+t3)%2 == 0) {
            S0(((-25*t2+29*t3-22)/2), (25*t2+24), 29*t3, ((-25*t2+29*t3-22)/2) + (25*t2+24) - 1);;
          }
        }
        if ((t1 == 0) && (t2 >= ceild(29*t3-26,25))) {
          for (t5=max(max(1,25*t2),29*t3-2);t5<=min(min(n-2,25*t2+24),29*t3+26);t5++) {
            S0(2, t5, (t5+2), 2 + t5 - 1);;
          }
        }
        for (t4=max(max(3,ceild(-25*t2+29*t3-21,2)),19*t1);t4<=min(min(min(floord(29*t3+1,2),floord(-25*t2+n-22,2)),floord(-25*t2+29*t3+2,2)),19*t1+18);t4++) {
          S0(t4, (29*t3-2*t4+2), 29*t3, t4 + (29*t3-2*t4+2) - 1);;
          for (t5=29*t3-2*t4+3;t5<=25*t2+24;t5++) {
            lbv=max(29*t3,t4+t5);
            ubv=2*t4+t5-3;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              S0(t4, t5, t6, -t4 + t6 + 1);;
              S0(t4, t5, t6, t4 + t5 - 1);;
            }
            S0(t4, t5, (2*t4+t5-2), t4 + t5 - 1);;
          }
        }
        if (29*t3 == n) {
          for (t4=max(max(3,ceild(-25*t2+n-21,2)),19*t1);t4<=min(min(floord(n+1,2),floord(-25*t2+n+2,2)),19*t1+18);t4++) {
            if (n%29 == 0) {
              S0(t4, (-2*t4+n+2), n, t4 + (-2*t4+n+2) - 1);;
            }
            for (t5=-2*t4+n+3;t5<=min(25*t2+24,-t4+n);t5++) {
              if (n%29 == 0) {
                S0(t4, t5, n, -t4 + n + 1);;
              }
              if (n%29 == 0) {
                S0(t4, t5, n, t4 + t5 - 1);;
              }
            }
          }
        }
        if ((t1 <= floord(29*t3+2,38)) && (t1 >= ceild(29*t3-34,38)) && (t2 == 0) && (t3 >= 2) && (t3 <= floord(n-24,29))) {
          if (t3%2 == 0) {
            for (t5=1;t5<=24;t5++) {
              lbv=29*t3;
              ubv=29*t3+t5-1;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(((29*t3+2)/2), t5, t6, -((29*t3+2)/2) + t6 + 1);;
                S0(((29*t3+2)/2), t5, t6, ((29*t3+2)/2) + t5 - 1);;
              }
              S0(((29*t3+2)/2), t5, (29*t3+t5), ((29*t3+2)/2) + t5 - 1);;
            }
          }
        }
        if (t3 <= floord(n-1,29)) {
          for (t4=max(max(3,ceild(-25*t2+n-21,2)),19*t1);t4<=min(min(floord(29*t3+1,2),floord(-25*t2+29*t3+2,2)),19*t1+18);t4++) {
            S0(t4, (29*t3-2*t4+2), 29*t3, t4 + (29*t3-2*t4+2) - 1);;
            for (t5=29*t3-2*t4+3;t5<=-2*t4+n+2;t5++) {
              lbv=max(29*t3,t4+t5);
              ubv=2*t4+t5-3;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(t4, t5, t6, -t4 + t6 + 1);;
                S0(t4, t5, t6, t4 + t5 - 1);;
              }
              S0(t4, t5, (2*t4+t5-2), t4 + t5 - 1);;
            }
            for (t5=-2*t4+n+3;t5<=min(25*t2+24,-t4+n);t5++) {
              lbv=max(29*t3,t4+t5);
              ubv=n;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(t4, t5, t6, -t4 + t6 + 1);;
                S0(t4, t5, t6, t4 + t5 - 1);;
              }
            }
          }
        }
        if (t3 <= floord(n-28,29)) {
          for (t4=max(max(3,ceild(-25*t2+29*t3+3,2)),19*t1);t4<=min(floord(-25*t2+29*t3+6,2),19*t1+18);t4++) {
            for (t5=max(1,25*t2);t5<=25*t2+24;t5++) {
              lbv=max(29*t3,t4+t5);
              ubv=2*t4+t5-3;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(t4, t5, t6, -t4 + t6 + 1);;
                S0(t4, t5, t6, t4 + t5 - 1);;
              }
              S0(t4, t5, (2*t4+t5-2), t4 + t5 - 1);;
            }
          }
        }
        if (t3 >= ceild(n-27,29)) {
          for (t4=max(max(3,ceild(-25*t2+29*t3+3,2)),19*t1);t4<=min(floord(-25*t2+n-22,2),19*t1+18);t4++) {
            for (t5=max(1,25*t2);t5<=25*t2+24;t5++) {
              lbv=max(29*t3,t4+t5);
              ubv=2*t4+t5-3;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(t4, t5, t6, -t4 + t6 + 1);;
                S0(t4, t5, t6, t4 + t5 - 1);;
              }
              S0(t4, t5, (2*t4+t5-2), t4 + t5 - 1);;
            }
          }
        }
        if ((t1 <= floord(29*t3+2,38)) && (t1 >= ceild(29*t3-34,38)) && (t2 == 0) && (t3 <= floord(n-1,29)) && (t3 >= max(2,ceild(n-23,29)))) {
          if (t3%2 == 0) {
            for (t5=1;t5<=-29*t3+n;t5++) {
              lbv=29*t3;
              ubv=29*t3+t5-1;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(((29*t3+2)/2), t5, t6, -((29*t3+2)/2) + t6 + 1);;
                S0(((29*t3+2)/2), t5, t6, ((29*t3+2)/2) + t5 - 1);;
              }
              S0(((29*t3+2)/2), t5, (29*t3+t5), ((29*t3+2)/2) + t5 - 1);;
            }
            for (t5=-29*t3+n+1;t5<=24;t5++) {
              lbv=29*t3;
              ubv=n;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(((29*t3+2)/2), t5, t6, -((29*t3+2)/2) + t6 + 1);;
                S0(((29*t3+2)/2), t5, t6, ((29*t3+2)/2) + t5 - 1);;
              }
            }
          }
        }
        for (t4=max(max(max(3,ceild(-25*t2+n-21,2)),ceild(-25*t2+29*t3+3,2)),19*t1);t4<=min(min(min(floord(n+1,2),floord(-25*t2+n+2,2)),19*t1+18),29*t3-n+30);t4++) {
          for (t5=max(1,25*t2);t5<=-2*t4+n+2;t5++) {
            lbv=max(29*t3,t4+t5);
            ubv=2*t4+t5-3;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              S0(t4, t5, t6, -t4 + t6 + 1);;
              S0(t4, t5, t6, t4 + t5 - 1);;
            }
            S0(t4, t5, (2*t4+t5-2), t4 + t5 - 1);;
          }
          for (t5=-2*t4+n+3;t5<=min(25*t2+24,-t4+n);t5++) {
            lbv=max(29*t3,t4+t5);
            ubv=n;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              S0(t4, t5, t6, -t4 + t6 + 1);;
              S0(t4, t5, t6, t4 + t5 - 1);;
            }
          }
        }
        for (t4=max(max(max(ceild(-25*t2+n-21,2),ceild(-25*t2+29*t3+3,2)),19*t1),29*t3-n+31);t4<=min(min(min(floord(n+1,2),floord(-25*t2+n+2,2)),floord(-25*t2+29*t3+6,2)),19*t1+18);t4++) {
          for (t5=max(1,25*t2);t5<=-2*t4+n+2;t5++) {
            lbv=max(29*t3,t4+t5);
            ubv=2*t4+t5-3;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              S0(t4, t5, t6, -t4 + t6 + 1);;
              S0(t4, t5, t6, t4 + t5 - 1);;
            }
            S0(t4, t5, (2*t4+t5-2), t4 + t5 - 1);;
          }
          for (t5=-2*t4+n+3;t5<=25*t2+24;t5++) {
            lbv=max(29*t3,t4+t5);
            ubv=n;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              S0(t4, t5, t6, -t4 + t6 + 1);;
              S0(t4, t5, t6, t4 + t5 - 1);;
            }
          }
        }
        if (t3 <= floord(n-28,29)) {
          for (t4=max(max(3,ceild(-25*t2+29*t3+7,2)),19*t1);t4<=min(min(floord(29*t3+29,2),floord(-25*t2+29*t3+30,2)),19*t1+18);t4++) {
            for (t5=max(1,25*t2);t5<=29*t3-2*t4+30;t5++) {
              lbv=max(29*t3,t4+t5);
              ubv=2*t4+t5-3;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(t4, t5, t6, -t4 + t6 + 1);;
                S0(t4, t5, t6, t4 + t5 - 1);;
              }
              S0(t4, t5, (2*t4+t5-2), t4 + t5 - 1);;
            }
            for (t5=29*t3-2*t4+31;t5<=min(25*t2+24,29*t3-t4+28);t5++) {
              lbv=max(29*t3,t4+t5);
              ubv=29*t3+28;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(t4, t5, t6, -t4 + t6 + 1);;
                S0(t4, t5, t6, t4 + t5 - 1);;
              }
            }
          }
        }
        if ((t1 <= floord(n+2,38)) && (t1 >= ceild(n-34,38)) && (t2 == 0) && (t3 >= ceild(3*n-58,58))) {
          if (n%2 == 0) {
            for (t5=1;t5<=min(24,floord(n-2,2));t5++) {
              lbv=max(ceild(2*t5+n+2,2),29*t3);
              ubv=n;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(((n+2)/2), t5, t6, -((n+2)/2) + t6 + 1);;
                S0(((n+2)/2), t5, t6, ((n+2)/2) + t5 - 1);;
              }
            }
          }
        }
        if ((t1 <= floord(n+2,38)) && (t1 >= ceild(n-34,38)) && (t2 == 0) && (t3 <= floord(3*n-60,58)) && (t3 >= ceild(n-4,29))) {
          if (n%2 == 0) {
            for (t5=1;t5<=24;t5++) {
              lbv=29*t3;
              ubv=n;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(((n+2)/2), t5, t6, -((n+2)/2) + t6 + 1);;
                S0(((n+2)/2), t5, t6, ((n+2)/2) + t5 - 1);;
              }
            }
          }
        }
        if (t3 >= ceild(n-27,29)) {
          for (t4=max(max(ceild(-25*t2+29*t3+7,2),19*t1),29*t3-n+31);t4<=min(min(floord(n+1,2),floord(-25*t2+n+2,2)),19*t1+18);t4++) {
            for (t5=max(1,25*t2);t5<=-2*t4+n+2;t5++) {
              lbv=max(29*t3,t4+t5);
              ubv=2*t4+t5-3;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(t4, t5, t6, -t4 + t6 + 1);;
                S0(t4, t5, t6, t4 + t5 - 1);;
              }
              S0(t4, t5, (2*t4+t5-2), t4 + t5 - 1);;
            }
            for (t5=-2*t4+n+3;t5<=29*t3-2*t4+30;t5++) {
              lbv=max(29*t3,t4+t5);
              ubv=n;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(t4, t5, t6, -t4 + t6 + 1);;
                S0(t4, t5, t6, t4 + t5 - 1);;
              }
            }
            for (t5=29*t3-2*t4+31;t5<=min(25*t2+24,-t4+n);t5++) {
              lbv=max(29*t3,t4+t5);
              ubv=n;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(t4, t5, t6, -t4 + t6 + 1);;
                S0(t4, t5, t6, t4 + t5 - 1);;
              }
            }
          }
        }
        for (t4=max(ceild(-25*t2+n+3,2),19*t1);t4<=min(min(min(n-1,19*t1+18),-25*t2+n),29*t3-n+30);t4++) {
          for (t5=max(1,25*t2);t5<=min(25*t2+24,-t4+n);t5++) {
            lbv=max(29*t3,t4+t5);
            ubv=n;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              S0(t4, t5, t6, -t4 + t6 + 1);;
              S0(t4, t5, t6, t4 + t5 - 1);;
            }
          }
        }
        for (t4=max(max(ceild(-25*t2+n+3,2),19*t1),29*t3-n+31);t4<=min(floord(-25*t2+29*t3+6,2),19*t1+18);t4++) {
          for (t5=max(1,25*t2);t5<=25*t2+24;t5++) {
            lbv=max(29*t3,t4+t5);
            ubv=n;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              S0(t4, t5, t6, -t4 + t6 + 1);;
              S0(t4, t5, t6, t4 + t5 - 1);;
            }
          }
        }
        if ((t1 <= floord(n+2,38)) && (t1 >= ceild(n-34,38)) && (t2 == 0) && (t3 <= min(floord(n-5,29),floord(3*n-60,58))) && (t3 >= ceild(n-27,29))) {
          if (n%2 == 0) {
            for (t5=1;t5<=29*t3-n+28;t5++) {
              lbv=max(ceild(2*t5+n+2,2),29*t3);
              ubv=n;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(((n+2)/2), t5, t6, -((n+2)/2) + t6 + 1);;
                S0(((n+2)/2), t5, t6, ((n+2)/2) + t5 - 1);;
              }
            }
            for (t5=29*t3-n+29;t5<=min(24,floord(n-2,2));t5++) {
              lbv=max(ceild(2*t5+n+2,2),29*t3);
              ubv=n;
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(((n+2)/2), t5, t6, -((n+2)/2) + t6 + 1);;
                S0(((n+2)/2), t5, t6, ((n+2)/2) + t5 - 1);;
              }
            }
          }
        }
        for (t4=max(max(max(ceild(-25*t2+n+3,2),ceild(-25*t2+29*t3+7,2)),19*t1),29*t3-n+31);t4<=min(min(floord(29*t3+29,2),floord(-25*t2+29*t3+30,2)),19*t1+18);t4++) {
          for (t5=max(1,25*t2);t5<=29*t3-2*t4+30;t5++) {
            lbv=max(29*t3,t4+t5);
            ubv=n;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              S0(t4, t5, t6, -t4 + t6 + 1);;
              S0(t4, t5, t6, t4 + t5 - 1);;
            }
          }
          for (t5=29*t3-2*t4+31;t5<=min(25*t2+24,-t4+n);t5++) {
            lbv=max(29*t3,t4+t5);
            ubv=n;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              S0(t4, t5, t6, -t4 + t6 + 1);;
              S0(t4, t5, t6, t4 + t5 - 1);;
            }
          }
        }
        if ((t1 <= floord(29*t3+30,38)) && (t1 >= ceild(29*t3-6,38)) && (t2 == 0) && (t3 <= floord(2*n-32,29))) {
          if (t3%2 == 0) {
            for (t5=1;t5<=min(min(24,floord(29*t3+26,2)),floord(-29*t3+2*n-30,2));t5++) {
              lbv=max(ceild(29*t3+2*t5+30,2),29*t3);
              ubv=min(n,29*t3+28);
#pragma ivdep
#pragma vector always
              for (t6=lbv;t6<=ubv;t6++) {
                S0(((29*t3+30)/2), t5, t6, -((29*t3+30)/2) + t6 + 1);;
                S0(((29*t3+30)/2), t5, t6, ((29*t3+30)/2) + t5 - 1);;
              }
            }
          }
        }
        for (t4=max(ceild(-25*t2+29*t3+31,2),19*t1);t4<=min(min(min(min(n-1,19*t1+18),-25*t2+n),29*t3+27),-25*t2+29*t3+28);t4++) {
          for (t5=max(1,25*t2);t5<=min(min(25*t2+24,-t4+n),29*t3-t4+28);t5++) {
            lbv=max(29*t3,t4+t5);
            ubv=min(n,29*t3+28);
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              S0(t4, t5, t6, -t4 + t6 + 1);;
              S0(t4, t5, t6, t4 + t5 - 1);;
            }
          }
        }
      }
    }
  }
}
  double execution_time = omp_get_wtime() - start;
  printf("normal: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(d, n, 0);
  deallocateMatrix(c, n + 1);
  deallocateMatrix(d, n + 1);
}

void computeDYN1(int** matrix, int n) {
  int** c = allocateMatrix(n + 1);
  int** d = allocateMatrix(n + 1);
  int i, j;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      c[i][j] = matrix[i][j];
  double start = omp_get_wtime();
  for (int c0 = 2; c0 < n; c0 += 1)
#pragma omp parallel for private(c1, c2, c0)
    for (int c1 = 1; c1 <= n - c0; c1 += 1)
      for (int c2 = c0 + c1; c2 <= min(n, 2 * c0 + c1 - 2); c2 += 1) {
        if (2 * c0 + c1 >= c2 + 3)
          S0(c0, c1, c2, -c0 + c2 + 1);
        S0(c0, c1, c2, c0 + c1 - 1);
      }
  double execution_time = omp_get_wtime() - start;
  printf("parallel: %lf\n", execution_time);
  write_results(n, execution_time);
  printMatrix(d, n, 1);
  deallocateMatrix(c, n + 1);
  deallocateMatrix(d, n + 1);
}

void computeDYN2(int** matrix, int n) {
  int** c = allocateMatrix(n + 1);
  int** d = allocateMatrix(n + 1);
  int i, j;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      c[i][j] = matrix[i][j];
  double start = omp_get_wtime();
  int tile_size = 2;
  for (int c0 = 0; c0 < floord(n, tile_size); c0 += 1)
    for (int c1 = 0; c1 < min(-c0 + n / tile_size, (n + 1) / tile_size - 1); c1 += 1)
      for (int c2 = max(c0 + c1, c1 + 1); c2 <= min(tile_size * c0 + c1 + 1, (n + 1) / tile_size - 1); c2 += 1)
        for (int c3 = max(tile_size * c0 + 1, -c1 + c2 + 1); c3 <= min(tile_size * c0 + 2, -tile_size * c1 + tile_size * c2 + 1); c3 += 1)
#pragma omp parallel for 
          for (int c4 = tile_size * c1 + 1; c4 <= min(min(tile_size * c1 + 2, n - c3), tile_size * c2 - c3 + 2); c4 += 1)
            for (int c5 = max(tile_size * c2 + 1, c3 + c4); c5 <= min(min(n, tile_size * c2 + 2), tile_size * c3 + c4 - 2); c5 += 1) {
              if (tile_size * c3 + c4 >= c5 + 3)
                S0(c3, c4, c5, -c3 + c5 + 1);
              S0(c3, c4, c5, c3 + c4 - 1);
            }
  double execution_time = omp_get_wtime() - start;
  printf("tiles: %lf\n", execution_time);
  write_results(n, execution_time, '\n');
  printMatrix(d, n, 2);
  deallocateMatrix(c, n + 1);
  deallocateMatrix(d, n + 1);
}

void printMatrix(int** matrix, int N, int fileno) {
  char filename[10];
  sprintf_s(filename, "nontiled%d", fileno);
  FILE* f;
  fopen_s(&f, filename, "wt");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++)
      fprintf(f, "%d ", matrix[i][j]);
    fprintf(f, "\n");
  }
  fclose(f);
}

int** allocateMatrix(int N) {
  int** t = (int**)malloc(sizeof(int*) * N);
  for (int i = 0; i < N; i++) {
    t[i] = (int*)malloc(sizeof(int) * N);
  }
  return t;
}

void deallocateMatrix(int **t, int N) {
  for (int i = 0; i < N; i++) {
    free(t[i]);
  }
  free(t);
}

void write_results(int n, double execution_time, char end_char)
{
  FILE* f;
  fopen_s(&f, "results.txt", "at");
  fprintf(f, "%d:%lf%c", n, execution_time, end_char);
  fclose(f);
}
void write_results(int n, double execution_time)
{
  write_results(n, execution_time, ';');
}

int main(void) {
  const int ZMAX = 120;
  int** graph = allocateMatrix(ZMAX);
  int g[4][4] = { {1, 1, 0, 1}, {0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1} };
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      graph[i][j] = g[i][j];
  for (int i = 0; i < ZMAX; i++)
    graph[i][i] = 1;
  int N = 110;
  while (N < ZMAX)
  {

  //printMatrix(graph, 6, 9);
    computeDYN0(graph, N);
    computeDYN1(graph, N);
    computeDYN2(graph, N);
    N += 10;
  }
  deallocateMatrix(graph, ZMAX);
  return 0;
}

