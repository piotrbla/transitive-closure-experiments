#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <omp.h>
#include <math.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))

double ** Q; double ** Q1; double ** Qbp; double ** Qbp1;
double ** Pbp; double ** Pu; double ** M; int CHECK_VALID = 1;

int Ebp = 0; // Energy weight of base pair  -2, -1, 0, 1, 2
int RT = 1; // 'Normalized' temperature 1,2,3,4,5
float ERT;
int l = 1; //minimum loop length 0-5
int delta = 1;  // Base pair weighting  1-5

unsigned char * RNA;  //only ACGU
int N; int DIM;

#include "mem.h"

int paired(int i, int j) {
   char nt1 = RNA[i];
   char nt2 = RNA[j];
         if ((nt1 == 'A' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'A') ||
             (nt1 == 'G' && nt2 == 'C') || (nt1 == 'C' && nt2 == 'G') ||
             (nt1 == 'G' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'G')){

            return 1;}
         else
            return 0;
}

int main(int argc, char *argv[]){
    int num_proc=1;
    int ll,p,q;
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    int t1, t2, t3, t4, t5, t6;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;

    ERT = exp((float)-Ebp/(float)RT);
    //printf("ERT: %.f %.18f \n", ERT, ERT);
 
    srand(time(NULL));
    if(argc > 1)
        num_proc = atoi(argv[1]);
    int kind=1;
    N = 8;
    DIM = 12;
    if(argc > 2)
        N = atoi(argv[2]);
    DIM = N+10;
    if(argc > 3)
        kind = atoi(argv[3]);
    omp_set_num_threads(num_proc);
    //printf(" -exp(Ebp/RT) = %5.3f\n", ERT);
    RNA =  (unsigned char*) malloc(DIM * sizeof(unsigned char*));  //read from FASTA file
    rand_seq(RNA, N);
//    printf("Sequence: ");
//    for(i=0; i<N; i++)
//       printf("%c", RNA[i]);
//    printf("\n\n");
    Q = memd();Q1 = memd();
    Qbp = memd(); Qbp1 = memd();
    Pbp = memd(); Pu = memd(); M = memd();

    rna_array_init(Q, 0.4, 0.4);
    rna_array_init(Q1, 0.4, 0.4);
    rna_array_init(Qbp, 0.5, 0.5);
    rna_array_init(Qbp1, 0.5, 0.5);
    rna_array_init(Pbp, 0, 0);
    rna_array_init(Pu, 0, 0);
    rna_array_init(M, 0, 0);
//for(i=0; i<N; i++)
  //for(j=0; j<N; j++)
    // printf("%.6f %.6f -- %d %d\n", Q[i][j], Q1[i][j], i, j);

    double start = omp_get_wtime();
    //  compute the partition functions Q and Qbp
    if (kind == 1 || CHECK_VALID)
    {
      int assignment_counter = 0;
      int add_counter = 0;
#pragma scop
      if (N >= 1 && l >= 0 && l <= 5)
        for (int i = N - 1; i >= 0; i--)
        {
          for (int j = i + 1; j < N; j++)
          {
            Q1[i][j] = Q1[i][j - 1];
            assignment_counter++;
            for (int k = 0; k < j - i - 1; k++) //zamienione l na 1
            {
              Qbp1[k + i][j] = Q1[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
              Q1[i][j] += Q1[i][k + i] * Qbp1[k + i][j];
              add_counter++;
              // printf("%.f\n", Q1[i][j]);
            }
            //write_result_matrix(Q1, N);
          }
        }
#pragma endscop
      //printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }




    if (kind == 11 )
    {
      int assignment_counter = 0;
      int add_counter = 0;
#pragma scop
      if (N >= 2 && l >= 0 && l <= 5)
        Q[N - 2][N - 1] = Q[N - 2][N - 2];
      for (int i = N - 1; i >= 0; i--)
      {
        for (int j = i + 1; j < N; j++)
        {
          for (int k = 0; k < j - i - 1; k++)
          {
            Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
            add_counter++;
            // printf("%.f\n", Q1[i][j]);
          }
          if (j < N - 1)
          {
            Q[i][j + 1] = Q[i][j];
            assignment_counter++;
          }
          //write_result_matrix(Q, N);
        }
      }
#pragma endscop
      //printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }
    if (kind == 2) // pluto
    {
        /* Start of CLooG code */
        if ((N >= 2) && (l >= 0) && (l <= 5)) {
          for (t1=1;t1<=N-1;t1++) {
            for (t2=0;t2<=floord(t1-1,16);t2++) {
              for (t3=0;t3<=t2;t3++) {
                if ((t1 >= l+1) && (t2 == 0) && (t3 == 0)) {
                  Qbp[0 +0][t1] = Q[0 +0 +1][t1-1] * ERT * paired(0 +0,t1-1);;
                  Q[0][t1] = Q[0][t1-1];;
                  Q[0][t1] += Q[0][0 +0] * Qbp[0 +0][t1];;
                }
                if (t3 == 0) {
                  for (t4=max(1,16*t2);t4<=min(16*t2+15,t1-l-1);t4++) {
                    Qbp[0 +t4][t1] = Q[0 +t4+1][t1-1] * ERT * paired(0 +t4,t1-1);;
                    Q[t4][t1] = Q[t4][t1-1];;
                    Q[t4][t1] += Q[t4][0 +t4] * Qbp[0 +t4][t1];;
                    for (t5=1;t5<=min(15,t4);t5++) {
                      Qbp[t5+(t4-t5)][t1] = Q[t5+(t4-t5)+1][t1-1] * ERT * paired(t5+(t4-t5),t1-1);;
                      Q[(t4-t5)][t1] += Q[(t4-t5)][t5+(t4-t5)] * Qbp[t5+(t4-t5)][t1];;
                    }
                  }
                }
                if (t3 == 0) {
                  for (t4=max(16*t2,t1-l);t4<=min(t1-1,16*t2+15);t4++) {
                    Q[t4][t1] = Q[t4][t1-1];;
                  }
                }
                if (t3 >= 1) {
                  for (t4=16*t2;t4<=min(16*t2+15,t1-l-1);t4++) {
                    for (t5=16*t3;t5<=min(t4,16*t3+15);t5++) {
                      Qbp[t5+(t4-t5)][t1] = Q[t5+(t4-t5)+1][t1-1] * ERT * paired(t5+(t4-t5),t1-1);;
                      Q[(t4-t5)][t1] += Q[(t4-t5)][t5+(t4-t5)] * Qbp[t5+(t4-t5)][t1];;
                    }
                  }
                }
              }
            }
          }
        }
        /* End of CLooG code */
    }
    if(kind==3) // tile corr
{
printf("traco\n");
         if (N >= 10 && l >= 0 && l <= 5)
  for( c1 = 1; c1 < N + (N - 2) / 16; c1 += 1)
    #pragma omp parallel for schedule(dynamic, 1) shared(c1) private(c3,c4,c5,c9,c11)
    for( c3 = max(0, -N + c1 + 1); c3 <= (c1 - 1) / 17; c3 += 1)
      for( c4 = 0; c4 <= 1; c4 += 1) {
        if (c4 == 1) {
          for( c5 = 0; c5 <= c3; c5 += 1)
            for( c9 = N - c1 + 17 * c3; c9 <= min(N - 1, N - c1 + 17 * c3 + 15); c9 += 1) {
              if (c5 == c3 && c1 + c9 >= N + 17 * c3 + 1)
                Q[(N-c1+c3-1)][c9] = Q[(N-c1+c3-1)][c9-1];
              if (c5 == c3 && c1 + c9 >= N + 17 * c3 + 1)
                for( c11 = 0; c11 < 16 * c3; c11 += 1)
                  Q[(N-c1+c3-1)][c9] += Q[(N-c1+c3-1)][c11+(N-c1+c3-1)] * Qbp[c11+(N-c1+c3-1)][c9];
              for( c11 = 16 * c5; c11 <= min(16 * c5 + 15, -N + c1 - c3 + c9); c11 += 1) {
                Qbp[c11+(N-c1+c3-1)][c9] = Q[c11+(N-c1+c3-1)+1][c9-1] * ERT * paired(c11+(N-c1+c3-1),c9-1);
                if (c5 == c3) {
                  Q[(N-c1+c3-1)][c9] += Q[(N-c1+c3-1)][c11+(N-c1+c3-1)] * Qbp[c11+(N-c1+c3-1)][c9];
                } else if (c1 + c9 == N + 17 * c3) {
                  Q[(N-c1+c3-1)][(N-c1+17*c3)] += Q[(N-c1+c3-1)][c11+(N-c1+c3-1)] * Qbp[c11+(N-c1+c3-1)][(N-c1+17*c3)];
                }
              }
            }
        } else {
          Q[(N-c1+c3-1)][(N-c1+17*c3)] = Q[(N-c1+c3-1)][(N-c1+17*c3)-1];
        }
      }



}
    if (kind==4)
{
  printf("mccIfk==0\n");
  {
#pragma scop
    for (int i = N - 1; i >= 0; i--)
    {
      for (int j = i + 1; j < N; j++)
      {
        for (int k = 0; k < j - i - l; k++)
        {
          if (k == 0)
            Q[i][j] = Q[i][j - 1];
          Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
        }
      }
    }
#pragma endscop
  }
}
    if (kind == 5)
    {
      printf("mccNoIfMod\n"); // pamietac o zmianie Q1 na Q
      {
        //#pragma scop
        //     for(int i=N-1; i>=0; i--){
        //         for(int j=i+1; j<N; j++){
        // //printf("%.f\n", Q1[i][j]);
        //         for(int k=0; k<j-i-l; k++){
        //             if (k==0)
        //                 Q[i][j] = Q[i][j-1];
        //             Q[i][j] +=  Q[i][k+i] * Q[k+i+1][j-1] * ERT * paired(k+i,j-1);
        // //             printf("%.f\n", Q1[i][j]);
        //         }

        //         }
        //     }

        for (int c0 = 1; c0 < N; c0 += 1)
        {
          for (int c2 = 0; c2 < N - c0; c2 += 1)
            Q[c2][c0 + c2] = Q[c2][c0 + c2 - 1];
          // S_0(c2, c0 + c2, 0);
          for (int c1 = (c0 + 1) / 2; c1 < c0; c1 += 1)
            for (int c2 = 0; c2 < N - c0; c2 += 1)
            {
              if (2 * c1 >= c0 + 1)
                Q[c2][c0 + c2] += Q[c2][c0 - c1 - 1 + c2] * Q[c0 - c1 - 1 + c2 + 1][c0 + c2 - 1] * ERT * paired(c0 - c1 - 1 + c2, c0 + c2 - 1);
              // S_1(c2, c0 + c2, c0 - c1 - 1);
              Q[c2][c0 + c2] += Q[c2][c1 - 1 + c2] * Q[c1 - 1 + c2 + 1][c0 + c2 - 1] * ERT * paired(c1 - 1 + c2, c0 + c2 - 1);
              // S_1(c2, c0 + c2, c1 - 1);
            }
        }
        ;
        //#pragma endscop
      }
    }
    if (kind == 51)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod51\n"); // pamietac o zmianie Q1 na Q
      {
        for (int c2 = 0; c2 < N - 1; c2 += 1)
        {
          Q[c2][c2+1] = Q[c2][c2];
          assignment_counter++;
        }
        for (int c0 = 1; c0 < N; c0 += 1)
        {
          for (int c1 = (c0 + 1) / 2; c1 < c0; c1 += 1)
            for (int c2 = 0; c2 < N - c0; c2 += 1)
            {
              if (2 * c1 >= c0 + 1)
              {
                Q[c2][c0 + c2] += Q[c2][c0 - c1 - 1 + c2] * Q[c0 - c1 - 1 + c2 + 1][c0 + c2 - 1] * ERT * paired(c0 - c1 - 1 + c2, c0 + c2 - 1);
                add_counter++;
              }
              Q[c2][c0 + c2] += Q[c2][c1 - 1 + c2] * Q[c1 - 1 + c2 + 1][c0 + c2 - 1] * ERT * paired(c1 - 1 + c2, c0 + c2 - 1);
              add_counter++;
            }
            for (int c2 = 0; c2 < N - c0 - 1; c2 += 1)
            {
              Q[c2][c0 + 1 + c2] = Q[c2][c0 + c2];
              assignment_counter++;
            }
        }
        printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
        //#pragma endscop
      }
    }
    if (kind==52)
{
      int assignment_counter = 0;
      int add_counter = 0;
if (N >= 2) {
  for (int w0 = 0; w0 <= (3 * N - 3) / 32; w0 += 1) {
    #pragma omp parallel for
    for (int h0 = max(0, w0 - (N + 31) / 32 + 1); h0 <= min((N - 1) / 16, w0 - (w0 + 1) / 3); h0 += 1) {
      if (N >= 3) {
        for (int t0 = max(4 * h0, 4 * w0 - 2 * h0); t0 <= min(min(4 * w0 - 2 * h0 + 5, (N - 1) / 4), 2 * h0 + (N - 2) / 8 + 2); t0 += 1) {
          if (h0 == 0 && 4 * w0 + 3 >= t0) {
            for (int i3 = 8 * t0; i3 <= min(N - 2, 8 * t0 + 7); i3 += 1) {
              Q[i3][i3 + 1] = Q[i3][i3];
            }
          }
          for (int i2 = max(max(max(max(max(1, 16 * h0), -N + 32 * w0 - 32 * h0 + 3), -32 * w0 + 32 * h0 + 8 * t0 - 31), -N + 8 * t0 + 1), -N + 4 * t0 + N / 2 + 2); i2 <= min(min(min(N - 1, 16 * h0 + 15), 4 * t0 + 3), -32 * w0 + 32 * h0 + 8 * t0 + 7); i2 += 1) {
            for (int i3 = max(max(32 * w0 - 32 * h0, i2), 8 * t0 - i2); i3 <= min(min(min(N - 1, 32 * w0 - 32 * h0 + 31), N + i2 - 3), 8 * t0 - i2 + 7); i3 += 1) {
              for (int i5 = (i2 + 1) / 2; i5 < i2; i5 += 1) {
                if (2 * i5 >= i2 + 1) {
                  Q[-i2 + i3][i3] += (((Q[-i2 + i3][i3 - i5 - 1] * Q[i3 - i5][i3 - 1]) * ERT) * paired((i3 - i5 - 1), (i3 - 1)));
                }
                Q[-i2 + i3][i3] += (((Q[-i2 + i3][-i2 + i3 + i5 - 1] * Q[-i2 + i3 + i5][i3 - 1]) * ERT) * paired((-i2 + i3 + i5 - 1), (i3 - 1)));
              }
              if (N >= i3 + 2) {
                Q[-i2 + i3][i3 + 1] = Q[-i2 + i3][i3];
              }
            }
          }
        }
      } else {
        Q[0][1] = Q[0][0];
      }
    }
  }
}
for (int i2 = 0; i2 < N - 1; i2 += 1) {
  assignment_counter++;
}
for (int i2 = 1; i2 < N - 1; i2 += 1) {
  for (int i3 = 0; i3 < N - i2 - 1; i3 += 1) {
    assignment_counter++;
  }
}
for (int i2 = 2; i2 < N; i2 += 1) {
  for (int i3 = (i2 + 1) / 2; i3 < i2; i3 += 1) {
    for (int i5 = 0; i5 < N - i2; i5 += 1) {
      if (2 * i3 >= i2 + 1) {
        add_counter++;
      }
      add_counter++;
    }
  }
}
        printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);

}
//Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);;
//Q[i][j] = Q[i][j - 1];
    if (kind == 53)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod53\n"); // pamietac o zmianie Q1 na Q
      {
        if (N >= 2 && l >= 0 && l <= 5)
        {
          Q[N - 2][N - 1] = Q[N - 2][N - 2];//niepotrzebne w ogóle
        }
        for (int i = N - 1; i >= 0; i--)
        {
          if (i > 0)
          {
            Q[i - 1][N - 1] = Q[i - 1][N - 2];
            assignment_counter++;
          }
        }

        for (int c0 = -N + 3; c0 <= min(0, -N + l + 1); c0 += 1)
          for (int c1 = -c0 + 1; c1 < N - 1; c1 += 1)
          {
            Q[-c0][c1] = Q[-c0][c1 - 1];
            //S_1(-c0, c1);
            assignment_counter++;
          }
        for (int c0 = max(-N + 2, -N + l + 2); c0 <= 0; c0 += 1)
          for (int c1 = -c0 + 1; c1 < N; c1 += 1)
          {
            for (int c3 = 0; c3 < -l + c0 + c1; c3 += 1)
            {
              //Q[i][ j] += Q[  i] [k +  i] *  Q[k +  i + 1]   [j - 1] * ERT * paired(k +    i, j - 1);
              Q[-c0][c1] += Q[-c0][c3 + -c0] * Q[c3 + -c0 + 1][c1 - 1] * ERT * paired(c3 + -c0, c1 - 1);
              //S_0(-c0, c1, c3);
              add_counter++;
            }
            if (N >= c1 + 2)
            {
              Q[-c0][c1] = Q[-c0][c1 - 1];
              //S_1(-c0, c1);
              assignment_counter++;
            }
          }

        printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
        //#pragma endscop
      }
    }
    if (kind == 54)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod54\n"); // pamietac o zmianie Q1 na Q
      {
        for (int i = N - 1; i >= 0; i--)
        {
          if (i > 0)
          {
            Q[i - 1][N - 1] = Q[i - 1][N - 2];
            assignment_counter++;
          }
        }

        for (int c0 = 1; c0 < N - 1; c0 += 1)
          for (int c1 = 0; c1 < N - c0 - 1; c1 += 1)
            for (int c2 = c0 + c1 + 1; c2 <= min(N - 1, 2 * c0 + c1); c2 += 1)
            {
              if (2 * c0 + c1 >= c2 + 1)
              {
                Q[c1][c2] += Q[c1][-c0 - c1 + c2 - 1 + c1] * Q[-c0 - c1 + c2 - 1 + c1 + 1][c2 - 1] * ERT * paired(-c0 - c1 + c2 - 1 + c1, c2 - 1);
                //S_0(c1, c2, -c0 - c1 + c2 - 1);
                //Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
                add_counter++;
              }
              //Q[c1][c2] += Q[c1][c0 - 1 + c1] * Q[c0 - 1 + c1 + 1][c2 - 1] * ERT * paired(c0 - 1 + c1, c2 - 1);
              //S_0(c1, c2, c0 - 1);
              add_counter++;
              if (N >= c0 + c1 + 3 && c2 == c0 + c1 + 1)
              {
                Q[c1][c0 + c1 + 1] = Q[c1][c0 + c1 + 1 - 1];
                //S_1(c1, c0 + c1 + 1);
                assignment_counter++;
              }
            }
      printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
      //#pragma endscop
    }
}
//Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);;
//Q[i][j] = Q[i][j - 1];
if (kind == 55)
{
  int assignment_counter = 0;
  int add_counter = 0;
  printf("mccNoIfMod55\n"); // pamietac o zmianie Q1 na Q
  {
    for (int i = N - 1; i >= 0; i--)
    {
      if (i > 0)
      {
        Q[i - 1][N - 1] = Q[i - 1][N - 2];
        assignment_counter++;
      }
    }

    
    for (int c0 = 1; c0 < N - 1; c0 += 1)
      for (int c1 = -N + c0 + 2; c1 <= 0; c1 += 1)
        for (int c2 = c0 - c1 + 1; c2 <= min(N - 1, 2 * c0 - c1); c2 += 1)
        {
          if (2 * c0 >= c1 + c2 + 1)
          {
            //S_0(-c1, c2, -c0 + c1 + c2 - 1);
            Q[-c1][c2] += Q[-c1][-c0 + c1 + c2 - 1 + -c1] * Q[-c0 + c1 + c2 - 1 + -c1 + 1][c2 - 1] * ERT * paired(-c0 + c1 + c2 - 1 + -c1, c2 - 1);;
            add_counter++;
          }
          Q[-c1][c2] += Q[-c1][c0 - 1 + -c1] * Q[c0 - 1 + -c1 + 1][c2 - 1] * ERT * paired(c0 - 1 + -c1, c2 - 1);;
          //S_0(-c1, c2, c0 - 1);
//Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);;
          add_counter++;
          if (N + c1 >= c0 + 3 && c1 + c2 == c0 + 1)
          {
            //S_1(-c1, c0 - c1 + 1);
            Q[-c1][c0 - c1 + 1] = Q[-c1][c0 - c1 + 1 - 1];
            assignment_counter++;
          }
        }
    printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
    //#pragma endscop
  }
}
if (kind == 56)
{
  int assignment_counter = 0;
  int add_counter = 0;
  printf("mccNoIfMod55\n"); // pamietac o zmianie Q1 na Q
  {
    for (int i = N - 1; i >= 0; i--)
    {
      if (i > 0)
      {
        Q[i - 1][N - 1] = Q[i - 1][N - 2];
        assignment_counter++;
      }
    }

    for (int c0 = -N + 3; c0 <= 0; c0 += 1)
      for (int c1 = -c0 + 1; c1 < N; c1 += 1)
      {
        for (int c3 = 0; c3 < c0 + c1 - 1; c3 += 1)
        {
          Q[-c0][c1] += Q[-c0][c3 + -c0] * Q[c3 + -c0 + 1][c1 - 1] * ERT * paired(c3 + -c0, c1 - 1);;
          // S_0(-c0, c1, c3);
          add_counter++;
        }
        if (N >= c1 + 2)
        {
          assignment_counter++;
          Q[-c0][c1] = Q[-c0][c1 - 1];
          // S_1(-c0, c1);
        }
      }
    printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
    //#pragma endscop
  }
}
    if (kind == 57 )
    {
      int assignment_counter = 0;
      int add_counter = 0;
  printf("mccNoIfMod57\n"); // pamietac o zmianie Q1 na Q
    for (int i = N - 1; i >= 0; i--)
    {
      if (i > 0)
      {
        Q[i - 1][N - 1] = Q[i - 1][N - 2];
        assignment_counter++;
      }
    }

#pragma scop
        for (int i = N - 1; i >= 0; i--)
        {
          for (int j = i + 1; j < N; j++)
          {
            for (int k = 0; k < j - i - 1; k++)
            {
              Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
              add_counter++;
              // printf("%.f\n", Q1[i][j]);
            }
            if (j<N-1)
            {
              Q[i][j+1] = Q[i][j];
              assignment_counter++;
            }  
          }

        }
#pragma endscop
      printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    if (kind == 571)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod571\n"); // pamietac o zmianie Q1 na Q
      for (int i = N - 1; i >= 0; i--)
      {
        if (i > 0)
        {
          Q[i - 1][N - 1] = Q[i - 1][N - 2];
          assignment_counter++;
        }
      }

#pragma scop
      for (int i = N - 1; i >= 0; i--)
      {
        for (int j = i + 1; j < N; j++)
        {
          for (int k = 0; k <= j - i - 1; k++)
          {
            if (k < j - i - 1)
            {
              add_counter++;
              Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
            } // S_0

            if (j < N - 1 && k == j - i - 2 && k >= 0)
            {
              {
                Q[i][j + 1] = Q[i][j];
                assignment_counter++;
              } // S_1
            }
          }
        }

#pragma endscop
      }
        printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    if (kind == 572)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod572\n"); // pamietac o zmianie Q1 na Q
      for (int i = N - 1; i >= 0; i--)
      {
        if (i > 0)
        {
          Q[i - 1][N - 1] = Q[i - 1][N - 2];
          assignment_counter++;
        }
      }

#pragma scop
      for (int c0 = 1; c0 < N; c0 += 1)
        for (int c1 = 0; c1 < N - c0; c1 += 1)
          for (int c2 = c0 + c1; c2 <= min(N - 1, 2 * c0 + c1); c2 += 1)
          {
            if (c2 >= c0 + c1 + 1 && 2 * c0 + c1 >= c2 + 1)
            {
              if (-c0 - c1 + c2 - 1 < c2 - c1 - 1)
              {
                add_counter++;
                Q[c1][c2] += Q[c1][-c0 - c1 + c2 - 1 + c1] * Q[-c0 - c1 + c2 - 1 + c1 + 1][c2 - 1] * ERT * paired(-c0 - c1 + c2 - 1 + c1, c2 - 1);
              } // S_0

              if (c2 < N - 1 && -c0 - c1 + c2 - 1 == c2 - c1 - 1 && -c0 - c1 + c2 - 1 >= 0)
              {
                {
                  Q[c1][c2 + 1] = Q[c1][c2];
                  assignment_counter++;
                } // S_1
              }
              // S(c1, c2, -c0 - c1 + c2 - 1);
            }
            if (c0 - 1 < c2 - c1 - 1)
            {
              add_counter++;
              Q[c1][c2] += Q[c1][c0 - 1 + c1] * Q[c0 - 1 + c1 + 1][c2 - 1] * ERT * paired(c0 - 1 + c1, c2 - 1);
            } // S_0

            if (c2 < N - 1 && c0 - 1 == c2 - c1 - 1 && c0 - 1 >= 0)
            {
              {
                Q[c1][c2 + 1] = Q[c1][c2];
                assignment_counter++;
              } // S_1
            }
            // S(c1, c2, c0 - 1);
          }

#pragma endscop
      printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    if (kind == 573)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod573\n"); // pamietac o zmianie Q1 na Q
      for (int i = N - 1; i >= 0; i--)
      {
        if (i > 0)
        {
          Q[i - 1][N - 1] = Q[i - 1][N - 2];
          assignment_counter++;
        }
      }

#pragma scop
      for (int i = N - 1; i >= 0; i--)
      {
        for (int j = i + 1; j < N; j++)
        {
          for (int k = 0; k <= j - i - 1; k++)
          {
            if (k < j - i - 1)
            {
              add_counter++;
              Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
            } // S_0

            if (j < N - 1 && k == j - i - 1 && k >= 0)
            // if (j < N - 1 && k == j - i - 2 && k >= 0)
            {
              {
                Q[i][j + 1] = Q[i][j];
                assignment_counter++;
              } // S_1
            }
          }
        }

#pragma endscop
      }
        printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    // Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);;
    // Q[i][j] = Q[i][j - 1];
    // Q[i][j+1] = Q[i][j]; // sprawdzic, to pomodyfikacji
    
    if (kind == 58)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod58\n"); // pamietac o zmianie Q1 na Q
      for (int i = N - 1; i >= 0; i--)
      {
        if (i > 0)
        {
          Q[i - 1][N - 1] = Q[i - 1][N - 2];
          assignment_counter++;
        }
      }

#pragma scop
      for (int c0 = -N + 3; c0 <= 0; c0 += 1)
        for (int c1 = -c0 + 1; c1 < N; c1 += 1)
        {
          for (int c3 = 0; c3 < c0 + c1 - 1; c3 += 1)
          {
            Q[-c0][c1] += Q[-c0][c3 + -c0] * Q[c3 + -c0 + 1][c1 - 1] * ERT * paired(c3 + -c0, c1 - 1);;
            //S_0(-c0, c1, c3);
            add_counter++;
          }
          if (N >= c1 + 2)
          {
            //Q[-c0][c1] = Q[-c0][c1 - 1];
            Q[-c0][c1+1] = Q[-c0][c1];
            //S_1(-c0, c1);
            assignment_counter++;
          }
        }

#pragma endscop
      printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }
    if (kind == 59)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod59\n"); // pamietac o zmianie Q1 na Q
      for (int i = N - 1; i >= 0; i--)
      {
        if (i > 0)
        {
          Q[i - 1][N - 1] = Q[i - 1][N - 2];
          assignment_counter++;
        }
      }

#pragma scop
      for (int c0 = 1; c0 < N - 1; c0 += 1)

        for (int c1 = 0; c1 < N - c0 - 1; c1 += 1)
          for (int c2 = c0 + c1 + 1; c2 <= min(N - 1, 2 * c0 + c1); c2 += 1)
          {
            if (2 * c0 + c1 >= c2 + 1)
            { 
              Q[c1][c2] += Q[c1][-c0 - c1 + c2 - 1 + c1] * Q[-c0 - c1 + c2 - 1 + c1 + 1][c2 - 1] * ERT * paired(-c0 - c1 + c2 - 1 + c1, c2 - 1);;
              // S_0(c1, c2, -c0 - c1 + c2 - 1);
              add_counter++;
            }
            Q[c1][c2] += Q[c1][c0 - 1 + c1] * Q[c0 - 1 + c1 + 1][c2 - 1] * ERT * paired(c0 - 1 + c1, c2 - 1);
            // S_0(c1, c2, c0 - 1);
            add_counter++;
            if (N >= c0 + c1 + 3 && c2 == c0 + c1 + 1)
            { 
              //Q[c1][c0 + c1 + 1] = Q[c1][c0 + c1 + 1 - 1];
              Q[c1][c0 + c1 + 1+1] = Q[c1][c0 + c1 + 1];
              // S_1(c1, c0 + c1 + 1);
              assignment_counter++;
            }
          }

#pragma endscop
      printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }
    if (kind == 591)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod591\n"); // pamietac o zmianie Q1 na Q
      for (int i =0; i <= N-3; i++)
      {
        Q[i][i+2] = Q[i][i+1];
        assignment_counter++;
      }
#pragma scop

      for (int i = N - 1; i >= 0; i--)
      {
        for (int j = i + 1; j < N; j++)
        {
          for (int k = 0; k < j - i - 1; k++)
          {
            Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1); // S_0
            add_counter++;
            if (j < N - 1 && k == j - i - 2 && k >= 0)
            {
              Q[i][j + 1] = Q[i][j];
              assignment_counter++;
            } // S_1
          }
        }
      }

      

      

#pragma endscop
      printf("assignemento_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    if (kind == 592)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod592\n"); // pamietac o zmianie Q1 na Q
      for (int i = 0; i <= N - 3; i++)
      {
        Q[i][i + 2] = Q[i][i + 1];
        assignment_counter++;
      }
#pragma scop

      for (int c0 = -N + 3; c0 <= 0; c0 += 1)
        for (int c1 = -c0 + 2; c1 < N; c1 += 1)
        {
          for (int c2 = 0; c2 < c0 + c1 - 1; c2 += 1)
          {
            add_counter++;
            Q[-c0][c1] += Q[-c0][c2 + -c0] * Q[c2 + -c0 + 1][c1 - 1] * ERT * paired(c2 + -c0, c1 - 1);
            // S_0(-c0, c1, c2);
          }
          if (N >= c1 + 1)
          {
            assignment_counter++;
            Q[-c0][c1 + 1] = Q[-c0][c1];
            // S_1(-c0, c1, c0 + c1 - 2);
          }
        }

#pragma endscop
      printf("assignemento_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    if (kind == 593)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod593\n"); // pamietac o zmianie Q1 na Q
      for (int i = 0; i <= N - 3; i++)
      {
        Q[i][i + 2] = Q[i][i + 1];
        assignment_counter++;
      }
#pragma scop

      for (int c0 = 1; c0 < N - 1; c0 += 1)
        for (int c1 = -N + c0 + 2; c1 <= 0; c1 += 1)
          for (int c2 = c0 - c1 + 1; c2 <= min(N - 1, 2 * c0 - c1); c2 += 1)
          {
            if (2 * c0 >= c1 + c2 + 1)
            {
              Q[-c1][c2] += Q[-c1][-c0 + c1 + c2 - 1 + -c1] * Q[-c0 + c1 + c2 - 1 + -c1 + 1][c2 - 1] * ERT * paired(-c0 + c1 + c2 - 1 + -c1, c2 - 1); // S_0
              add_counter++;
              if (c2 < N - 1 &&-c0 + c1 + c2 - 1 == c2 - (-c1) - 2 && -c0 + c1 + c2 - 1 >= 0)
              {
                {
                  Q[-c1][c2 + 1] = Q[-c1][c2];
                  assignment_counter++;
                } // S_1
              }
              // S(-c1, c2, -c0 + c1 + c2 - 1);
            }
            {
              Q[-c1][c2] += Q[-c1][c0 - 1 + -c1] * Q[c0 - 1 + -c1 + 1][c2 - 1] * ERT * paired(c0 - 1 + -c1, c2 - 1); // S_0
              add_counter++;
              if (c2 < N - 1 && c0 - 1 == c2 - (-c1) - 2 && c0 - 1 >= 0)
              {
                {
                  Q[-c1][c2 + 1] = Q[-c1][c2];
                  assignment_counter++;
                } // S_1
              }
              // S(-c1, c2, c0 - 1);
            }
          }

#pragma endscop
      printf("assignemento_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    if (kind ==594)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod592\n"); // pamietac o zmianie Q1 na Q
      for (int i = 0; i <= N - 3; i++)
      {
        Q[i][i + 2] = Q[i][i + 1];
        assignment_counter++;
      }

      //   if (2 * c0 + c1 >= c2 + 2)
      //     S(c1, c2, -c0 - c1 + c2);
      //   if (c2 >= c0 + c1 + 1)
      //     S(c1, c2, c0 - 1);
      // }
      for (int c0 = 2; c0 < N; c0 += 1)
        for (int c1 = 0; c1 < N - c0; c1 += 1)
          for (int c2 = c0 + c1; c2 < min(N, 2 * c0 + c1); c2 += 1)
          {
            {
              if (2 * c0 + c1 >= c2 + 2)
              {
                add_counter++;
                Q[c1][c2] += Q[c1][-c0 - c1 + c2 + c1] * Q[-c0 - c1 + c2 + c1 + 1][c2 - 1] * ERT * paired(-c0 - c1 + c2 + c1, c2 - 1); // S_0
                if (c2 < N - 1 && (-c0 - c1 + c2) == (c2 - c1 - 2) && -c0 - c1 + c2 >= 0 || c2-c1==2)
                {
                  assignment_counter++;
                  Q[c1][c2 + 1] = Q[c1][c2]; // S_1
                }
                // S(c1, c2, -c0 - c1 + c2 - 1);
              }
              if (c2 >= c0 + c1 + 1)
              {
                add_counter++;
                Q[c1][c2] += Q[c1][c0 - 1 + c1] * Q[c0 - 1 + c1 + 1][c2 - 1] * ERT * paired(c0 - 1 + c1, c2 - 1); // S_0
                if (c2 < N - 1 && (c0 - 1) == (c2 - c1 - 2) && c0 - 1 >= 0 || c2-c1==2)
                {
                  assignment_counter++;
                  Q[c1][c2 + 1] = Q[c1][c2]; // S_1
                }

                // S(c1, c2, c0 - 1);
              }
            }
          }
          printf("assignemento_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }




    if (kind ==595)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod595\n"); // pamietac o zmianie Q1 na Q
      // for (int i = 0; i <= N - 3; i++)
      // {
      //   Q[i][i + 2] = Q[i][i + 1];
      //   assignment_counter++;
      // }
      for (int i = 0; i <= N - 1; i++)
      {
        Q[i][i + 1] = Q[i][i];
        assignment_counter++;
      }
      for (int c0 = 1; c0 < N - 1; c0 += 1)
        for (int c1 = 0; c1 < N - c0 - 1; c1 += 1)
          for (int c2 = c0 + c1 + 1; c2 <= min(N - 1, 2 * c0 + c1); c2 += 1)
          {
            if (2 * c0 + c1 >= c2 + 1)
              {
                add_counter++;
                Q[c1][c2] += Q[c1][-c0 - c1 + c2 + c1] * Q[-c0 - c1 + c2 + c1 + 1][c2 - 1] * ERT * paired(-c0 - c1 + c2 + c1, c2 - 1); // S_0
                if (c2 < N - 1 && (-c0 - c1 + c2) == (c2 - c1 - 2) && -c0 - c1 + c2 >= 0 || c2-c1==2)
                {
                  assignment_counter++;
                  Q[c1][c2 + 1] = Q[c1][c2]; // S_1
                }
                //S(c1, c2, -c0 - c1 + c2 - 1);
              }
                add_counter++;
                Q[c1][c2] += Q[c1][c0 - 1 + c1] * Q[c0 - 1 + c1 + 1][c2 - 1] * ERT * paired(c0 - 1 + c1, c2 - 1); // S_0
                if (c2 < N - 1 && (c0 - 1) == (c2 - c1 - 2) && c0 - 1 >= 0 || c2-c1==2)
                {
                  assignment_counter++;
                  Q[c1][c2 + 1] = Q[c1][c2]; // S_1
                }
                //S(c1, c2, c0 - 1);
          }
          printf("assignemento_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    if (kind == 596)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod596\n"); // pamietac o zmianie Q1 na Q
#pragma scop
      for (int i = N - 1; i >= 0; i--)
      {
        for (int j = i + 1; j <= i + 1; j++)
        {

          if (j < N - 1)
          {
            Q[i][j + 1] = Q[i][j];
            assignment_counter++;
          }
        }
      }

      for (int i = N - 1; i >= 0; i--)
      {
        for (int j = i + 2; j < N; j++)
        {
          for (int k = 0; k < j - i - 1; k++)
          {
            Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
            add_counter++;
          }
          if (j < N - 1)
          {
            Q[i][j + 1] = Q[i][j];
            assignment_counter++;
          }
        }
      }

#pragma endscop
      printf("assignemento_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    if (kind == 597)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mccNoIfMod597\n"); // pamietac o zmianie Q1 na Q
      for (int i = N - 1; i >= 0; i--)
      {
        if (i > 0)
        {
          Q[i - 1][N - 1] = Q[i - 1][N - 2];
          assignment_counter++;
        }
      } 

      for (int i = N - 1; i >= 0; i--)
      {
        for (int j = i + 1; j <= i + 1; j++)
        {

          if (j < N - 1)
          {
            Q[i][j + 1] = Q[i][j];
            assignment_counter++;
          }
        }
      }
#pragma scop
      for (int c0 = 1; c0 < N - 1; c0 += 1)
        for (int c1 = 0; c1 < N - c0 - 1; c1 += 1)
          for (int c2 = c0 + c1 + 1; c2 <= min(N - 1, 2 * c0 + c1); c2 += 1)
          {
            if (2 * c0 + c1 >= c2 + 1)
            {
              //S_0(c1, c2, -c0 - c1 + c2 - 1);
              add_counter++;
              Q[c1][c2] += Q[c1][-c0 - c1 + c2 - 1 + c1] * Q[-c0 - c1 + c2 - 1 + c1 + 1][c2 - 1] * ERT * paired(-c0 - c1 + c2 - 1 + c1, c2 - 1);
            }
            add_counter++;
            Q[c1][c2] += Q[c1][c0 - 1 + c1] * Q[c0 - 1 + c1 + 1][c2 - 1] * ERT * paired(c0 - 1 + c1, c2 - 1);
            //S_0(c1, c2, c0 - 1);
            if (c2 == c0 + c1 + 1)
            {
              assignment_counter++;
              Q[c1][c0 + c1 + 1 + 1] = Q[c1][c0 + c1 + 1];
              //S_1(c1, c0 + c1 + 1);
            }
          }
          //       Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
          //       Q[i][j + 1] = Q[i][j];

#pragma endscop
      printf("assignemento_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    if (kind == 598)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      //printf("mccNoIfMod598\n"); 
  
      for (int i = N - 1; i >= 0; i--)
      {
        for (int j = i + 1; j <= i + 1; j++)
        {
          if (j < N - 1)
          {
            Q[i][j + 1] = Q[i][j];
            assignment_counter++;
          }
          //write_result_matrix(Q, N);
        }
      }
#pragma scop
      for (int c0 = 2; c0 < N; c0 += 1)
        for (int c1 = 0; c1 < N - c0; c1 += 1)
        {
          for (int c2 = c0 + c1; c2 < min(N, 2 * c0 + c1 - 1); c2 += 1)
          {
            if (2 * c0 + c1 >= c2 + 3)
            {
              add_counter++;
              Q[c1][c2] += Q[c1][-c0 - c1 + c2 + c1] * Q[-c0 - c1 + c2 + c1 + 1][c2 - 1] * ERT * paired(-c0 - c1 + c2 + c1, c2 - 1);
              // S_0(c1, c2, -c0 - c1 + c2);
            }
            add_counter++;
            Q[c1][c2] += Q[c1][c0 - 2 + c1] * Q[c0 - 2 + c1 + 1][c2 - 1] * ERT * paired(c0 - 2 + c1, c2 - 1);
            // S_0(c1, c2, c0 - 2);
            if (c2 == c0 + c1)
            {
              Q[c1][c0 + c1 + 1] = Q[c1][c0 + c1];
              assignment_counter++;
              // S_1(c1, c0 + c1);
            }
            //       Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
            //       Q[i][j + 1] = Q[i][j];

          }
          //write_result_matrix(Q, N);
        }
#pragma endscop
        //printf("assignemento_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    if (kind == 599)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      // printf("mccNoIfMod599\n");

      for (int i = N - 1; i >= 0; i--)
      {
        for (int j = i + 1; j <= i + 1; j++)
        {
          if (j < N - 1)
          {
            Q[i][j + 1] = Q[i][j];
            assignment_counter++;
          }
        }
      }
#pragma scop
      for (int c0 = 2; c0 < N; c0 += 1)
        for (int c1 = -N + c0 + 1; c1 <= 0; c1 += 1)
          for (int c2 = c0 - c1; c2 < min(N, 2 * c0 - c1 - 1); c2 += 1)
          {
            if (2 * c0 >= c1 + c2 + 3)
            {
              add_counter++;
              Q[-c1][c2] += Q[-c1][-c0 + c1 + c2 -c1] * Q[-c0 + c1 + c2 -c1 + 1][c2 - 1] * ERT * paired(-c0 + c1 + c2 -c1, c2 - 1);
              //S_0(-c1, c2, -c0 + c1 + c2);
            }
            add_counter++;
            Q[-c1][c2] += Q[-c1][c0 - 2 -c1] * Q[c0 - 2 -c1 + 1][c2 - 1] * ERT * paired(c0 - 2 -c1, c2 - 1);
            //S_0(-c1, c2, c0 - 2);
            if (c1 + c2 == c0)
            {
              assignment_counter++;
              Q[-c1][c0 - c1 + 1] = Q[-c1][c0 - c1];
              //S_1(-c1, c0 - c1);
            }
          }

          //       Q[i][j] += Q[i][k + i] * Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
          //       Q[i][j + 1] = Q[i][j];

#pragma endscop
          //printf("assignemento_counter: %d add counter: %d \n", assignment_counter, add_counter);
    }

    if (kind == 6)
    {
      printf("mccNoIfDapt\n"); // pamietac o zmianie Q1 na Q
      {
        for (int i0 = 1; i0 < N; i0 += 1)
        {
          for (int h1 = 0; h1 <= (N - i0 - 1) / 16; h1 += 1)
          {
            for (int i2 = 16 * h1; i2 <= min(N - i0 - 1, 16 * h1 + 15); i2 += 1)
            {
              Q[i2][i0 + i2] = Q[i2][i0 + i2 - 1];
            }
          }
#pragma omp parallel for
          for (int ph0 = 0; ph0 <= (N - i0 - 1) / 16; ph0 += 1)
          {
            for (int h1 = (i0 + 1) / 32; h1 <= (i0 - 1) / 16; h1 += 1)
            {
              for (int i2 = max(16 * h1, (i0 + 1) / 2); i2 <= min(i0 - 1, 16 * h1 + 15); i2 += 1)
              {
                for (int i3 = 16 * ph0; i3 <= min(N - i0 - 1, 16 * ph0 + 15); i3 += 1)
                {
                  if (2 * i2 >= i0 + 1)
                  {
                    Q[i3][i0 + i3] += (((Q[i3][i0 - i2 + i3 - 1] * Q[i0 - i2 + i3][i0 + i3 - 1]) * ERT) * paired((i0 - i2 + i3 - 1), (i0 + i3 - 1)));
                  }
                  Q[i3][i0 + i3] += (((Q[i3][i2 + i3 - 1] * Q[i2 + i3][i0 + i3 - 1]) * ERT) * paired((i2 + i3 - 1), (i0 + i3 - 1)));
                }
              }
            }
          }
        }
      }
    }

    if (kind == 61)
    {
      printf("mccNoIfDapt2\n"); // pamietac o zmianie Q1 na Q
      {
        if (l >= 0 && l <= 5)
        {
          for (int w0 = max(floord(-N + 2, 16), floord(-N - 14, 32)); w0 <= floord(N - 1, 32); w0 += 1)
          {
#pragma omp parallel for
            for (int h0 = max(-((N + 13) / 16), w0 - (N + 31) / 32 + 1); h0 <= min(min(0, w0), 2 * w0 + 2); h0 += 1)
            {
              for (int t0 = max(0, 4 * w0 - 2 * h0); t0 <= min(min(min(4 * w0 - 4 * h0 + 3, 4 * w0 - 2 * h0 + 5), 2 * h0 + (N + 6) / 8 + 1), (N - 1) / 8); t0 += 1)
              {
                for (int i0 = max(max(max(max(-N + 2, -32 * w0 + 32 * h0 - 30), 16 * h0), -32 * w0 + 32 * h0 + 8 * t0 - 31), -N + 8 * t0 + 1); i0 <= min(min(0, 16 * h0 + 15), -32 * w0 + 32 * h0 + 8 * t0 + 7); i0 += 1)
                {
                  for (int i1 = max(max(32 * w0 - 32 * h0, 8 * t0 - i0), -i0 + 1); i1 <= min(min(N - 1, 32 * w0 - 32 * h0 + 31), 8 * t0 - i0 + 7); i1 += 1)
                  {
                    Q[-i0][i1] = Q[-i0][i1 - 1];
                    for (int i2 = 0; i2 < -l + i0 + i1; i2 += 1)
                    {
                      Q[-i0][i1] += (((Q[-i0][-i0 + i2] * Q[-i0 + i2 + 1][i1 - 1]) * ERT) * paired((-i0 + i2), (i1 - 1)));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (kind == 62)
    {
      printf("mccNoIfDapt2\n"); // pamietac o zmianie Q1 na Q
      {
        for (int w0 = 0; w0 <= floord(3 * N - 3, 32); w0 += 1) {
  #pragma omp parallel for
  for (int h0 = max(0, w0 - (N + 31) / 32 + 1); h0 <= min((N - 1) / 16, w0 - (w0 + 1) / 3); h0 += 1) {
    for (int t0 = max(4 * h0, 4 * w0 - 2 * h0); t0 <= min(min(4 * w0 - 2 * h0 + 5, (N - 1) / 4), 2 * h0 + (N + 6) / 8 + 1); t0 += 1) {
      for (int i0 = max(max(max(1, 16 * h0), -32 * w0 + 32 * h0 + 8 * t0 - 31), -N + 8 * t0 + 1); i0 <= min(min(min(N - 1, 16 * h0 + 15), 4 * t0 + 3), -32 * w0 + 32 * h0 + 8 * t0 + 7); i0 += 1) {
        for (int i1 = max(max(32 * w0 - 32 * h0, i0), 8 * t0 - i0); i1 <= min(min(N - 1, 32 * w0 - 32 * h0 + 31), 8 * t0 - i0 + 7); i1 += 1) {
          Q[-i0 + i1][i1] = Q[-i0 + i1][i1 - 1];
          for (int i2 = (i0 + 1) / 2; i2 < i0; i2 += 1) {
            if (2 * i2 >= i0 + 1) {
              Q[-i0 + i1][i1] += (((Q[-i0 + i1][i1 - i2 - 1] * Q[i1 - i2][i1 - 1]) * ERT) * paired((i1 - i2 - 1), (i1 - 1)));
            }
            Q[-i0 + i1][i1] += (((Q[-i0 + i1][-i0 + i1 + i2 - 1] * Q[-i0 + i1 + i2][i1 - 1]) * ERT) * paired((-i0 + i1 + i2 - 1), (i1 - 1)));
          }
        }
      }
    }
  }
}

      }
    }


    if (kind == 7)
    {
// Q[i][j] = Q[i][j-1];
// Q[i][j] +=  Q[i][k+i] * Q[k+i+1][j-1] * ERT * paired(k+i,j-1);

      printf("mccIfPerfectMod\n"); // pamietac o zmianie Q1 na Q
      {
        for (int c0 = 1; c0 < N - 1; c0 += 1)
        {
          for (int c1 = 0; c1 < N - c0 - 1; c1 += 1)
          {
            if (c0 >= 2)
            {
              Q[c1][c0 + c1 + 1] = Q[c1][c0 + c1 + 1-1];
              //S_0(c1, c0 + c1 + 1, 0);
            }
            else
            {
              for (int c2 = c1 + 1; c2 <= c1 + 2; c2 += 1)
                Q[c1][c2] = Q[c1][c2-1];
                //S_0(c1, c2, 0);
            }
            for (int c2 = c0 + c1 + 1; c2 <= min(N - 1, 2 * c0 + c1); c2 += 1)
            {
              if (2 * c0 + c1 >= c2 + 1)
                Q[c1][c2] +=  Q[c1][-c0 - c1 + c2 - 1 +c1] * Q[-c0 - c1 + c2 - 1+c1+1][c2-1] * ERT * paired(-c0 - c1 + c2 - 1+c1, c2-1); // 
                // Q[i][j] +=  Q[i][k+i] * Q[k+i+1][j-1] * ERT * paired(k+i,j-1);

                //S_1(c1, c2, -c0 - c1 + c2 - 1);
              Q[c1][c2] +=  Q[c1][c0 - 1+c1] * Q[c0 - 1+c1+1][c2-1] * ERT * paired(c0 - 1+c1,c2-1);
              //S_1(c1, c2, c0 - 1);
            }
          }
          if (c0 == 1)
            Q[N-2][N-1] = Q[N - 2][N - 1-1];
            //S_0(N - 2, N - 1, 0);
        }
        if (N == 2)
          Q[0][1] = Q[0][1-1];
          //S_0(0, 1, 0);
      }
    }

    if (kind == 8)
    {
      int assignment_counter = 0;
      int add_counter = 0;
      printf("mcc3DModPerfectBadsss\n"); // pamietac o zmianie Q1 na Q
      {
        {
  for (int c1 = 0; c1 < N - 1; c1 += 1) {
    //S_0(c1, c1 + 1, 0, 0);
    assignment_counter++;
    if (N >= c1 + 3)
      add_counter++;
      //S_1(c1, c1 + 2, 0, 1);
  }
  for (int c0 = 2; c0 < N - 1; c0 += 1)
    for (int c1 = 0; c1 < N - c0 - 1; c1 += 1)
      for (int c2 = c0 + c1 + 1; c2 <= min(N - 1, 2 * c0 + c1); c2 += 1) {
        if (2 * c0 + c1 >= c2 + 1)
          add_counter++;
          //S_1(c1, c2, -c0 - c1 + c2 - 1, 1);
        //S_1(c1, c2, c0 - 1, 1);
        add_counter++;
      }
}
               
        printf("assignement_counter: %d add counter: %d \n", assignment_counter, add_counter);
      }
    }

    if (kind == 9 )
    {
      printf("mccIfPerfectModDapt\n");
        for (int w0 = -((N + 13) / 16); w0 <= (N - 1) / 16; w0 += 1) {
    #pragma omp parallel for
    for (int h0 = max(0, w0); h0 <= (N + 16 * w0 + 14) / 32; h0 += 1) {
      for (int i0 = max(1, 16 * h0); i0 <= min(min(N - 2, N + 16 * w0 - 16 * h0 + 13), 16 * h0 + 15); i0 += 1) {
        for (int i1 = max(0, -16 * w0 + 16 * h0 - 15); i1 <= min(-16 * w0 + 16 * h0, N - i0 - 1); i1 += 1) {
          for (int i2 = i0 + i1; i2 <= min(N - 1, 2 * i0 + i1); i2 += 1) {
            if (2 * i0 + i1 >= i2 + 1) {
              Q[i1][i2] = Q[i1][i2 - 1];
            }
            Q[i1][i2] = Q[i1][i2 - 1];
            if (i2 >= i0 + i1 + 1 && 2 * i0 + i1 >= i2 + 1) {
              Q[i1][i2] += (((Q[i1][-i0 - i1 + i2 - 1 + i1] * Q[-i0 + i2][i2 - 1]) * ERT) * paired((-i0 + i2 - 1), (i2 - 1)));
            }
            if (i2 >= i0 + i1 + 1) {
              Q[i1][i2] += (((Q[i1][i0 + i1 - 1] * Q[i0 + i1][i2 - 1]) * ERT) * paired((i0 + i1 - 1), (i2 - 1)));
            }
          }
        }
      }
      if (h0 >= w0 + 1 && 32 * h0 + 1 >= N + 16 * w0) {
        Q[-16 * w0 + 16 * h0 - 15][N - 1] = Q[-16 * w0 + 16 * h0 - 15][N - 2];
        Q[-16 * w0 + 16 * h0 - 15][N - 1] = Q[-16 * w0 + 16 * h0 - 15][N - 2];
      } else if (16 * w0 + 16 >= N && h0 == w0) {
        Q[0][N - 1] = Q[0][N - 2];
        Q[0][N - 1] = Q[0][N - 2];
      }
    }
  }

    }

    //    FILE *f = fopen("client.data", "wb");
    //    fwrite(Q, sizeof(double), sizeof(double)*DIM*DIM, f);



    double stop = omp_get_wtime();
    //printf("%.4f\n",stop - start);
  char bufferQ[256];
  char bufferQ1[256];
  int error_counter = 0;

  if (CHECK_VALID == 1)
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) // TODO: float to str->strncmp
      {
        //printf("Q:%g, Q1:%g i:%d, j:%d \n", Q[i][j], Q1[i][j], i, j);
        if (fabs(Q[i][j] - Q1[i][j]) > 1e-9)
        {
          sprintf(bufferQ, "%g", Q[i][j]);
          sprintf(bufferQ1, "%g", Q1[i][j]);
          // printf("Q:%s, Q1:%s i:%d, j:%d \n", bufferQ, bufferQ1, i, j);
          // printf("%g \n", fabs(Q[i][j] - Q1[i][j]));
          if (strcmp(bufferQ1, bufferQ) != 0)
          {
            if(error_counter<1)
              printf("First error %.18f %.18f -- %d %d\n", Q[i][j], Q1[i][j], i, j);
            // printf("error %s %s -- str      ", bufferQ, bufferQ1);
            // printf("error %g %g -- g      ", Q[i][j], Q1[i][j]);
            // printf("diff  %g\n", Q[i][j] - Q1[i][j]);
            error_counter++;
            //exit(0);
          }
        }
      }
    printf("Error count: %d\n", error_counter);
    return 0;

}
