#include <math.h>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))

bool is_files_equal(const char * , const char * ) ;
void print_files_equal(const char * , const char * );
void write_results_full(int, double, char);
void write_results(int, double );
void write_results_last(int, double );
int** allocate_matrix(int) ;

int **get_full_copy(int ** , int);
int* allocate_vector(int);
int *get_vector_copy(int *, int);

void deallocate_matrix(int **, int ) ;
void print_vector(int* , int );
void print_matrix(int**, int);

/*
  *******************************************************************
  *   Kernel 6 -- general linear recurrence equations
  *******************************************************************
  *    DO  6  L= 1,Loop
  *    DO  6  i= 2,n
  *    DO  6  k= 1,i-1
  *        W(i)= W(i)  + B(i,k) * W(i-k)
  *  6 CONTINUE
*/

void OBST_01_Base(int _PB_N, int **input_path)
{
  int** path = get_full_copy(input_path, _PB_N);

//  double start = omp_get_wtime();

  for (int k = 0; k < _PB_N; k++)
  {
    for (int i = 0; i < _PB_N; i++)
    {
      for (int j = 0; j < _PB_N; j++)
      {
        path[i][j] = path[i][j] < path[i][k] + path[k][j] ? path[i][j] : path[i][k] + path[k][j];
      }
    }
  }
  //double execution_time = omp_get_wtime() - start;

  // printf("BAS_01: %lf\n", execution_time);
  // write_results(_PB_N, execution_time);
  print_matrix(path, _PB_N);
  deallocate_matrix(path, _PB_N);
  return;
}

void OBST_02_Mod(int _PB_N, int **input_path)
{
  int** path = get_full_copy(input_path, _PB_N);
  if (_PB_N >= 2)
    for (int c0 = 0; c0 < _PB_N; c0 += 1)
    {
      if (_PB_N >= c0 + 2)
      {
        for (int c1 = 0; c1 < _PB_N; c1 += 1)
        {
          if (c0 >= 1 && c1 >= c0 + 2 && _PB_N >= c1 + 2)
          {
            path[c0][c1] = path[c0][c1] < path[c0][c0] + path[c0][c1] ? path[c0][c1] : path[c0][c0] + path[c0][c1];
          }
          else if (_PB_N >= c0 + 3 && c1 + 1 == _PB_N)
          {
            for (int c2 = 0; c2 < c0; c2 += 1)
            {
              if (c2 + 1 == c0)
                path[c0][_PB_N - 1] = path[c0][_PB_N - 1] < path[c0][c0] + path[c0][_PB_N - 1] ? path[c0][_PB_N - 1] : path[c0][c0] + path[c0][_PB_N - 1];
              path[c0 + 1][c0] = path[c0 + 1][c0] < path[c0 + 1][c2] + path[c2][c0] ? path[c0 + 1][c0] : path[c0 + 1][c2] + path[c2][c0];
            }
          }
          else
          {
            if (c0 == 0 && c1 == 0)
              path[0][0] = path[0][0] < path[0][0] + path[0][0] ? path[0][0] : path[0][0] + path[0][0];
            if (c1 + 1 >= c0 && c0 + 1 >= c1)
              for (int c2 = 0; c2 < min(c0, c1); c2 += 1)
              {
                if (c1 >= c0)
                {
                  if (c1 == c0 + 1 && c2 + 1 == c0)
                    path[c0][c0 + 1] = path[c0][c0 + 1] < path[c0][c0] + path[c0][c0 + 1] ? path[c0][c0 + 1] : path[c0][c0] + path[c0][c0 + 1];
                  if (c1 == c0 + 1)
                  {
                    for (int c4 = 0; c4 < c0; c4 += 1)
                      path[c0 + 1][c4] = path[c0 + 1][c4] < path[c0 + 1][c2] + path[c2][c4] ? path[c0 + 1][c4] : path[c0 + 1][c2] + path[c2][c4];
                    if (c0 + 2 == _PB_N)
                      path[_PB_N - 1][_PB_N - 2] = path[_PB_N - 1][_PB_N - 2] < path[_PB_N - 1][c2] + path[c2][_PB_N - 2] ? path[_PB_N - 1][_PB_N - 2] : path[_PB_N - 1][c2] + path[c2][_PB_N - 2];
                  }
                  else
                  {
                    if (c2 + 1 == c0)
                      path[c0][c0] = path[c0][c0] < path[c0][c0] + path[c0][c0] ? path[c0][c0] : path[c0][c0] + path[c0][c0];
                    for (int c4 = c0 + 1; c4 < _PB_N; c4 += 1)
                      path[c0][c4] = path[c0][c4] < path[c0][c2] + path[c2][c4] ? path[c0][c4] : path[c0][c2] + path[c2][c4];
                  }
                }
                else
                  path[c0][c0] = path[c0][c0] < path[c0][c2] + path[c2][c0] ? path[c0][c0] : path[c0][c2] + path[c2][c0];
              }
            if (c0 >= c1 + 1)
            {
              if (c0 >= c1 + 2)
              {
                path[c0][c1] = path[c0][c1] < path[c0][c0] + path[c0][c1] ? path[c0][c1] : path[c0][c0] + path[c0][c1];
              }
              else
                for (int c4 = c0 - 1; c4 <= c0; c4 += 1)
                  path[c0][c4] = path[c0][c4] < path[c0][2 * c0 - c4 - 1] + path[2 * c0 - c4 - 1][c4] ? path[c0][c4] : path[c0][2 * c0 - c4 - 1] + path[2 * c0 - c4 - 1][c4];
            }
          }
          for (int c5 = c0 + 1; c5 < _PB_N; c5 += 1)
            path[c0][c1] = path[c0][c1] < path[c0][c5] + path[c5][c1] ? path[c0][c1] : path[c0][c5] + path[c5][c1];
          if (c0 == 0 && c1 == 0)
            for (int c4 = 1; c4 < _PB_N; c4 += 1)
              path[0][c4] = path[0][c4] < path[0][0] + path[0][c4] ? path[0][c4] : path[0][0] + path[0][c4];
          if (_PB_N == 2 && c0 == 0 && c1 == 1)
          {
            path[1][0] = path[1][0] < path[1][0] + path[0][0] ? path[1][0] : path[1][0] + path[0][0];
          }
          else if (_PB_N >= c0 + 3 && c1 + 1 == _PB_N)
          {
            path[c0 + 1][c0] = path[c0 + 1][c0] < path[c0 + 1][c0] + path[c0][c0] ? path[c0 + 1][c0] : path[c0 + 1][c0] + path[c0][c0];
          }
          else if (c0 >= 1 && c1 == c0 + 1)
          {
            for (int c4 = 0; c4 < c0; c4 += 1)
              path[c0 + 1][c4] = path[c0 + 1][c4] < path[c0 + 1][c0] + path[c0][c4] ? path[c0 + 1][c4] : path[c0 + 1][c0] + path[c0][c4];
            if (c0 + 2 == _PB_N)
              path[_PB_N - 1][_PB_N - 2] = path[_PB_N - 1][_PB_N - 2] < path[_PB_N - 1][_PB_N - 2] + path[_PB_N - 2][_PB_N - 2] ? path[_PB_N - 1][_PB_N - 2] : path[_PB_N - 1][_PB_N - 2] + path[_PB_N - 2][_PB_N - 2];
          }
        }
      }
      else
      {
        if (_PB_N == 2)
        {
          for (int c4 = 0; c4 <= 1; c4 += 1)
            path[1][c4] = path[1][c4] < path[1][-c4 + 1] + path[-c4 + 1][c4] ? path[1][c4] : path[1][-c4 + 1] + path[-c4 + 1][c4];
        }
        else
          for (int c1 = 0; c1 < _PB_N - 1; c1 += 1)
          {
            if (c1 + 2 == _PB_N)
              for (int c2 = 0; c2 < _PB_N - 2; c2 += 1)
                path[_PB_N - 1][_PB_N - 1] = path[_PB_N - 1][_PB_N - 1] < path[_PB_N - 1][c2] + path[c2][_PB_N - 1] ? path[_PB_N - 1][_PB_N - 1] : path[_PB_N - 1][c2] + path[c2][_PB_N - 1];
            if (_PB_N >= c1 + 3)
            {
              path[_PB_N - 1][c1] = path[_PB_N - 1][c1] < path[_PB_N - 1][_PB_N - 1] + path[_PB_N - 1][c1] ? path[_PB_N - 1][c1] : path[_PB_N - 1][_PB_N - 1] + path[_PB_N - 1][c1];
            }
            else
              for (int c4 = _PB_N - 2; c4 < _PB_N; c4 += 1)
                path[_PB_N - 1][c4] = path[_PB_N - 1][c4] < path[_PB_N - 1][2 * _PB_N - c4 - 3] + path[2 * _PB_N - c4 - 3][c4] ? path[_PB_N - 1][c4] : path[_PB_N - 1][2 * _PB_N - c4 - 3] + path[2 * _PB_N - c4 - 3][c4];
          }
        path[_PB_N - 1][_PB_N - 1] = path[_PB_N - 1][_PB_N - 1] < path[_PB_N - 1][_PB_N - 1] + path[_PB_N - 1][_PB_N - 1] ? path[_PB_N - 1][_PB_N - 1] : path[_PB_N - 1][_PB_N - 1] + path[_PB_N - 1][_PB_N - 1];
      }
    } 
  print_matrix(path, _PB_N);
  deallocate_matrix(path, _PB_N);
  return;
}


void OBST_03_Mod(int _PB_N, int **input_path)
{
  int** path = get_full_copy(input_path, _PB_N);
  if (_PB_N >= 2) {
  for (int c0 = 0; c0 < _PB_N; c0 += 1) {
    if (_PB_N >= c0 + 2) {
      for (int c1 = 0; c1 < _PB_N; c1 += 1) {
        if (c0 >= 1 && c1 >= c0 + 2 && _PB_N >= c1 + 2) {
          path[c0][c1] = path[c0][c1] < path[c0][c0] + path[c0][c1] ? path[c0][c1] : path[c0][c0] + path[c0][c1];;
        } else if (_PB_N >= c0 + 3 && c1 + 1 == _PB_N) {
          for (int c2 = 0; c2 < c0; c2 += 1) {
            if (c2 + 1 == c0)
              path[c0][_PB_N - 1] = path[c0][_PB_N - 1] < path[c0][c0] + path[c0][_PB_N - 1] ? path[c0][_PB_N - 1] : path[c0][c0] + path[c0][_PB_N - 1];;
            path[c0 + 1][c0] = path[c0 + 1][c0] < path[c0 + 1][c2] + path[c2][c0] ? path[c0 + 1][c0] : path[c0 + 1][c2] + path[c2][c0];;
          }
        } else {
          if (c0 == 0 && c1 == 0)
            path[0][0] = path[0][0] < path[0][0] + path[0][0] ? path[0][0] : path[0][0] + path[0][0];;
          if (c1 + 1 >= c0 && c0 + 1 >= c1)
            for (int c2 = 0; c2 < min(c0, c1); c2 += 1) {
              if (c1 >= c0) {
                if (c1 == c0 + 1 && c2 + 1 == c0)
                  path[c0][c0 + 1] = path[c0][c0 + 1] < path[c0][c0] + path[c0][c0 + 1] ? path[c0][c0 + 1] : path[c0][c0] + path[c0][c0 + 1];;
                if (c1 == c0 + 1) {
                  for (int c4 = 0; c4 < c0; c4 += 1)
                    path[c0 + 1][c4] = path[c0 + 1][c4] < path[c0 + 1][c2] + path[c2][c4] ? path[c0 + 1][c4] : path[c0 + 1][c2] + path[c2][c4];;
                  if (c0 + 2 == _PB_N)
                    path[_PB_N - 1][_PB_N - 2] = path[_PB_N - 1][_PB_N - 2] < path[_PB_N - 1][c2] + path[c2][_PB_N - 2] ? path[_PB_N - 1][_PB_N - 2] : path[_PB_N - 1][c2] + path[c2][_PB_N - 2];;
                } else {
                  if (c2 + 1 == c0)
                    path[c0][c0] = path[c0][c0] < path[c0][c0] + path[c0][c0] ? path[c0][c0] : path[c0][c0] + path[c0][c0];;
                  for (int c4 = c0 + 1; c4 < _PB_N; c4 += 1)
                    path[c0][c4] = path[c0][c4] < path[c0][c2] + path[c2][c4] ? path[c0][c4] : path[c0][c2] + path[c2][c4];;
                }
              } else
                path[c0][c0] = path[c0][c0] < path[c0][c2] + path[c2][c0] ? path[c0][c0] : path[c0][c2] + path[c2][c0];;
            }
          if (c0 >= c1 + 1) {
            if (c0 >= c1 + 2) {
              path[c0][c1] = path[c0][c1] < path[c0][c0] + path[c0][c1] ? path[c0][c1] : path[c0][c0] + path[c0][c1];;
            } else
              for (int c4 = c0 - 1; c4 <= c0; c4 += 1)
                path[c0][c4] = path[c0][c4] < path[c0][2 * c0 - c4 - 1] + path[2 * c0 - c4 - 1][c4] ? path[c0][c4] : path[c0][2 * c0 - c4 - 1] + path[2 * c0 - c4 - 1][c4];;
          }
        }
        for (int c5 = c0 + 1; c5 < _PB_N; c5 += 1)
          path[c0][c1] = path[c0][c1] < path[c0][c5] + path[c5][c1] ? path[c0][c1] : path[c0][c5] + path[c5][c1];;
        if (c0 == 0 && c1 == 0)
          for (int c4 = 1; c4 < _PB_N; c4 += 1)
            path[0][c4] = path[0][c4] < path[0][0] + path[0][c4] ? path[0][c4] : path[0][0] + path[0][c4];;
        if (_PB_N >= c0 + 3 && c1 + 1 == _PB_N) {
          path[c0 + 1][c0] = path[c0 + 1][c0] < path[c0 + 1][c0] + path[c0][c0] ? path[c0 + 1][c0] : path[c0 + 1][c0] + path[c0][c0];;
        } else if (c1 == c0 + 1)
          for (int c4 = 0; c4 < c0; c4 += 1)
            path[c0 + 1][c4] = path[c0 + 1][c4] < path[c0 + 1][c0] + path[c0][c4] ? path[c0 + 1][c4] : path[c0 + 1][c0] + path[c0][c4];;
      }
    } else {
      if (_PB_N == 2) {
        for (int c4 = 0; c4 <= 1; c4 += 1)
          path[1][c4] = path[1][c4] < path[1][-c4 + 1] + path[-c4 + 1][c4] ? path[1][c4] : path[1][-c4 + 1] + path[-c4 + 1][c4];;
      } else
        for (int c1 = 0; c1 < _PB_N - 1; c1 += 1) {
          if (c1 + 2 == _PB_N)
            for (int c2 = 0; c2 < _PB_N - 2; c2 += 1)
              path[_PB_N - 1][_PB_N - 1] = path[_PB_N - 1][_PB_N - 1] < path[_PB_N - 1][c2] + path[c2][_PB_N - 1] ? path[_PB_N - 1][_PB_N - 1] : path[_PB_N - 1][c2] + path[c2][_PB_N - 1];;
          if (_PB_N >= c1 + 3) {
            path[_PB_N - 1][c1] = path[_PB_N - 1][c1] < path[_PB_N - 1][_PB_N - 1] + path[_PB_N - 1][c1] ? path[_PB_N - 1][c1] : path[_PB_N - 1][_PB_N - 1] + path[_PB_N - 1][c1];;
          } else
            for (int c4 = _PB_N - 2; c4 < _PB_N; c4 += 1)
              path[_PB_N - 1][c4] = path[_PB_N - 1][c4] < path[_PB_N - 1][2 * _PB_N - c4 - 3] + path[2 * _PB_N - c4 - 3][c4] ? path[_PB_N - 1][c4] : path[_PB_N - 1][2 * _PB_N - c4 - 3] + path[2 * _PB_N - c4 - 3][c4];;
        }
      path[_PB_N - 1][_PB_N - 1] = path[_PB_N - 1][_PB_N - 1] < path[_PB_N - 1][_PB_N - 1] + path[_PB_N - 1][_PB_N - 1] ? path[_PB_N - 1][_PB_N - 1] : path[_PB_N - 1][_PB_N - 1] + path[_PB_N - 1][_PB_N - 1];;
    }
  }
  path[_PB_N - 1][_PB_N - 2] = path[_PB_N - 1][_PB_N - 2] < path[_PB_N - 1][_PB_N - 2] + path[_PB_N - 2][_PB_N - 2] ? path[_PB_N - 1][_PB_N - 2] : path[_PB_N - 1][_PB_N - 2] + path[_PB_N - 2][_PB_N - 2];;
}
  print_matrix(path, _PB_N);
  deallocate_matrix(path, _PB_N);
  return;
}

void OBST_04_Mod(int _PB_N, int **input_path)
{
  int** path = get_full_copy(input_path, _PB_N);
 {
  if (_PB_N >= 3)
    path[0][0] = path[0][0] < path[0][0] + path[0][0] ? path[0][0] : path[0][0] + path[0][0];;
  if (_PB_N >= 2) {
    if (_PB_N == 2)
      path[0][0] = path[0][0] < path[0][0] + path[0][0] ? path[0][0] : path[0][0] + path[0][0];;
    path[_PB_N - 1][_PB_N - 2] = path[_PB_N - 1][_PB_N - 2] < path[_PB_N - 1][_PB_N - 2] + path[_PB_N - 2][_PB_N - 2] ? path[_PB_N - 1][_PB_N - 2] : path[_PB_N - 1][_PB_N - 2] + path[_PB_N - 2][_PB_N - 2];;
  }
  for (int c0 = 0; c0 < 111 * _PB_N - 210; c0 += 1) {
    for (int c1 = max(0, floord(-10 * _PB_N + c0 + 9, 101) + 1); c1 <= min(_PB_N - 2, c0 / 101); c1 += 1) {
      if ((c0 - c1) % 10 == 0) {
        for (int c3 = c1 + 1; c3 < _PB_N; c3 += 1)
          path[c1][(c0 - 101 * c1) / 10] = path[c1][(c0 - 101 * c1) / 10] < path[c1][c3] + path[c3][(c0 - 101 * c1) / 10] ? path[c1][(c0 - 101 * c1) / 10] : path[c1][c3] + path[c3][(c0 - 101 * c1) / 10];;
        if (c0 == 0 && c1 == 0) {
          for (int c2 = 1; c2 < _PB_N; c2 += 1)
            path[0][c2] = path[0][c2] < path[0][0] + path[0][c2] ? path[0][c2] : path[0][0] + path[0][c2];;
        } else if (c1 >= 1 && 10 * _PB_N + 101 * c1 >= c0 + 110 && (c0 - c1) % 10 == 0) {
          path[c1][((c0 - 101 * c1) / 10) + 10] = path[c1][((c0 - 101 * c1) / 10) + 10] < path[c1][c1] + path[c1][((c0 - 101 * c1) / 10) + 10] ? path[c1][((c0 - 101 * c1) / 10) + 10] : path[c1][c1] + path[c1][((c0 - 101 * c1) / 10) + 10];;
          if (111 * c1 >= c0 + 110 && (91 * c0 - c1 + 10) % 100 == 0) {
            path[c1][c1] = path[c1][c1] < path[c1][(c0 - 11 * c1 + 10) / 100] + path[(c0 - 11 * c1 + 10) / 100][c1] ? path[c1][c1] : path[c1][(c0 - 11 * c1 + 10) / 100] + path[(c0 - 11 * c1 + 10) / 100][c1];;
          } else if (111 * c1 >= c0 + 100 && (91 * c0 - c1) % 100 == 0)
            for (int c2 = c1 + 1; c2 < _PB_N; c2 += 1)
              path[c1][c2] = path[c1][c2] < path[c1][(c0 - 11 * c1) / 100] + path[(c0 - 11 * c1) / 100][c2] ? path[c1][c2] : path[c1][(c0 - 11 * c1) / 100] + path[(c0 - 11 * c1) / 100][c2];;
        }
      }
      if (111 * c1 >= c0 + 11) {
        if ((91 * c0 - c1 + 91) % 100 == 0)
          for (int c2 = 0; c2 < c1 - 1; c2 += 1)
            path[c1][c2] = path[c1][c2] < path[c1][(c0 - 11 * c1 + 1) / 100] + path[(c0 - 11 * c1 + 1) / 100][c2] ? path[c1][c2] : path[c1][(c0 - 11 * c1 + 1) / 100] + path[(c0 - 11 * c1 + 1) / 100][c2];;
        if (c0 + 11 >= 10 * _PB_N + c1 && (90 * _PB_N + c0 - c1 + 11) % 100 == 0)
          path[c1][c1 - 1] = path[c1][c1 - 1] < path[c1][(-10 * _PB_N + c0 - c1 + 11) / 100] + path[(-10 * _PB_N + c0 - c1 + 11) / 100][c1 - 1] ? path[c1][c1 - 1] : path[c1][(-10 * _PB_N + c0 - c1 + 11) / 100] + path[(-10 * _PB_N + c0 - c1 + 11) / 100][c1 - 1];;
      } else if (c0 + 11 >= 10 * _PB_N + c1 && 10 * _PB_N + 101 * c1 >= c0 + 111 && (90 * _PB_N + c0 - c1 + 11) % 100 == 0)
        path[c1][c1 - 1] = path[c1][c1 - 1] < path[c1][(-10 * _PB_N + c0 - c1 + 11) / 100] + path[(-10 * _PB_N + c0 - c1 + 11) / 100][c1 - 1] ? path[c1][c1 - 1] : path[c1][(-10 * _PB_N + c0 - c1 + 11) / 100] + path[(-10 * _PB_N + c0 - c1 + 11) / 100][c1 - 1];;
    }
    if (c0 + 101 >= 101 * _PB_N) {
      if ((_PB_N - c0 + 9) % 10 == 0) {
        path[_PB_N - 1][(-101 * _PB_N + c0 + 201) / 10] = path[_PB_N - 1][(-101 * _PB_N + c0 + 201) / 10] < path[_PB_N - 1][_PB_N - 1] + path[_PB_N - 1][(-101 * _PB_N + c0 + 201) / 10] ? path[_PB_N - 1][(-101 * _PB_N + c0 + 201) / 10] : path[_PB_N - 1][_PB_N - 1] + path[_PB_N - 1][(-101 * _PB_N + c0 + 201) / 10];;
        if ((11 * _PB_N - c0 + 79) % 100 == 0)
          path[_PB_N - 1][_PB_N - 1] = path[_PB_N - 1][_PB_N - 1] < path[_PB_N - 1][(-11 * _PB_N + c0 + 21) / 100] + path[(-11 * _PB_N + c0 + 21) / 100][_PB_N - 1] ? path[_PB_N - 1][_PB_N - 1] : path[_PB_N - 1][(-11 * _PB_N + c0 + 21) / 100] + path[(-11 * _PB_N + c0 + 21) / 100][_PB_N - 1];;
      }
      if ((11 * _PB_N - c0 + 88) % 100 == 0) {
        for (int c2 = 0; c2 < _PB_N - 2; c2 += 1)
          path[_PB_N - 1][c2] = path[_PB_N - 1][c2] < path[_PB_N - 1][(-11 * _PB_N + c0 + 12) / 100] + path[(-11 * _PB_N + c0 + 12) / 100][c2] ? path[_PB_N - 1][c2] : path[_PB_N - 1][(-11 * _PB_N + c0 + 12) / 100] + path[(-11 * _PB_N + c0 + 12) / 100][c2];;
        if (111 * _PB_N >= c0 + 312)
          path[_PB_N - 1][_PB_N - 2] = path[_PB_N - 1][_PB_N - 2] < path[_PB_N - 1][(-11 * _PB_N + c0 + 12) / 100] + path[(-11 * _PB_N + c0 + 12) / 100][_PB_N - 2] ? path[_PB_N - 1][_PB_N - 2] : path[_PB_N - 1][(-11 * _PB_N + c0 + 12) / 100] + path[(-11 * _PB_N + c0 + 12) / 100][_PB_N - 2];;
      }
    }
    for (int c1 = max(c0 / 101 + 1, floord(-10 * _PB_N + c0 + 8, 101) + 2); c1 <= min(_PB_N - 1, (c0 + 10) / 11); c1 += 1) {
      if (c0 + 100 >= 101 * c1 && (c0 - 11 * c1 + 100) % 100 >= 10 && (c0 - c1) % 10 == 0)
        path[c1][((c0 - 101 * c1) / 10) + 10] = path[c1][((c0 - 101 * c1) / 10) + 10] < path[c1][c1] + path[c1][((c0 - 101 * c1) / 10) + 10] ? path[c1][((c0 - 101 * c1) / 10) + 10] : path[c1][c1] + path[c1][((c0 - 101 * c1) / 10) + 10];;
      if ((91 * c0 - c1 + 10) % 100 == 0)
        path[c1][c1] = path[c1][c1] < path[c1][(c0 - 11 * c1 + 10) / 100] + path[(c0 - 11 * c1 + 10) / 100][c1] ? path[c1][c1] : path[c1][(c0 - 11 * c1 + 10) / 100] + path[(c0 - 11 * c1 + 10) / 100][c1];;
      if ((91 * c0 - c1 + 91) % 100 == 0)
        for (int c2 = 0; c2 < c1 - 1; c2 += 1)
          path[c1][c2] = path[c1][c2] < path[c1][(c0 - 11 * c1 + 1) / 100] + path[(c0 - 11 * c1 + 1) / 100][c2] ? path[c1][c2] : path[c1][(c0 - 11 * c1 + 1) / 100] + path[(c0 - 11 * c1 + 1) / 100][c2];;
      if (110 * _PB_N + c1 >= c0 + 311 && c0 + 11 >= 10 * _PB_N + c1 && (90 * _PB_N + c0 - c1 + 11) % 100 == 0)
        path[c1][c1 - 1] = path[c1][c1 - 1] < path[c1][(-10 * _PB_N + c0 - c1 + 11) / 100] + path[(-10 * _PB_N + c0 - c1 + 11) / 100][c1 - 1] ? path[c1][c1 - 1] : path[c1][(-10 * _PB_N + c0 - c1 + 11) / 100] + path[(-10 * _PB_N + c0 - c1 + 11) / 100][c1 - 1];;
      if ((91 * c0 - c1) % 100 == 0) {
        if (c0 + 100 >= 101 * c1)
          path[c1][((c0 - 101 * c1) / 10) + 10] = path[c1][((c0 - 101 * c1) / 10) + 10] < path[c1][c1] + path[c1][((c0 - 101 * c1) / 10) + 10] ? path[c1][((c0 - 101 * c1) / 10) + 10] : path[c1][c1] + path[c1][((c0 - 101 * c1) / 10) + 10];;
        for (int c2 = c1 + 1; c2 < _PB_N; c2 += 1)
          path[c1][c2] = path[c1][c2] < path[c1][(c0 - 11 * c1) / 100] + path[(c0 - 11 * c1) / 100][c2] ? path[c1][c2] : path[c1][(c0 - 11 * c1) / 100] + path[(c0 - 11 * c1) / 100][c2];;
      }
    }
  }
}
  print_matrix(path, _PB_N);
  deallocate_matrix(path, _PB_N);
  return;
}





#define PERFORMANCE_TEST 1

void make_work_one_size(const int ZMAX)
{
  int** input_b = allocate_matrix(ZMAX);
  for (int i = 0; i < ZMAX; i++)
    for (int j = 0; j < ZMAX; j++)
      input_b[i][j] = 0;
  const char* seqTest = "1234432432123412";
#if PERFORMANCE_TEST==1
  for (int i = 0; i < ZMAX; i++)
    for (int j = 0; j < ZMAX; j++)
      input_b[i][j] = rand()%9+1;
#else
  for (int i = 0; i < ZMAX; i++)
    for (int j = 0; j < ZMAX; j++)
      input_b[i][j] = seqTest[i]-'0';
      
#endif
  OBST_01_Base(ZMAX, input_b);
  OBST_02_Mod(ZMAX, input_b);
  OBST_03_Mod(ZMAX, input_b);
  OBST_04_Mod(ZMAX, input_b);
  print_matrix(input_b, ZMAX);
  deallocate_matrix(input_b, ZMAX);
}


int main()
{
  //TODO: change int to double
#if PERFORMANCE_TEST==1
    const int ZMAX = 55;
    for (int z = 50 ; z <ZMAX ; z+=100)
      make_work_one_size(z);
#else
  const int ZMAX = 16;
  make_work_one_size(ZMAX);
#endif 

  print_files_equal("resMat_1", "resMat_2");
  print_files_equal("resMat_1", "resMat_3");
  print_files_equal("resMat_1", "resMat_4");
  print_files_equal("resMat_2", "resMat_3");
  print_files_equal("resMat_2", "resMat_4");
  
}

bool is_files_equal(const char * filename_template, const char * filename_compared) {
  std::ifstream f1(filename_template, std::ifstream::binary | std::ifstream::ate);
  std::ifstream f2(filename_compared, std::ifstream::binary | std::ifstream::ate);
  if (f1.fail() || f2.fail()) {
    return false;
  }
  if (f1.tellg() != f2.tellg()) {
    return false;
  }
  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator <char >(f1.rdbuf()),
    std::istreambuf_iterator <char >(),
    std::istreambuf_iterator <char >(f2.rdbuf()));
}

void print_files_equal(const char* filename_template, const char* filename_compared)
{
  
  std::cout << "Template: " << filename_template << " compared: " << filename_compared ;
  if (is_files_equal(filename_template, filename_compared))
    std::cout<<" : OK\n";
  else
    std::cout << ": ERROR\n";
}

void write_results_full(int n, double execution_time, char end_char)
{
  FILE* f = fopen("results.txt", "at");
  fprintf(f, "%d:%lf%c", n, execution_time, end_char);
  fclose(f);
}

void write_results(int n, double execution_time)
{
  write_results_full(n, execution_time, ';');
}

void write_results_last(int n, double execution_time)
{
  write_results_full(n, execution_time, '\n');
}

int** allocate_matrix(int N) {
  int** t = (int**)malloc(sizeof(int*) * N);
  for (int i = 0; i < N; i++) {
    t[i] = (int*)malloc(sizeof(int) * N);
  }
  return t;
}

int **get_full_copy(int ** table, int N)
{
  int **S = allocate_matrix(N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      S[i][j] = table[i][j];
  return S;
}

int* allocate_vector(int N) {
  int* t = (int*)malloc(sizeof(int) * N);
  return t;
}

int *get_vector_copy(int *table, int N)
{
  int *S = allocate_vector(N);
  for (int i = 0; i < N; i++)
      S[i] = table[i];
  return S;
}


void deallocate_matrix(int **t, int N) {
  for (int i = 0; i < N; i++) {
    free(t[i]);
  }
  free(t);
}

void print_vector(int* vector, int N) {
  static int fileno=1;
  char filename[12];
  sprintf(filename, "resVec_%d", fileno);
  FILE* f = fopen(filename, "wt");
  for (int i = 0; i < N; i++) {
    fprintf(f, "%d ", vector[i]);
    fprintf(f, "\n");
  }
  fclose(f);
  fileno++;
}

void print_matrix(int** matrix, int N) {
  static int fileno=1;
  char filename[10];
  sprintf(filename, "resMat_%d", fileno);
  FILE* f = fopen(filename, "wt");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++)
      fprintf(f, "%d ", matrix[i][j]);
    fprintf(f, "\n");
  }
  fclose(f);
  fileno++;
}
