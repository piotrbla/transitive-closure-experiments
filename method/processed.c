# 1 "toPreprocess.c"
# 1 "<built-in>"
# 1 "<command-line>"
# 31 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 32 "<command-line>" 2
# 1 "toPreprocess.c"


void OBST_02_Mod(int _PB_N, int **path)
{

if (_PB_N >= 2)
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
        if (_PB_N == 2 && c0 == 0 && c1 == 1) {
          path[1][0] = path[1][0] < path[1][0] + path[0][0] ? path[1][0] : path[1][0] + path[0][0];;
        } else if (_PB_N >= c0 + 3 && c1 + 1 == _PB_N) {
          path[c0 + 1][c0] = path[c0 + 1][c0] < path[c0 + 1][c0] + path[c0][c0] ? path[c0 + 1][c0] : path[c0 + 1][c0] + path[c0][c0];;
        } else if (c0 >= 1 && c1 == c0 + 1) {
          for (int c4 = 0; c4 < c0; c4 += 1)
            path[c0 + 1][c4] = path[c0 + 1][c4] < path[c0 + 1][c0] + path[c0][c4] ? path[c0 + 1][c4] : path[c0 + 1][c0] + path[c0][c4];;
          if (c0 + 2 == _PB_N)
            path[_PB_N - 1][_PB_N - 2] = path[_PB_N - 1][_PB_N - 2] < path[_PB_N - 1][_PB_N - 2] + path[_PB_N - 2][_PB_N - 2] ? path[_PB_N - 1][_PB_N - 2] : path[_PB_N - 1][_PB_N - 2] + path[_PB_N - 2][_PB_N - 2];;
        }
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
}

void OBST_03_Mod(int _PB_N, int **path)
{
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
}

void OBST_04_Mod(int _PB_N, int **path)
{
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
}
