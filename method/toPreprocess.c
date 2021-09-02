#define S_0(i,j,k)  path[i][j] = path[i][j] < path[i][k] + path[k][j] ? path[i][j] : path[i][k] + path[k][j];

void OBST_02_Mod(int _PB_N, int **path)
{

if (_PB_N >= 2)
  for (int c0 = 0; c0 < _PB_N; c0 += 1) {
    if (_PB_N >= c0 + 2) {
      for (int c1 = 0; c1 < _PB_N; c1 += 1) {
        if (c0 >= 1 && c1 >= c0 + 2 && _PB_N >= c1 + 2) {
          S_0(c0, c1, c0);
        } else if (_PB_N >= c0 + 3 && c1 + 1 == _PB_N) {
          for (int c2 = 0; c2 < c0; c2 += 1) {
            if (c2 + 1 == c0)
              S_0(c0, _PB_N - 1, c0);
            S_0(c0 + 1, c0, c2);
          }
        } else {
          if (c0 == 0 && c1 == 0)
            S_0(0, 0, 0);
          if (c1 + 1 >= c0 && c0 + 1 >= c1)
            for (int c2 = 0; c2 < min(c0, c1); c2 += 1) {
              if (c1 >= c0) {
                if (c1 == c0 + 1 && c2 + 1 == c0)
                  S_0(c0, c0 + 1, c0);
                if (c1 == c0 + 1) {
                  for (int c4 = 0; c4 < c0; c4 += 1)
                    S_0(c0 + 1, c4, c2);
                  if (c0 + 2 == _PB_N)
                    S_0(_PB_N - 1, _PB_N - 2, c2);
                } else {
                  if (c2 + 1 == c0)
                    S_0(c0, c0, c0);
                  for (int c4 = c0 + 1; c4 < _PB_N; c4 += 1)
                    S_0(c0, c4, c2);
                }
              } else
                S_0(c0, c0, c2);
            }
          if (c0 >= c1 + 1) {
            if (c0 >= c1 + 2) {
              S_0(c0, c1, c0);
            } else
              for (int c4 = c0 - 1; c4 <= c0; c4 += 1)
                S_0(c0, c4, 2 * c0 - c4 - 1);
          }
        }
        for (int c5 = c0 + 1; c5 < _PB_N; c5 += 1)
          S_0(c0, c1, c5);
        if (c0 == 0 && c1 == 0)
          for (int c4 = 1; c4 < _PB_N; c4 += 1)
            S_0(0, c4, 0);
        if (_PB_N == 2 && c0 == 0 && c1 == 1) {
          S_0(1, 0, 0);
        } else if (_PB_N >= c0 + 3 && c1 + 1 == _PB_N) {
          S_0(c0 + 1, c0, c0);
        } else if (c0 >= 1 && c1 == c0 + 1) {
          for (int c4 = 0; c4 < c0; c4 += 1)
            S_0(c0 + 1, c4, c0);
          if (c0 + 2 == _PB_N)
            S_0(_PB_N - 1, _PB_N - 2, _PB_N - 2);
        }
      }
    } else {
      if (_PB_N == 2) {
        for (int c4 = 0; c4 <= 1; c4 += 1)
          S_0(1, c4, -c4 + 1);
      } else
        for (int c1 = 0; c1 < _PB_N - 1; c1 += 1) {
          if (c1 + 2 == _PB_N)
            for (int c2 = 0; c2 < _PB_N - 2; c2 += 1)
              S_0(_PB_N - 1, _PB_N - 1, c2);
          if (_PB_N >= c1 + 3) {
            S_0(_PB_N - 1, c1, _PB_N - 1);
          } else
            for (int c4 = _PB_N - 2; c4 < _PB_N; c4 += 1)
              S_0(_PB_N - 1, c4, 2 * _PB_N - c4 - 3);
        }
      S_0(_PB_N - 1, _PB_N - 1, _PB_N - 1);
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
          S_0(c0, c1, c0);
        } else if (_PB_N >= c0 + 3 && c1 + 1 == _PB_N) {
          for (int c2 = 0; c2 < c0; c2 += 1) {
            if (c2 + 1 == c0)
              S_0(c0, _PB_N - 1, c0);
            S_0(c0 + 1, c0, c2);
          }
        } else {
          if (c0 == 0 && c1 == 0)
            S_0(0, 0, 0);
          if (c1 + 1 >= c0 && c0 + 1 >= c1)
            for (int c2 = 0; c2 < min(c0, c1); c2 += 1) {
              if (c1 >= c0) {
                if (c1 == c0 + 1 && c2 + 1 == c0)
                  S_0(c0, c0 + 1, c0);
                if (c1 == c0 + 1) {
                  for (int c4 = 0; c4 < c0; c4 += 1)
                    S_0(c0 + 1, c4, c2);
                  if (c0 + 2 == _PB_N)
                    S_0(_PB_N - 1, _PB_N - 2, c2);
                } else {
                  if (c2 + 1 == c0)
                    S_0(c0, c0, c0);
                  for (int c4 = c0 + 1; c4 < _PB_N; c4 += 1)
                    S_0(c0, c4, c2);
                }
              } else
                S_0(c0, c0, c2);
            }
          if (c0 >= c1 + 1) {
            if (c0 >= c1 + 2) {
              S_0(c0, c1, c0);
            } else
              for (int c4 = c0 - 1; c4 <= c0; c4 += 1)
                S_0(c0, c4, 2 * c0 - c4 - 1);
          }
        }
        for (int c5 = c0 + 1; c5 < _PB_N; c5 += 1)
          S_0(c0, c1, c5);
        if (c0 == 0 && c1 == 0)
          for (int c4 = 1; c4 < _PB_N; c4 += 1)
            S_0(0, c4, 0);
        if (_PB_N >= c0 + 3 && c1 + 1 == _PB_N) {
          S_0(c0 + 1, c0, c0);
        } else if (c1 == c0 + 1)
          for (int c4 = 0; c4 < c0; c4 += 1)
            S_0(c0 + 1, c4, c0);
      }
    } else {
      if (_PB_N == 2) {
        for (int c4 = 0; c4 <= 1; c4 += 1)
          S_0(1, c4, -c4 + 1);
      } else
        for (int c1 = 0; c1 < _PB_N - 1; c1 += 1) {
          if (c1 + 2 == _PB_N)
            for (int c2 = 0; c2 < _PB_N - 2; c2 += 1)
              S_0(_PB_N - 1, _PB_N - 1, c2);
          if (_PB_N >= c1 + 3) {
            S_0(_PB_N - 1, c1, _PB_N - 1);
          } else
            for (int c4 = _PB_N - 2; c4 < _PB_N; c4 += 1)
              S_0(_PB_N - 1, c4, 2 * _PB_N - c4 - 3);
        }
      S_0(_PB_N - 1, _PB_N - 1, _PB_N - 1);
    }
  }
  S_0(_PB_N - 1, _PB_N - 2, _PB_N - 2);
} 
}

void OBST_04_Mod(int _PB_N, int **path)
{
 {
  if (_PB_N >= 3)
    S_0(0, 0, 0);
  if (_PB_N >= 2) {
    if (_PB_N == 2)
      S_0(0, 0, 0);
    S_0(_PB_N - 1, _PB_N - 2, _PB_N - 2);
  }
  for (int c0 = 0; c0 < 111 * _PB_N - 210; c0 += 1) {
    for (int c1 = max(0, floord(-10 * _PB_N + c0 + 9, 101) + 1); c1 <= min(_PB_N - 2, c0 / 101); c1 += 1) {
      if ((c0 - c1) % 10 == 0) {
        for (int c3 = c1 + 1; c3 < _PB_N; c3 += 1)
          S_0(c1, (c0 - 101 * c1) / 10, c3);
        if (c0 == 0 && c1 == 0) {
          for (int c2 = 1; c2 < _PB_N; c2 += 1)
            S_0(0, c2, 0);
        } else if (c1 >= 1 && 10 * _PB_N + 101 * c1 >= c0 + 110 && (c0 - c1) % 10 == 0) {
          S_0(c1, ((c0 - 101 * c1) / 10) + 10, c1);
          if (111 * c1 >= c0 + 110 && (91 * c0 - c1 + 10) % 100 == 0) {
            S_0(c1, c1, (c0 - 11 * c1 + 10) / 100);
          } else if (111 * c1 >= c0 + 100 && (91 * c0 - c1) % 100 == 0)
            for (int c2 = c1 + 1; c2 < _PB_N; c2 += 1)
              S_0(c1, c2, (c0 - 11 * c1) / 100);
        }
      }
      if (111 * c1 >= c0 + 11) {
        if ((91 * c0 - c1 + 91) % 100 == 0)
          for (int c2 = 0; c2 < c1 - 1; c2 += 1)
            S_0(c1, c2, (c0 - 11 * c1 + 1) / 100);
        if (c0 + 11 >= 10 * _PB_N + c1 && (90 * _PB_N + c0 - c1 + 11) % 100 == 0)
          S_0(c1, c1 - 1, (-10 * _PB_N + c0 - c1 + 11) / 100);
      } else if (c0 + 11 >= 10 * _PB_N + c1 && 10 * _PB_N + 101 * c1 >= c0 + 111 && (90 * _PB_N + c0 - c1 + 11) % 100 == 0)
        S_0(c1, c1 - 1, (-10 * _PB_N + c0 - c1 + 11) / 100);
    }
    if (c0 + 101 >= 101 * _PB_N) {
      if ((_PB_N - c0 + 9) % 10 == 0) {
        S_0(_PB_N - 1, (-101 * _PB_N + c0 + 201) / 10, _PB_N - 1);
        if ((11 * _PB_N - c0 + 79) % 100 == 0)
          S_0(_PB_N - 1, _PB_N - 1, (-11 * _PB_N + c0 + 21) / 100);
      }
      if ((11 * _PB_N - c0 + 88) % 100 == 0) {
        for (int c2 = 0; c2 < _PB_N - 2; c2 += 1)
          S_0(_PB_N - 1, c2, (-11 * _PB_N + c0 + 12) / 100);
        if (111 * _PB_N >= c0 + 312)
          S_0(_PB_N - 1, _PB_N - 2, (-11 * _PB_N + c0 + 12) / 100);
      }
    }
    for (int c1 = max(c0 / 101 + 1, floord(-10 * _PB_N + c0 + 8, 101) + 2); c1 <= min(_PB_N - 1, (c0 + 10) / 11); c1 += 1) {
      if (c0 + 100 >= 101 * c1 && (c0 - 11 * c1 + 100) % 100 >= 10 && (c0 - c1) % 10 == 0)
        S_0(c1, ((c0 - 101 * c1) / 10) + 10, c1);
      if ((91 * c0 - c1 + 10) % 100 == 0)
        S_0(c1, c1, (c0 - 11 * c1 + 10) / 100);
      if ((91 * c0 - c1 + 91) % 100 == 0)
        for (int c2 = 0; c2 < c1 - 1; c2 += 1)
          S_0(c1, c2, (c0 - 11 * c1 + 1) / 100);
      if (110 * _PB_N + c1 >= c0 + 311 && c0 + 11 >= 10 * _PB_N + c1 && (90 * _PB_N + c0 - c1 + 11) % 100 == 0)
        S_0(c1, c1 - 1, (-10 * _PB_N + c0 - c1 + 11) / 100);
      if ((91 * c0 - c1) % 100 == 0) {
        if (c0 + 100 >= 101 * c1)
          S_0(c1, ((c0 - 101 * c1) / 10) + 10, c1);
        for (int c2 = c1 + 1; c2 < _PB_N; c2 += 1)
          S_0(c1, c2, (c0 - 11 * c1) / 100);
      }
    }
  }
}   
}