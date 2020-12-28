#include "utils.h"
int MIN(int x, int y)
{
  int z;
  z= x<y?x:y;
  return z;
}

int MAX(int x, int y)
{
  int z;
  z= x>y?x:y;
  return z;
}

long double factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


