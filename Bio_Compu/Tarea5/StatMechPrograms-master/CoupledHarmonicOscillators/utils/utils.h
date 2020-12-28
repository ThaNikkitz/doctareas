#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__
#include <math.h>
#include <algorithm>
// Some simple functions
int MIN(int x, int y);
int MAX(int x, int y);
long double factorial (int n);

template<typename T>
T const& max3(T const& a, T const& b, T const& c)
{
   using std::max;
   return max(max(a,b),c); // non-qualified max allows ADL
}


#endif
