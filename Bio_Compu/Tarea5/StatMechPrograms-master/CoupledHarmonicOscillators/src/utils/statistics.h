#ifndef __STATISTICS_H_INCLUDED__
#define __STATISTICS_H_INCLUDED__

#include <iostream>
#include <vector>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_histogram.h>
#include <stdio.h>
#include <stdlib.h>

class Statistics 
{
 public:
  Statistics ();
  ~Statistics();
  static void Histogram (double min, double max, int BinsNumber, std::vector<double> &array);
  static void Histogram (double min, double max, int BinsNumber, std::vector<double> &array, std::string outname);
  static void Average (std::vector<double> &array);
  static double min(std::vector<double> &array);
  static double max(std::vector<double> &array);
 private:
};

#endif

