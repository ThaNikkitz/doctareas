#include "statistics.h"

Statistics::Statistics()
{  
}

Statistics::~Statistics()
{
}

void Statistics::Histogram (double min, double max, int BinsNumber,std::vector<double> &array)
{
  gsl_histogram * h = gsl_histogram_alloc (BinsNumber);
  gsl_histogram_set_ranges_uniform (h, min, max);
  for (std::vector<double>::const_iterator it = array.begin() ; it != array.end(); ++it) 
  {
    gsl_histogram_increment (h, *it);
  }
  gsl_histogram_scale (h, (1/gsl_histogram_sum(h)));
  gsl_histogram_fprintf (stdout, h, "%15.5g", "%15.5g");
  gsl_histogram_free (h);
}

void Statistics::Histogram (double min, double max, int BinsNumber,std::vector<double> &array, std::string outname)
{
  FILE *out;
  out = fopen(outname.c_str(), "a");
  gsl_histogram * h = gsl_histogram_alloc (BinsNumber);
  gsl_histogram_set_ranges_uniform (h, min, max);
  for (std::vector<double>::const_iterator it = array.begin() ; it != array.end(); ++it) 
  {
    gsl_histogram_increment (h, *it);
  }
  gsl_histogram_scale (h, (1/gsl_histogram_sum(h)));
  gsl_histogram_fprintf (out, h, "%15.5g", "%15.5g");
  gsl_histogram_free (h);
  fclose(out);
}

void Statistics::Average(std::vector<double> &array)  
{
  double sum = 0.0;
  double average = 0.0;
  double sd = 0.0;  
  int NumElements = array.size();
  for (std::vector<double>::const_iterator it = array.begin() ; it != array.end(); ++it)
    sum +=*it;
  average = sum/NumElements;
  sum=0.0;
  for (std::vector<double>::const_iterator it = array.begin() ; it != array.end(); ++it)
    sum += (*it - average)*(*it - average);
  sd =sqrt(sum/NumElements);
  std::cout<<std::setw(10)<<"#Average ="<<std::setw(7)<<std::setprecision(3)<<average<<std::endl;
  std::cout<<std::setw(10)<<"#SD ="<<std::setw(7)<<std::setprecision(3)<<sd<<std::endl;
        
}

double Statistics::min(std::vector<double> &array) 
{
  double x,temp;
  x= array[0];
 for (std::vector<double>::const_iterator it = array.begin() ; it != array.end(); ++it)
 {
   temp=*it;
   if (temp<x)
     x=temp;
 }
 return x;
}

double Statistics::max(std::vector<double> &array) 
{
  double x,temp;
  x= array[0];
 for (std::vector<double>::const_iterator it = array.begin() ; it != array.end(); ++it)
 {
   temp=*it;
   if (temp>x)
     x=temp;
 }
 return x;
}
