#ifndef _IDEAL_GAS_S_H_INCLUDED
#define _IDEAL_GAS_S_H_INCLUDED 

//#include "ran_uniform.h"
#include <time.h>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "../Eigen/Dense"
#include <boost/multi_array.hpp>


// 3D array typedefs
typedef boost::multi_array<int, 3> int3dArray;
typedef int3dArray::index array_index;

class IdealGasSys
{
 public:
  IdealGasSys();
  IdealGasSys(int & particles, int & maxEner, int &totEnergy);
  ~IdealGasSys();
  int GetNumPart() const;
  void SetNumPart(int & particles);
  int GetTotEnergy() const;
  int GetEnergySys() const;
  int GetTotMoves() const;
  int GetTotPos() const;
  int GetAcceptedMoves() const;
  Eigen::MatrixXi GetMatrixParticles ()const;
  int3dArray GetPositionHist ()const;
  void SetPositionHist(int & maxEner);
  void initParticles ();
  void initParticles (int & particles);
  void FillPart();
  void SamplePos();
  void SetTotEnergy(int & totEnergy);
  void SetEnergySys(int & particles);
  void SetEnergySys(double & particles);
  void CheckEnergiesInt() const;
  bool CheckEnergies() const ;
  void CalcEnergy() const;
  void cycle(int & Ncycles,int & equilibration, long int  seed = 0);
  void WriteFrame (std::ofstream &OutFrame, int &frame);
  void WriteEnergy(std::ofstream &OutFrame, int &frame);
  void WriteEneDist(std::ofstream & OutDist);
 private:
  int NumPart;
  int TotEnergy;
  int EnergySys;
  unsigned long  int TotMoves;
  unsigned long  int TotPos;
  unsigned long  int AcceptedMoves;
  // define a (Nx3) matrix with eigen
  Eigen::MatrixXi MatrixParticles;
  int3dArray PositionHist;
  //int3dArray PositionHist (extents[0][0][0]);
  
  
};

#endif
