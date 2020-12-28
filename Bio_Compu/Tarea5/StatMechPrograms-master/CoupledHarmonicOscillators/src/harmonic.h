#ifndef _COUPLEDHARMONICSYS_H_INCLUDED_
#define _COUPLEDHARMONICSYS_H_INCLUDED_ 

//#include "ran_uniform.h"
#include "utils.h"
#include <time.h>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>

class CoupledHarmonicSys
{
 public:
  CoupledHarmonicSys();
  CoupledHarmonicSys(int oscillators, std::string ensemble,int entropy,int totalE);
  CoupledHarmonicSys(int oscillators, std::string ensemble, double beta);
  ~CoupledHarmonicSys();
  int GetNumOscil() const;
  void SetNumOscil(int oscillators);
  int GetTotEnergy() const;
  void SetTotEnergy(int energy);
  void SetTotEnergy(double energy);
  double GetBeta() const;
  void SetBeta(double beta);
  double GetEntropy() const;
  void SetEntropy(bool entropy);
  std::string GetEnsemble() const;
  void SetEnsemble (std::string ensemble );
  //std::vector<int> GetArrayOscillators () const;
  std::vector<double> GetArrayOscillators ()const ;
  void initHarmonicOscill ();
  void FillOcill();
  void CheckEnergiesInit() const;
  void CheckEnergiesEnd() const;
  void NVEcycle();
  void NVTcycle();
  double EntropyCalc();
 private:
  int NumOscil;
  int TotEnergy;
  double Beta;
  double Entropy;
  const gsl_rng_type * T;
  gsl_rng * r;
  std::string Ensemble;
  std::vector<int> ArrayOscillators;
  
};

#endif
