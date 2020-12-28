#ifndef __INIT_H_INCLUDED__
#define __INIT_H_INCLUDED__

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
//#include "harmonic.h"

//class CoupledHarmonicSys;

//void initHarmonicOscill (CoupledHarmonicSys &rHarmonic);
void initNumberCycles (int  &NumberOfCycles);
void initNeq (int  &Neq, int NumberOfCycles);
std::string initOutHist();
std::string initOutEntropy();

#endif
