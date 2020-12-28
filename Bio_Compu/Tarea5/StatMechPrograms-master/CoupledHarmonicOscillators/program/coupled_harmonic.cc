#include "../src/init.h"
#include "../src/harmonic.h"
#include  "../src/statistics.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include "../src/args/Arguments.h"
#include <iomanip> 
using namespace args;
using namespace std;

enum{NVE=1,NVT=2};
int debug_level = 0;

int main(int argc, char **argv) {
   Argument_List knowns;
   knowns << "oscillators" << "cycles"<<"equilibration"
	  <<"beta"<<"ensemble"<<"entropy"
	  <<"CalcEntropy"<<"TotEnergy";
  string usage = "# " + string(argv[0]);
  usage += "\n\t@oscillators   <number of oscillators>\n";
  usage += "\t@cycles   <number of cycles>\n";
  usage += "\t@TotEnergy   <Total initial energy of the system>\n";
  usage += "\t[@ensemble   <NVE=0 or NVT=1: default NVE>]\n";
  usage += "\t[@CalcEntropy   <calculate entropy : 0:no, 1:eq, 2:no eq, default 0>]\n";
  usage += "\t[@beta   <proportionality constant for NVT ensemble>]\n";
  usage += "\t@equilibtation   <equlibration steps>\n";
    
  try {
    Arguments args(argc, argv, knowns, usage);
    
    int oscillators;
    if (args.count("oscillators") > 0) {
      oscillators = args.getValue<int>("oscillators");
    } else {
      throw gromos::Exception("coupled_harmonic", "no oscillators specified (@oscillators)");
    }


    int TotEnergy;
    if (args.count("TotEnergy") > 0) {
      TotEnergy = args.getValue<int>("TotEnergy");
    } else {
      throw gromos::Exception("coupled_harmonic", "no oscillators specified (@TotEnergy)");
    }

  
    int cycles;
    if (args.count("cycles") > 0) {
      cycles = args.getValue<int>("cycles");
    } else {
      throw gromos::Exception("coupled_harmonic", "no cycles specified (@cycles)");
    }


    int equilibration;
    if (args.count("equilibration") > 0) {
      equilibration = args.getValue<int>("equilibration");
    } else {
      throw gromos::Exception("coupled_harmonic", "no equilibration steps specified (@equilibration)");
    }

    
    std::string ensemble;
    {
      bool tmp;
      if (args.count("ensemble") > 0) 
	tmp = args.getValue<int>("ensemble");
      if (!tmp)
	ensemble = "NVE";
      else
	ensemble = "NVT";
    }
    
    double beta=0.0;
    if (ensemble.compare("NVT") == 0) 
      if (args.count("beta") > 0) {
	beta = args.getValue<double>("beta");
      } else {
	throw gromos::Exception("coupled_harmonic", "NVT ensemble but no beta (@beta)");
      }
  

    int CalcEntropy =0; 
    if (args.count("CalcEntropy") > 0)
      CalcEntropy= args.getValue<int>("CalcEntropy");
  
     
     CoupledHarmonicSys Harmonic(oscillators,ensemble,CalcEntropy,TotEnergy);
    
    
    
    if (beta >0)
      Harmonic.SetBeta(beta);
    

    //set output variables
    std::string HistOutName="hist.dat";
    std::string EntropyOutName="entropy.dat";
    
    std::vector<double>Statistics;
    
    std::ofstream OutEntropy;
    if (Harmonic.GetEntropy() > 0) {
      OutEntropy.open(EntropyOutName.c_str());
      OutEntropy<<std::setw(15)<<"#Cycle"<<std::setw(15)<<"Entropy(log #)"<<std::endl;
    }
      
    // Check which ensemble
    int Choice;
    if ((Harmonic.GetEnsemble()).compare("NVE") == 0)
      Choice =NVE;
    if ((Harmonic.GetEnsemble()).compare("NVT") == 0)
      Choice =NVT;
 
  
    //loop over all cycles
    
    for(int i=0;i<cycles;i++)
      {
	switch(Choice)
	  {
	    // NVE ensemble
	    // choose 2 different levels at random
	    // exchange energy at random
	    // Compute entropy (NlogN - sum(nlogn))
	  case NVE:
	    Harmonic.NVEcycle();
	    if ((Harmonic.GetEntropy()) > 0) 
	      {
		//std::ofstream out(EntropyOutName.c_str() ,std::ios::app);
		OutEntropy<<std::setw(15)<<i<<std::setw(15)<<Harmonic.EntropyCalc()<<std::endl;
		//std::cout<<"Entropy:"<<Harmonic.EntropyCalc()<<std::endl;
	      }
	    break;
	    // NVT ensemble
	    // Boltzmann sampling at random
	  case NVT:
	    Harmonic.NVTcycle();
	    break;
	  }
	// sample first oscillator 
	if(i>equilibration)
	  Statistics.push_back((Harmonic.GetArrayOscillators())[0]);
      }

    if ((Harmonic.GetEnsemble()).compare("NVT") == 0) 
      Harmonic.SetTotEnergy((Harmonic.GetArrayOscillators())[0]);
    // Compute final results
    
    Harmonic.CheckEnergiesEnd();
    std::cout<<"First Oscillator"<< std::endl; 
    Statistics::Average(Statistics);
    std::cout<<"#Min Energy = "<<Statistics::min(Statistics)<< std::endl;
    std::cout<<"#Max Energy = "<<Statistics::max(Statistics)<< std::endl;
    
    std::ofstream OutHist(HistOutName.c_str());
    OutHist<<std::setw(15)<<"#Range"<<std::setw(15)<<""<<std::setw(15)<<"Probability"<<std::endl;
    OutHist.close();
    Statistics::Histogram(Statistics::min(Statistics),Statistics::max(Statistics),Statistics::max(Statistics),Statistics,HistOutName);
    std::cout<<"#Energy Histogram is stored in "<<HistOutName<< std::endl;
    if ((Harmonic.GetEntropy()) > 0)
      std::cout<<"#Entropy vs cycles is stored in "<<EntropyOutName<< std::endl;
    
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
    return 0;
}
