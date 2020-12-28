// args_Arguments.t.cc

#include <cassert>
#include <cstdlib>
#include "../src/args/Arguments.h"
#include <iostream>
#include "../src/system/IG_sys.h"
#include <iomanip> 

using namespace std;
using namespace args;

int debug_level = 0;

int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "particles" << "cycles" <<"MaxEnergy"<<"equilibration"
	 <<"TotEnergy"<<"RandomSeed";
  string usage = "# " + string(argv[0]);
  usage += "\n\t@particles   <number of particles>\n";
  usage += "\t@cycles   <number of cycles>\n";
  usage += "\t@MaxEnergy   <Maximum energy per dimension>\n";
  usage += "\t@equilibtation   <equlibration steps>\n";
  usage += "\t@TotEnergy <Total energy of the system + demon>\n";
  usage += "\t[@RandomSeed   <seed for random number generator, set it to get always the same numerical result>]\n";
  try {
    Arguments args(argc, argv, knowns, usage);
    int particles;
    if (args.count("particles") > 0) {
       particles = args.getValue<int>("particles");
    } else {
      throw gromos::Exception("IdealGas", "no particles specified (@particles)");
    }

    int cycles;
    if (args.count("cycles") > 0) {
      cycles = args.getValue<int>("cycles");
    } else {
      throw gromos::Exception("IdealGas", "no cycles specified (@cycles)");
    }

    int equilibration;
    if (args.count("equilibration") > 0) {
      equilibration = args.getValue<int>("equilibration");
    } else {
      throw gromos::Exception("IdealGas", "no equlibration steps specified (@equilibration)");
    }


    int MaxEnergy;
    if (args.count("MaxEnergy") > 0) {
      MaxEnergy = args.getValue<int>("MaxEnergy");
    } else {
      
      throw gromos::Exception("IdealGas", "no Max Energy specified (@MaxEnergy)");
    }
    
    int TotEnergy;
    if (args.count("TotEnergy") > 0) {
      TotEnergy = args.getValue<int>("TotEnergy");
    } else {
      
      throw gromos::Exception("IdealGas", "no Total Energy specified (@TotEnergy)");
    }
    

    
    long int seed;
    if (args.count("RandomSeed") > 0)
    {
      //forces an exact numerical result each time the code is run.
      seed = args.getValue<long int>("RandomSeed");
      if (seed < 0) 
      {
	throw gromos::Exception("IdealGas", "seed number lower than 0 (@RandomSeed)");
      } 
    }
      else 
    {
      // No exact numerical result each time the code is run.
      seed =-1;
      std::cout<<"\n\tWarning: random seed number."<<std::endl;
    }

    
    IdealGasSys system(particles,MaxEnergy,TotEnergy);
    if (system.CheckEnergies())
      throw gromos::Exception("IdealGas", "Total Energy lower or equal to system's energy ");
    
    std::cout<<"\n\tTotal Energy           : "<<system.GetTotEnergy()<<endl;
    std::cout<<"\tNumber of Particles    : "<<system.GetNumPart()<<endl;
    system.CheckEnergiesInt();
    //ad for loopwith the cycles
    ofstream et ( "et.dat" );
    //ofstream pt ("system_traj.xyz");
    et<<"#Energy trajectory for ideal quantum gas"<<std::endl; 
    et<<"#           Cycle         Energy            X              Y              Z"<<std::endl;
    et.precision(3);
    //pt.precision(3);
    if (seed < 0) 
    {
      srand((time(0))); // srand & time are built-in
      seed = random();
    }
    int TotCount =0; 
    for (int i=0;i<cycles;i++)
    {
      for (int j=0;j<1000;j++)
      {
	system.cycle(i,equilibration,seed);
	TotCount++;
	if (seed < 0)
	{
	  srand((time(0))); // srand & time are built-in
	  seed = random();  // set a random number form CPU internal clock
	// change seed number to randomize each cycle if no ramdom seed was occupied.
	// Forces an exact numerical result every time the code is run.
	} else
	  seed++;
	
	system.WriteEnergy(et,TotCount);
      }
      //system.WriteEnergy(et,i);
      //system.WriteFrame(pt,i);
    }
    
    et.close();
    //pt.close();
 
    
    ofstream ed("energy_distribution.dat");
    ed<<"#Energy distribution for ideal quantum gas"<<std::endl; 
    ed<<"#Energy       Count"<<std::endl;
    ed.precision(3);
    //std::cout<<equilibration<<std::endl;
    
    system.WriteEneDist(ed);
    ed.close();
    
    std::cout<<"\tAccepted moves         : "<<std::setprecision(3)
	     <<(float(system.GetAcceptedMoves())/float(system.GetTotMoves()))*100<<"%"<<std::endl;
    
    std::cout<<"\tFinal  Energy          : "<<system.GetEnergySys()<<std::endl;
    system.CalcEnergy();
    std::cout<<std::endl;
    
    

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
