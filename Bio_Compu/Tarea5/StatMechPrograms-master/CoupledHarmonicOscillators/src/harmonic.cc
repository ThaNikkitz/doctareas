#include "harmonic.h"

CoupledHarmonicSys::CoupledHarmonicSys():
  NumOscil (0),
  TotEnergy (0),
  Beta (0.0),
  Entropy(0.0),
  Ensemble ("NVE")
 {
   //InitializeRandomNumberGenerator(time(0l));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);		
 }


CoupledHarmonicSys::CoupledHarmonicSys(int oscillators, std::string ensemble,int entropy,int totalE):
  NumOscil (oscillators),
  Entropy(entropy),
  Ensemble (ensemble),
  TotEnergy(totalE)
 {
   std::cout<<"#Oscillators "<<NumOscil<<std::endl;
   std::cout<<"#Total Energy "<<TotEnergy<<std::endl;
   std::cout<<"#Ensemble "<<Ensemble<<std::endl;
   std::cout<<"#Entropy "<<Entropy<<std::endl;
   FillOcill();
   //initHarmonicOscill();
   CheckEnergiesInit();
   //InitializeRandomNumberGenerator(time(0l));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);		
 }

CoupledHarmonicSys::CoupledHarmonicSys(int oscillators, std::string ensemble, double beta):
  NumOscil (oscillators),
  Entropy(0.0),
  Beta(beta),
  Ensemble (ensemble)
 {
   FillOcill();
   //initHarmonicOscill();
   CheckEnergiesInit();
   //InitializeRandomNumberGenerator(time(0l));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);		
 }




CoupledHarmonicSys::~CoupledHarmonicSys()
{
}

int CoupledHarmonicSys::GetNumOscil() const
{
  return   NumOscil;
}

void CoupledHarmonicSys::SetNumOscil(int oscillators)
{
  NumOscil = oscillators;   
}

int CoupledHarmonicSys::GetTotEnergy() const
{
  return   TotEnergy;
}

void CoupledHarmonicSys::SetTotEnergy(int energy)
{
  int energy_int=(int)energy;
  TotEnergy = energy_int;
}

void CoupledHarmonicSys::SetTotEnergy(double energy)
{ 
  TotEnergy = energy;
}

double CoupledHarmonicSys::GetBeta() const
{
  return   Beta;
}

void CoupledHarmonicSys::SetBeta(double beta)
{
  Beta = beta;
}


double CoupledHarmonicSys::GetEntropy() const
{
  return   Entropy;
}

void CoupledHarmonicSys::SetEntropy(bool entropy)
{
  Entropy = entropy;
}


std::string CoupledHarmonicSys::GetEnsemble() const
{
  return   Ensemble;
}

void CoupledHarmonicSys::SetEnsemble(std::string ensemble)
{
  Ensemble = ensemble;
}

/*std::vector<int>  CoupledHarmonicSys::GetArrayOscillators () const
{
  return  ArrayOscillators;
}
*/


std::vector<double>  CoupledHarmonicSys::GetArrayOscillators () const
{
  std::vector<double> v_double(ArrayOscillators.begin(), ArrayOscillators.end()); 
 return v_double;
}

void CoupledHarmonicSys::FillOcill()
{
  int Utot = 0;
  int i =0;
  for (int j =0; j< NumOscil; j++)
    ArrayOscillators.push_back(0.0);
  if (Ensemble.compare("NVE")==0)
  {
    if (Entropy < 2)
    {
      do
      {
	if(i>=NumOscil) 
	  i=0;
	ArrayOscillators[i]++;
	i++;
	Utot++;
      }
      while(Utot!=TotEnergy);
    } else
      ArrayOscillators[0]=TotEnergy; 
  } 
}

void CoupledHarmonicSys::CheckEnergiesInit() const
{
  int Utot = 0;  
  for(int i=0;i<NumOscil;i++)
  {
    Utot+=ArrayOscillators[i];
  }
  if (Utot != TotEnergy)
  {
    std::cout<<"Error Initial energy =! Total energy"<<std::endl;
    exit(0);
  }
  else
    std::cout<<"\nInitial energy : "<<Utot<<std::endl; 
}
void CoupledHarmonicSys::CheckEnergiesEnd() const
{
  int Utot = 0;  
  for(int i=0;i<NumOscil;i++)
  {
    Utot+=ArrayOscillators[i];
  }
  if (Utot != TotEnergy)
    {
      std::cout<<"Error Final energy =! Total energy"<<std::endl;
      exit(0);
    }
  else
    std::cout<<"\nFinal energy : "<<Utot<<std::endl; 
}

void CoupledHarmonicSys::NVEcycle()
{ 
  int OscA,OscB,A,B;
 enum{UP=1,DOWN=-1};
  for(int i=0;i<NumOscil;i++)
  {
    do
    {
      OscA=NumOscil*gsl_rng_uniform (r);
      OscB=NumOscil*gsl_rng_uniform (r);
    } while(OscA==OscB);
    
    if(gsl_rng_uniform (r)< 0.5)
      A=UP,B=DOWN;
    else
      A=DOWN,B=UP;
    if(MIN(ArrayOscillators[OscA]+A,ArrayOscillators[OscB]+B)>=0)
    {
      ArrayOscillators[OscA]+=A;
      ArrayOscillators[OscB]+=B;
    }
  }
}

void CoupledHarmonicSys::NVTcycle()
{
  int OscA,A;
  enum{UP=1,DOWN=-1};
  for(int i=0;i<NumOscil;i++)
  { 
    if(gsl_rng_uniform (r)<0.5)
      A=UP;
    else
      A=DOWN;
    if((ArrayOscillators[0]+A>=0)&& gsl_rng_uniform (r)<exp(-Beta*A))
      ArrayOscillators[0]+=A;   
  }
}

double CoupledHarmonicSys::EntropyCalc()
{
  std::vector<int> OcuppNumberArray;
  double Nfact;
  double PitOcupp=0.0;
  Nfact=((NumOscil)*log(NumOscil)) -NumOscil;
  for(int i =0; i< TotEnergy; i++) { 
    OcuppNumberArray.push_back(0);
  }
  for (std::vector<int>::const_iterator it = ArrayOscillators.begin() ; it != ArrayOscillators.end(); ++it)
    OcuppNumberArray[*it]++;
  for (std::vector<int>::const_iterator it = OcuppNumberArray.begin() ; it != OcuppNumberArray.end(); ++it)
  {  
    if (*it != 0) 
      PitOcupp+= ((*it)*log(*it)-*it);
  }
  return (Nfact - PitOcupp);
}

void CoupledHarmonicSys::initHarmonicOscill ()  
{
  int NumberOfOscillators = 0;
  int entropy =0;
  char buf[256];
  std::string ensemble;
  int TotalEnergy = 0;
  double beta = 0.0;
  int i =0;
  int j =0;
  int k=0;
  
  while (i == 0)
    {
      std::cout<<"Number of oscillators?"<<std::endl;
      std::cin.getline(buf,256);
      //NumberOfOscillators = atoi(buf);
      NumberOfOscillators = strtol (buf,NULL,10);
      if (NumberOfOscillators == 0) 
	{
	  std::cout << "Number of oscillators needs to be a number > 0 :"<<std::endl;
	  std::cout<<"Try again...."<<std::endl;
	}
      else
	i++;
    }
  SetNumOscil(NumberOfOscillators);
  i=0;
  while (i == 0)
    {

      std::cout<<"Ensemble (NVE or NVT)?"<<std::endl;
      std::cin>>ensemble;
      std::cin.ignore();
      if (ensemble.compare("NVE") == 0)
	{
	  while (k==0)
	  {
	    SetEnsemble(ensemble);
	    std::cout<<"Calculate Entropy?"<<std::endl;
	    std::cout<<"(0)- No calculation "<<std::endl;
	    std::cout<<"(1)- From max calculation "<<std::endl;
	    std::cout<<"(>=2)- From min calculation "<<std::endl;
	    std::cin.getline(buf,256);
	    entropy  = strtol (buf,NULL,10);
	    if (entropy == 0)
	    {
	     std::cout << "Entropy needs to be a number > 0 :"<<std::endl;
	     std::cout<<"Try again...."<<std::endl;  
	    }
	     else
	       k++;
	  }
	  SetEntropy(entropy);
	    while (j == 0)
	      {
		std::cout<<"Total Energy?"<<std::endl;
		std::cin.getline(buf,256);
		TotalEnergy = atoi(buf);
		if (TotalEnergy == 0)
		  {
		    std::cout << "Total energy needs to be a number > 0 :"<<std::endl;
		    std::cout<<"Try again...."<<std::endl;
		  }
		else 
		  j++;
	      }
	  SetTotEnergy(TotalEnergy);
	  //std::cin.ignore();
	  i++;
	}
      
      if (ensemble.compare("NVT") == 0)
	{
	  SetEnsemble(ensemble);
	  SetNumOscil(1);
	  j=0;
	  while (j == 0)
	    {
	      std::cout<<"Beta (inverse T)?"<<std::endl;
	      std::cin.getline(buf,256);
	      //beta = atof(buf);
	      beta = strtod (buf,NULL);
	      if (beta == 0.0)
		{
		  std::cout << "Beta needs to be a float number > 0 :"<<std::endl;
		  std::cout<<"Try again...."<<std::endl;
		}
	      else 
		j++;
	    }
	  
	  SetBeta(beta);
	  i++;
	} 
      if ((ensemble.compare("NVT") != 0) && (ensemble.compare("NVE") != 0) )
	{
	  std::cout<<"Ensemble needs to be NVT or NVE"<<std::endl;
	  std::cout<<"Try again...."<<std::endl;
	  i=0;
	}	
    }
  std::cout<<"The system has "<<GetNumOscil()<<" oscillators"<<std::endl;
  std::cout<<"Ensemble:  "<<GetEnsemble()<<std::endl;
  if (ensemble.compare("NVE") == 0)
      std::cout<<"The total energy is  "<<GetTotEnergy()<<std::endl;
  else
    std::cout<<"Beta is equal to "<<GetBeta()<<std::endl;
  
}
