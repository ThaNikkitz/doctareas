#include "init.h"
/*
void initHarmonicOscill (CoupledHarmonicSys &rHarmonic)  
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
  rHarmonic.SetNumOscil(NumberOfOscillators);
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
	    rHarmonic.SetEnsemble(ensemble);
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
	  rHarmonic.SetEntropy(entropy);
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
	  rHarmonic.SetTotEnergy(TotalEnergy);
	  //std::cin.ignore();
	  i++;
	}
      
      if (ensemble.compare("NVT") == 0)
	{
	  rHarmonic.SetEnsemble(ensemble);
	  rHarmonic.SetNumOscil(1);
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
	  
	  rHarmonic.SetBeta(beta);
	  i++;
	} 
      if ((ensemble.compare("NVT") != 0) && (ensemble.compare("NVE") != 0) )
	{
	  std::cout<<"Ensemble needs to be NVT or NVE"<<std::endl;
	  std::cout<<"Try again...."<<std::endl;
	  i=0;
	}	
    }
  std::cout<<"The system has "<<rHarmonic.GetNumOscil()<<" oscillators"<<std::endl;
  std::cout<<"Ensemble:  "<<rHarmonic.GetEnsemble()<<std::endl;
  if (ensemble.compare("NVE") == 0)
      std::cout<<"The total energy is  "<<rHarmonic.GetTotEnergy()<<std::endl;
  else
    std::cout<<"Beta is equal to "<<rHarmonic.GetBeta()<<std::endl;
  
}
*/
void initNumberCycles (int  &NumberOfCycles) 
{
  char buf[256];
  int i =0;
   
  while (i == 0)
    {
      std::cout<<"Number of cycles ?"<<std::endl;
      std::cin.getline(buf,256);
      //NumberOfCycles = atoi(buf);
       NumberOfCycles = strtol (buf,NULL,10);
      if (NumberOfCycles == 0) 
	{
	  std::cout << "Number of cycles needs to be a number > 0 :"<<std::endl;
	  std::cout<<"Try again...."<<std::endl;
	}
      else
	i++;
    }
}

void initNeq (int  &Neq, int NumberOfCycles) 
{
  char buf[256];
  int i =0;
   
  while (i == 0)
    {
      std::cout<<"Number of equilibration steps ?"<<std::endl;
      std::cin.getline(buf,256);
      //Neq = atoi(buf);
       Neq = strtol (buf,NULL,10);
      if ((Neq == 0) || (Neq >= NumberOfCycles )) 
	{
	  std::cout << "Number of equlibration steps needs to be a number and 0 > Neq < NumberOfCycles   :"<<std::endl;
	  std::cout<<"Try again...."<<std::endl;
	}
      else
	i++;
    }
}

std::string  initOutHist()
{
  int i =0;
  std::string outname;
  std::string replace;
  while (i == 0) 
  {
    std::cout<<"Histogram outname:"<<std::endl;
    std::cin>>outname;
    std::cin.ignore();
    std::ifstream in(outname.c_str());
    if (in) 
    {
      std::cout<<outname<<" already exists, overwrite it ?(yes/no)"<<std::endl;
      std::cin>>replace;
      std::cin.ignore();
      if (replace.compare("yes") != 0)
      {
	i=0;
	in.close();
      }
      else
      {
	i++;
	in.close();
	std::ofstream out(outname.c_str());
	out<<std::setw(15)<<"#Range"<<std::setw(15)<<""<<std::setw(15)<<"Probability"<<std::endl;
	out.close();
      }
    }
    else
    {
      i++;
      in.close();
      std::ofstream out(outname.c_str());
      out<<std::setw(15)<<"#Range"<<std::setw(15)<<""<<std::setw(15)<<"Probability"<<std::endl;
      out.close();
    }
    return outname;
  }
  
    
}

std::string  initOutEntropy()
{
  int i =0;
  std::string outname;
  std::string replace;
  while (i == 0) 
  {
    std::cout<<"Entropy vs cycles outname:"<<std::endl;
    std::cin>>outname;
    std::cin.ignore();
    std::ifstream in(outname.c_str());
    if (in) 
    {
      std::cout<<outname<<" already exists, overwrite it ?(yes/no)"<<std::endl;
      std::cin>>replace;
      std::cin.ignore();
      if (replace.compare("yes") != 0)
      {
	i=0;
	in.close();
      }
      else
      {
	i++;
	in.close();
	std::ofstream out(outname.c_str());
	out<<std::setw(15)<<"#Cycle"<<std::setw(15)<<"Entropy(log #)"<<std::endl;
	out.close();
      }
    }
    else
    {
      i++;
      in.close();
      std::ofstream out(outname.c_str());
      out<<std::setw(15)<<"#Cycle"<<std::setw(15)<<"Entropy(log #)"<<std::endl; 
      out.close();
    }
    return outname;
  }
  
    
}
