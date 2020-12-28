#include "../utils/utils.h"
#include  <boost/range.hpp>
#include <iomanip>
#include "IG_sys.h"
#include <gsl/gsl_rng.h>

 IdealGasSys::IdealGasSys():
  NumPart (0),
  TotEnergy (0),
  EnergySys(0),
  TotMoves (0),
  TotPos (0),
  AcceptedMoves(0)
 {
 }

IdealGasSys::IdealGasSys(int &particles, int &maxEner, int &totEnergy):
  AcceptedMoves(0),
  TotPos (0),
  TotMoves(0),
  NumPart (particles),
  TotEnergy (totEnergy)
 {
    //Initialize the system;
    //SetNumPart(particles);
    //SetTotEnergy(totEnergy);
    SetEnergySys(particles);
    FillPart();
    //SetPositionHist(maxEner);
    SetPositionHist(maxEner);
 }

IdealGasSys::~IdealGasSys()
{
}


int IdealGasSys::GetNumPart() const
{
  return   NumPart;
}

void IdealGasSys::SetNumPart(int &particles)
{
  NumPart = particles;   
}

int IdealGasSys::GetTotEnergy() const
{
  return   TotEnergy;
}

int IdealGasSys::GetEnergySys() const
{
  return   EnergySys;
}

int IdealGasSys::GetTotMoves() const
{
  return   TotMoves;
}


int IdealGasSys::GetTotPos() const
{
  return   TotPos;
}


int IdealGasSys::GetAcceptedMoves() const
{
  return   AcceptedMoves;
}


void IdealGasSys::SetTotEnergy(int &totEnergy)
{
    TotEnergy = totEnergy;
}



void IdealGasSys::SetEnergySys(double &particles)
{
  int energy_int=(int)particles;
  EnergySys = 3*energy_int;
}

void IdealGasSys::SetEnergySys(int &particles)
{ 
  EnergySys = 3*particles;
}

Eigen::MatrixXi   IdealGasSys::GetMatrixParticles () const
{
  return MatrixParticles;
}

int3dArray IdealGasSys::GetPositionHist ()const
{
  
  return PositionHist;
}

void IdealGasSys::SetPositionHist(int &maxEner)
{
  int3dArray::extent_gen extents;
  PositionHist.resize(extents[maxEner][maxEner][maxEner]);
}



void IdealGasSys::FillPart()
{
  int Utot = 0;
  int i =0;
  MatrixParticles = (Eigen::MatrixXi::Constant(NumPart,3,0.0));
  do
     {
       if(i>=NumPart) 
	 i=0;
       MatrixParticles(i,0)++;
       MatrixParticles(i,1)++;
       MatrixParticles(i,2)++;
       i++;
       Utot+=3;
     }
      while(Utot!=EnergySys);
}

void IdealGasSys::SamplePos() 
{
  PositionHist[MatrixParticles(0,0)][MatrixParticles(0,1)][MatrixParticles(0,2)]++;
}

void IdealGasSys::CheckEnergiesInt() const
{
  int Utot = 0;  
  for(int i=0;i<NumPart;i++)
  {
    Utot+=MatrixParticles(i,0)*MatrixParticles(i,0);
    Utot+=MatrixParticles(i,1)*MatrixParticles(i,1);
    Utot+=MatrixParticles(i,2)*MatrixParticles(i,2);
  }
  if (Utot != EnergySys)
  {
    std::cout<<"Error,  energy system  internal =! energy system  energy"<<std::endl;
    exit(0);
  }
  else
    std::cout<<"\tInitial System energy  : "<<Utot<<std::endl; 
}

bool IdealGasSys::CheckEnergies() const
{
  if (TotEnergy<=3*NumPart)
    return true;
  else
    return false;
}

void IdealGasSys::CalcEnergy() const
{
  int Utot = 0;  
  for(int i=0;i<NumPart;i++)
  {
    Utot+=MatrixParticles(i,0)*MatrixParticles(i,0);
    Utot+=MatrixParticles(i,1)*MatrixParticles(i,1);
    Utot+=MatrixParticles(i,2)*MatrixParticles(i,2);
  }
  std::cout<<"\tSystem Energy          : "<<Utot<<std::endl; 
}


void IdealGasSys::cycle( int & Ncycles, int & equilibration, long int  seed )
{ 
  int particle,dimension,A,B,Uold,Unew,DeltaU;
  bool Lready;
  enum{PLUS=1,MINUS=-1};
  //InitializeRandomNumberGenerator(time(0l));
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r,seed);
  for(int i=0;i<NumPart*3;i++)
    
  {
    // Select Random Particle And Random X,Y,Z
    particle=NumPart*gsl_rng_uniform (r);
    
    dimension=3.0*gsl_rng_uniform (r);
    
    if(gsl_rng_uniform (r)< 0.5)
      A=PLUS;
    else
      A=MINUS;
        
    //Compute Energy Change
    Uold = MatrixParticles(particle,0)*MatrixParticles(particle,0) + 
           MatrixParticles(particle,1)*MatrixParticles(particle,1) + 
           MatrixParticles(particle,2)*MatrixParticles(particle,2);
    
    MatrixParticles(particle,dimension)+=A;
    
    Unew = MatrixParticles(particle,0)*MatrixParticles(particle,0) + 
           MatrixParticles(particle,1)*MatrixParticles(particle,1) + 
           MatrixParticles(particle,2)*MatrixParticles(particle,2);
    
    DeltaU = Unew - Uold;
    
  
    
    TotMoves++;
    // Accept Or Reject
    // Ground state is one
    if(MatrixParticles(particle,dimension) == 0)
    {
      Lready =false;
    }
    //System's energy  cannot be higher than the available energy 
    else if( (TotEnergy -(EnergySys + DeltaU) > 0))
    {
      Lready =true;
    }
    else 
    {
      Lready =false;
    }

    if(Lready) 
    {
      EnergySys += DeltaU;
      AcceptedMoves++;
    }
    else 
    {
	MatrixParticles(particle,dimension)-=A;
    }
  }
    if (Ncycles > equilibration )
    {
      TotPos++;
      if (max3<int>(MatrixParticles(0,0),MatrixParticles(0,1),MatrixParticles(0,2)) <= boost::size(PositionHist))
	SamplePos();
    }
    //}
  gsl_rng_free (r);
}	     
 	     
void IdealGasSys::WriteFrame (std::ofstream &OutFrame, int &frame) 
{
  OutFrame<<NumPart<<std::endl;
  OutFrame<<"Frame "<<frame<<std::endl;
    for(int i=0;i<NumPart;i++) 
  {
    OutFrame<<std::setw(10)<<"atom"<<i<<std::setw(10)<<MatrixParticles(i,0)<<std::setw(10)<<MatrixParticles(i,1)
	    <<std::setw(10)<<MatrixParticles(i,2)<<std::endl;
  }
}
  
void IdealGasSys::WriteEnergy(std::ofstream &Out, int &frame)
{
  
  Out<<std::setw(15)<<frame<<std::setw(15)<< MatrixParticles(0,0)*MatrixParticles(0,0) + 
     MatrixParticles(0,1)*MatrixParticles(0,1) + MatrixParticles(0,2)*MatrixParticles(0,2)
     <<std::setw(15)<<MatrixParticles(0,0)<<std::setw(15)<<MatrixParticles(0,1)
     <<std::setw(15)<<MatrixParticles(0,2)<<std::endl;

}


void IdealGasSys::WriteEneDist(std::ofstream &OutDist)
{
  for(int i=0;i<boost::size(PositionHist);i++)
    {
      for(int j=0;j<boost::size(PositionHist[0]);j++)
	{
	  for(int k=0;k<boost::size(PositionHist[0][0]);k++)
	    {
	      if (PositionHist[i][j][k] >= 0.5) 
		{
		  OutDist<<std::setw(10)<<i*i + j*j +k*k<<std::setw(10)<<PositionHist[i][j][k]<<std::endl;
		}
	    }
	}
    }
}
