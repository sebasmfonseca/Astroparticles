#include "Particle.h"

Particle::Particle(std::vector<double> x0,std::vector<double> d)
{
  if(x0.size()!=3 || d.size()!=3)
  {
    cout<<"Error"<<endl;
    exit(-1);
  }else
  {
    Position = x0;
    StartingPosition = x0;
    Direction = d;
    time =0;
  }


}

Particle::~Particle(){}

void Particle::PropagateSimple(double step)
{
  double NextHeight = StartingPosition[2]+2.9979e8*Direction[2]*step;

  while(NextHeight>0)
  {

    for(int i=0;i<3;i++)
    {
      Position[i] = Position[i]+Direction[i]*2.9979e8*step;
    }
    time+= step;
    NextHeight = Position[2]+2.9979e8*Direction[2]*step;
    
  }
  double auxstep = abs(Position[2]/(Direction[2]*2.9979e8));
  for(int i=0;i<3;i++)
  {
    Position[i] = Position[i]+Direction[i]*2.9979e8*auxstep;
  }
  time+= auxstep;


  return;
}
