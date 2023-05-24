#include "Generator.h"

Generator::Generator(int s) : TRandom(s)
{

}

Generator::~Generator()
{

}

std::vector<double> Generator::UniformDirection()
{

  std::vector<double> v(3);
  v[0]= Uniform(-1,1);
  v[1]= Uniform(-1,1);
  v[2]= Uniform(-1,0);
  double norm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  for(int i=0;i<3;i++)
  {
    v[i]/=norm;
  }
  return v;
}
