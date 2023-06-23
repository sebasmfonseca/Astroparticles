#include "Generator.h"

Generator::Generator(double x0, double l,double x,double r) : TRandom(time(0)),X0(x0),lambda(l),Xmax(x),Rmax(r)
{}

Generator::~Generator(){}

double Generator::GaisserHillas(double X)
{
  double m = (Xmax-X0)/lambda;
  double x = (X-X0)/lambda;

  return pow((x/m),m)*exp(m-x);
}

double Generator::GenerateDepth()
{
  double X = Uniform(0,1030);
  double y = Uniform(1);
  while(y>GaisserHillas(X))
  {
    X = Uniform(0.2*Xmax,5*Xmax);
    y = Uniform(1);

  }
  return X;
}

vector<double> Generator::GenerateDirection()
{
  double theta = Uniform(M_PI/2);
  return {sin(theta),-cos(theta)};
}
