#ifndef __generators__
#define __generators__

#include "TF1.h"
#include "tools.h"

class Generator : public TRandom{

public:

Generator(double,double,double,double);
~Generator();

double GaisserHillas(double);
double GenerateDepth();
double NormalHeight(double Xmax){return Gaus(Xmax,sqrt(Xmax));};
vector<double> GenerateDirection();


private:
  double X0;
  double lambda;
  double Xmax;
  double Rmax;
};


#endif
