#ifndef __generators__
#define __generators__

#include "tools.h"
#include "Particle.h"

class Generator : public TRandom{

public:

Generator(int s = time(0));
~Generator();
std::vector<double> UniformDirection();
double NormalHeight(double x,double s){return Gaus(x,s);};

private:

};


#endif
