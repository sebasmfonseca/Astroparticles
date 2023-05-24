#ifndef __Particle__
#define __Particle__

#include <vector>
#include <cmath>
#include <iostream>
#include "tools.h"

class Particle
{
public:
  Particle(std::vector<double>,std::vector<double>);
  ~Particle();

  void PropagateSimple(double); // Assumes Light speed;

  double GetTime(){return time;};
  double GetRadius(){return tools::Norm({Position[0],Position[1]});};
  vector<double> GetDirection(){return Direction;};

private:
  std::vector<double> StartingPosition;
  std::vector<double> Direction;
  std::vector<double> Position;
  double time;
};
#endif
