#ifndef __generators__
#define __generators__

#include "tools.h"
#include "Particle.h"
using namespace std;

class Generator{
public:
  Generator();
  Generator(int);
  vector<double> Generate_Vector(); //creates a random, normalized 3d vector
  double Generate_Angle(); //Generate angle between -M_PI and M_PI
  vector<double> Generate_Direction_From_Theta(double theta); //creates a random, normalized 3d vector, given the angle with the z axis
  double Random_Distribution(double ,double ,TF1* ); //generates a random number according to a distribution
  vector<double> Generate_Position(double,double,double); //Generates a muon's starting position
  vector<double> Random_Distribution_2D(TF1*,double,double,double,double,double); //Generates 2 random numbers according to a pdf
  double Generate_Photon_Energy(); //Generates photon energy according to an interpolation of a pdf, use tools::Interpolate_Photon_Spectrum to get the distribution
  int Generate_Photon_Number(double); //Generates an integer according to a Poisson distribution around the input
  Particle* Generate_CosmicMuon(vector<double>);
  Particle* Generate_Photon(vector<double>);
  TF1* GetMomentumDistribution(){return Momentum_Distribution;};
  double Uniform(double xmin,double xmax){return Random->Uniform(xmin,xmax);};
  double Generate_Photon_Step();
  TSpline3* GetDetector_Efficiency(){return Detector_Efficiency;};
  ~Generator();

private:
  TRandom *Random;
  TF1 *Momentum_Distribution;
  TSpline3 *Photon_Spectrum;
  TSpline3 *Detector_Efficiency;
};

#endif
