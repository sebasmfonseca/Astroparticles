#ifndef __Geometry__
#define __Geometry__

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include <iostream>

class Geometry{

public:

	// constructors and destructor
	Geometry();
	~Geometry();
	double GetDistance(){return Distance;};
	double GetHeight(){return Height;};
	double GetRadius(){return Radius;};
	TGeoManager* GetGeoManager(){return geom;};

	// Build Telescope methods
	void Build_MuonTelescope(double radius, double height, double distance,  double airgap, double althickness); // radius and height of scintillator and distance between scintillators

protected:
    double Radius;
	double Height;
	double Distance;
	double innerradius;
	double outerradius;
	double Thickness;
	double Airgap;
	double MaxHeight;
	double MaxRadius;


    TGeoManager* geom; //Manager pointer to geometry

};

#endif
