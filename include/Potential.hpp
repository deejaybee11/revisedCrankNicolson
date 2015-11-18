#ifndef _POTENTIAL_H
#define _POTENTIAL_H

#include <stdlib.h>

#include "mkl.h"

#include "SimulationData.hpp"

class Potential {
public:
	Potential(SimulationData &simData);
	~Potential();

	double *harmonicTrap = NULL;			//Pointer for standard harmonic trap
	double *specklePotential = NULL;			//Pointer for a laser speckle potential
	double correlationLengthX;				//Correlation length of the random potential in X
	double correlationLengthY;				//Correlation length of the random potential in Y
	double correlationLengthP_X;				//Scaled correlation length
	double correlationLengthP_Y;				//Scaled correlation length
	double correlationEnergy;				//"energy" of the correlation

	
private:
	double harmonicTrapStrength;				//Temporary value used when filling harmonic trap array
	double specklePotentialStrength;			//Temporary value used when filling speckle potential array
};

#endif    //    _POTENTIAL_H
