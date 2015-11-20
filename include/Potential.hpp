#ifndef _POTENTIAL_H
#define _POTENTIAL_H

#include <stdlib.h>

#include "mkl.h"

#include "SimulationData.hpp"
#include "WaveFunction.hpp"

class Potential {
public:
	Potential(SimulationData &simData);
	~Potential();

	double *harmonicTrap = NULL;				//Pointer for standard harmonic trap
	double *specklePotential = NULL;			//Pointer for a laser speckle potential
	double *nonLinearPotential = NULL;			//Pointer for the |psi|^2 term
	MKL_Complex16 *timeEvolutionOperator = NULL;		//Pointer to exp(-i * (V + |psi|^2) * dt / hbar)
	double correlationLengthX;				//Correlation length of the random potential in X
	double correlationLengthY;				//Correlation length of the random potential in Y
	double correlationLengthP_X;				//Scaled correlation length
	double correlationLengthP_Y;				//Scaled correlation length
	double correlationEnergy;				//"energy" of the correlation

	void computeNonlinearEnergy(SimulationData &simData, WaveFunction &psi);
	void assignTimeEvolutionOperator(SimulationData &simData, Potential &potentialData, bool trapOn);

	
private:
	double harmonicTrapStrength;				//Temporary value used when filling harmonic trap array
	double specklePotentialStrength;			//Temporary value used when filling speckle potential array
};

#endif    //    _POTENTIAL_H
