#ifndef _WAVE_FUNCTION_H
#define _WAVE_FUNCTION_H

#include <stdlib.h>

#include "mkl.h"
#include "SimulationData.hpp"


class WaveFunction {
public:
	
	WaveFunction(SimulationData &simData, double *harmonicTrap);
	~WaveFunction();
	
	void getAbs(int N);				//Calculates absolute square of input array and saves to WaveFunction.absPsi
	void getNorm(SimulationData &simData);		//Calculates the normalization value N for normalizing via Psi = Psi / N

	MKL_Complex16 *psi = NULL;			//Pointer to wavefunction of condensate
	double *absPsi = NULL;				//Pointer to absolute square of condensate
	MKL_Complex16 *tempPsi = NULL;			//Pointer to temporary wavefunction array used in PARDISO

private:

	double normPsi;					//Normalization coefficient N as given in getNorm()
};

#endif    //    _WAVE_FUNCTION_H
