#ifndef _TRIDIAGONAL_MATRICES_H
#define _TRIDIAGONAL_MATRICES_H

#include <stdlib.h>

#include "mkl.h"
#include "SimulationData.hpp"
#include "Potential.hpp"


class TridiagonalMatrices {
public:
	TridiagonalMatrices(SimulationData &simData, Potential &potentialData);
	~TridiagonalMatrices();
	
	void reassignMatrixValues(SimulationData &simData, Potential &potentialData);

	//Matrices A_t, B_t, t = x,y, in linear equation A_t * Psi = B_t * Psi
	MKL_Complex16 *aX = NULL;
	MKL_Complex16 *aY = NULL;
	MKL_Complex16 *bX = NULL;
	MKL_Complex16 *bY = NULL;

	//Arrays that hold the index at which each new row begins.
	MKL_INT *aRowsX = NULL;
	MKL_INT *aRowsY = NULL;
	MKL_INT *bRowsX = NULL;
	MKL_INT *bRowsY = NULL;

	//Arrays that hold the column index of each value in A, B
	MKL_INT *aColsX = NULL;
	MKL_INT *aColsY = NULL;
	MKL_INT *bColsX = NULL;
	MKL_INT *bColsY = NULL;
	
	double derivativeCoefficient;
private:
	int numMatrixElements;
};

#endif    //    _TRIDIAGONAL_MATRICES_H
