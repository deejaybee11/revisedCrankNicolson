#ifndef _SOLVE_H
#define _SOLVE_H

#include <stdlib.h>

#include "mkl.h"
#include "Potential.hpp"
#include "SimulationData.hpp"
#include "TridiagonalMatrices.hpp"
#include "WaveFunction.hpp"

class Solver {
public:
	Solver(SimulationData &simData, WaveFunction &psi, TridiagonalMatrices &matrices, Potential &potentialData);
	~Solver();
	
	//Initialises pt, given matrixType and iParameters
	void initialisePardiso(Solver &solver);
	//Analyses matrices input into solver
	void analysePardiso(Solver &solver, TridiagonalMatrices &matrices, WaveFunction &psi);
	//Numerical factorization of matrices
	void factorizePardiso(Solver &solver, TridiagonalMatrices &matrices, WaveFunction &psi);
	//Solve equation with iterative refinement
	void solvePardiso(Solver &solver, TridiagonalMatrices &matrices, WaveFunction &psi, SimulationData &simData, Potential &potentialData, bool isReal);
	//Release internal memory
	void clearPardiso(Solver &solver, TridiagonalMatrices &matrices, WaveFunction &psi);


private:
	MKL_INT *pt[64];				//Internal data structure for PARDISO solver
	MKL_INT maxFactorizations;			//Maximum number of matrix factorizations to store in memory
	MKL_INT mnum;					//Which factorized matrix to use in solver
	MKL_INT matrixType;				//Defines matrix type, in this solver matrixType = 13, Complex and non-symmetric 
							//future work will be non-symmetric
	MKL_INT phase;
	MKL_INT solverPhase;				//Which phase of PARDISO to perform
	MKL_INT nEquations;				//Number of equations in linear system AX = B
	MKL_INT *permutation = NULL;			//Either holds permutation vector size n, or elements used in partial solution
	MKL_INT numRightHandSides;			//Number of right hand sides to solve for
	MKL_INT iParameters[64];			//Used to pass parameters to PARDISO solver.
	MKL_INT messageLevel;				//messageLevel = 0 generates no output, 1 generates output
	MKL_INT error;					//Error values, reference table on software.intel.com
	MKL_Complex16 *alpha;				//Scaling coefficient for zimatcopy for equation alpha*op(A), where op = transpose
	MKL_Complex16 *beta;				//Scaling coefficient for zcsrmm for equation C := A*X + beta*C
	char matrixDescription[6];			//Describes the type of matrix input into zcsrmm
	const char N = 'N';
};

#endif    //    _SOLVE_H
