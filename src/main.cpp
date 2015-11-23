#include <stdlib.h>
#include <iostream>

#include "../include/SimulationData.hpp"
#include "../include/SaveData.hpp"
#include "../include/WaveFunction.hpp"
#include "../include/Potential.hpp"
#include "../include/TridiagonalMatrices.hpp"
#include "../include/Solve.hpp"
#include "mkl.h"



int main(){
	
	//Clear any old FITS files to avoid exceptions thrown when
	//attempting to save
	system("rm fits/*.fits");

	//Set Environment variables for OpenMP
	//Sets length of time a spawned thread stays alive
//	putenv("KMP_BLOCKTIME=infinite");
	//Sets parameters for parallelization
//	putenv("KMP_AFFINITY=granularity=fine,compact,norespect");
	//Sets number of threads to use to maximum number available
//	mkl_set_num_threads(mkl_get_max_threads());
	//Disables mkl's fast memory management to reduce CPU overhead
	mkl_disable_fast_mm();

	//SimulationData class instance.
	//SimlationData SimulationData(int num_x, int num_y)
	SimulationData simData(256, 256);	
	printf("SimulationData Constructed\n");
	
	//Potential class instance
	Potential potentialData(simData);
	
	//WaveFunction class instance.
	WaveFunction psi(simData, potentialData.harmonicTrap);
	printf("WaveFunction Constructed\n");
	psi.getAbs(simData.getN());
	psi.getNorm(simData);

	potentialData.computeNonlinearEnergy(simData, psi);
	potentialData.assignTimeEvolutionOperator(simData, potentialData, false, true);
	printf("PotentialData Constructed\n");
	
	psi.getAbs(simData.getN());
	system("rm thomasFermi.fits");
	saveFITS(psi.absPsi, "thomasFermi.fits", simData); 

	//TridiagonalMatrix class instance
	TridiagonalMatrices matrices(simData);
	printf("TridiagonalMatrices Constructed\n");

	//Solve for ground state through Normalized Gradient Flow
	Solver solver(simData, psi, matrices, potentialData);
	printf("Solver Constructed\n");
	
	
	//Initialise Pardiso parameters
	solver.initialisePardiso(solver);
	//Analyse input matrices
	solver.analysePardiso(solver, matrices, psi);
	//LU factorize matrices
	solver.factorizePardiso(solver, matrices, psi);

	printf("Beginning Ground State Solver\n");
	for (simData.currStep = 0; simData.currStep < simData.numSteps; ++simData.currStep) {

		solver.solvePardiso(solver, matrices, psi, simData, potentialData, false);
	}
	psi.getAbs(simData.getN());
	system("rm groundState.fits");
	saveFITS(psi.absPsi, "groundState.fits", simData);	
	//Reassign Matrix Values for Real Computation
	matrices.reassignMatrixValues(simData);	

	solver.initialisePardiso(solver);
	solver.analysePardiso(solver, matrices, psi);
	solver.factorizePardiso(solver, matrices, psi);

	//Solve for dynamics
	for (simData.currStep = 0; simData.currStep < simData.numSteps; ++simData.currStep) {

		solver.solvePardiso(solver, matrices, psi, simData, potentialData, true);
	}

	solver.clearPardiso(solver, matrices, psi);	


	return 0;
}

