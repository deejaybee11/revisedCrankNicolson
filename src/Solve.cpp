#include "../include/Solve.hpp"

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <string.h>

#include "mkl.h"
#include "../include/SaveData.hpp"

//Constructor
Solver::Solver(SimulationData &simData, WaveFunction &psi, TridiagonalMatrices &matrices, Potential &potentialData) {

	//pt must be initialised with zeros before inputting into pardisoinit. iParameters[0] MUST be set to zero before passing
	//into pardisoinit as it sets all default parameters
	for (int i = 0; i < 64; ++i) {
		this->pt[i] = 0;
	}
	this->iParameters[0] = 0;
	//Set PARDISO parameters
	this->maxFactorizations = 1;
	this->mnum = 1;
	//matrixType = 13 for complex non-symmetric
	this->matrixType = 13;
	this->numRightHandSides = simData.getNumY();
	this->nEquations = simData.getNumX();
	//Print PARDISO output, useful for debug purposes
	this->messageLevel = 0;
	this->error = 0;
	this->matrixDescription[0] = 'G';
	this->matrixDescription[1] = 'L';
	this->matrixDescription[2] = 'N';
	this->matrixDescription[3] = 'C';

	//Allocate memory
	try {
		this->permutation = (MKL_INT*)mkl_malloc(simData.getNumX() * sizeof(MKL_INT), 64);
		this->alpha = (MKL_Complex16*)mkl_malloc(sizeof(MKL_Complex16), 64);
		this->beta = (MKL_Complex16*)mkl_malloc(sizeof(MKL_Complex16), 64);
		if (permutation == NULL) {
			throw -1;
		}
		if (alpha == NULL) {
			throw -2;
		}
		if (beta == NULL) {
			throw -3;
		}
	}
	catch (int e) {
		if (e == -1) {
			std::cout << "EXCEPTION " << e << " PERMUTATION NOT ALLOCATED" << std::endl;
		}
		if (e == -2) {
			std::cout << "EXCEPTION " << e << " APLHA NOT ALLOCATED" << std::endl;
		}
		if (e == -3) {
			std::cout << "EXCEPTION " << e << " BETA NOT ALLOCATED" << std::endl;
		}
	}
	
	//Set alpha to 1
	this->alpha[0].real = 1.0;
	this->alpha[0].imag = 0.0;
	//Set beta to 0
	this->beta[0].real = 0.0;
	this->beta[0].imag = 0.0;
};

Solver::~Solver() {
	mkl_free(permutation);
	std::cout << "Solver Memory Cleared" << std::endl;
};

void Solver::initialisePardiso(Solver &solver) {
	std::cout << "Calling pardisoinit" << std::endl;
	pardisoinit(this->pt, &matrixType, this->iParameters);
	//iparm[5] does not store the solution in psi on output
	this->iParameters[5] = 0;
	//PARDISO checks sparse matrices for errors
	this->iParameters[26] = 1;
	//Performs zero based indexing
	this->iParameters[34] = 1;
	//Set CSR format
	this->iParameters[36] = 0;
}

void Solver::analysePardiso(Solver &solver, TridiagonalMatrices &matrices, WaveFunction &psi) {
	std::cout << "Analysing matrices" << std::endl;
	//phase = 11 performs analysis of the matrices
	this->phase = 11;
	pardiso_64(this->pt, &this->maxFactorizations, &this->mnum, &this->matrixType, &this->phase, &this->nEquations, matrices.aX, matrices.aRowsX, matrices.aColsX,
			this->permutation, &this->numRightHandSides, this->iParameters, &this->messageLevel, psi.tempPsi, psi.psi, &this->error);
	if (this->error != 0) {
		std::cout << "Analysis Stage Failed With Error " << this->error << ". Please Refine" << std::endl;
	}
}

void Solver::factorizePardiso(Solver &solver, TridiagonalMatrices &matrices, WaveFunction &psi) {
	std::cout << "Factorising Matrices" << std::endl;
	//phase 22 performs numerical factorization
	this->phase = 22;
	pardiso_64(this->pt, &this->maxFactorizations, &this->mnum, &this->matrixType, &this->phase, &this->nEquations, matrices.aX, matrices.aRowsX, matrices.aColsX,
			this->permutation, &this->numRightHandSides, this->iParameters, &this->messageLevel, psi.tempPsi, psi.psi, &this->error);
	if (this->error != 0) {
		std::cout << "Factorization Stage Failed With Error " << this->error << ". Please Refine" << std::endl;
	}

	this->phase = 33;	
}

void Solver::solvePardiso(Solver &solver, TridiagonalMatrices &matrices, WaveFunction &psi, SimulationData &simData, Potential &potentialData, bool isReal) {
	if (simData.currStep%simData.printSteps == 0) {
		if (isReal) {
			printf("Real Step %d out of %d\n", simData.currStep, simData.numSteps);
			std::string filename = "fits/psi" + std::to_string(simData.fileCount) + ".fits";
			psi.getAbs(simData.getN());
			saveFITS(psi.absPsi, filename.c_str(), simData); 
			simData.fileCount++;
		}
		else {
			printf("Imaginary Step %d out of %d\n", simData.currStep, simData.numSteps);
		}
	}

	//Performs multiplication of time evolution operator of potentials
	potentialData.computeNonlinearEnergy(simData, psi);
	if (isReal) {
		potentialData.assignTimeEvolutionOperator(simData, potentialData, false, true);
	}
	else {
		potentialData.assignTimeEvolutionOperator(simData, potentialData, true, false);
	}
	vzMul(simData.getN(), psi.psi, potentialData.timeEvolutionOperator, psi.psi);	

	//Matrix multiplication of B*PSI in X direction
	mkl_zcsrmm(&this->N, &this->nEquations, &this->nEquations, &this->nEquations, this->alpha, this->matrixDescription, matrices.bX, matrices.bColsX,
		       matrices.bRowsX, matrices.bRowsX+1, psi.psi, &this->nEquations, this->beta, psi.tempPsi, &this->nEquations);	
	
	//Pardiso Solver Step in X direction
	pardiso_64(this->pt, &this->maxFactorizations, &this->mnum, &this->matrixType, &this->phase, &this->nEquations, matrices.aX, matrices.aRowsX, matrices.aColsX,
			this->permutation, &this->numRightHandSides, this->iParameters, &this->messageLevel, psi.tempPsi, psi.psi, &this->error);
	if (this->error != 0) {
		throw std::invalid_argument("PARDISO SOLVER X");
	}

	//Transpose Psi so above step is performed in Y direction
	mkl_zimatcopy('R', 'T', simData.getNumY(), simData.getNumX(), alpha[0], psi.psi, simData.getNumX(), simData.getNumY());

	//Matrix Multiplication of B*PSI in Y direction
	mkl_zcsrmm(&this->N, &this->nEquations, &this->nEquations, &this->nEquations, this->alpha, this->matrixDescription, matrices.bY, matrices.bColsY,
			matrices.bRowsY, matrices.bRowsY+1, psi.psi, &this->nEquations, this->beta, psi.tempPsi, &this->nEquations);

	//Pardiso Solver Step in Y direction
	pardiso_64(this->pt, &this->maxFactorizations, &this->mnum, &this->matrixType, &this->phase, &this->nEquations, matrices.aY, matrices.aRowsY, matrices.aColsY,
			this->permutation, &this->numRightHandSides, this->iParameters, &this->messageLevel, psi.tempPsi, psi.psi, &this->error);
	if (this->error != 0) {
		throw std::invalid_argument("PARDISO SOLVER Y");
	}

	//Transpose Psi back so next step is performed in X direction
	mkl_zimatcopy('R', 'T', simData.getNumX(), simData.getNumY(), alpha[0], psi.psi, simData.getNumY(), simData.getNumX());

	//Normalize Psi
	psi.getAbs(simData.getN());
	psi.getNorm(simData);


}

void Solver::clearPardiso(Solver &solver, TridiagonalMatrices &matrices, WaveFunction &psi) {
	this->phase = -1;
	pardiso_64(this->pt, &this->maxFactorizations, &this->mnum, &this->matrixType, &this->phase, &this->nEquations, matrices.aX, matrices.aRowsX, matrices.aColsX,
			this->permutation, &this->numRightHandSides, this->iParameters, &this->messageLevel, psi.tempPsi, psi.psi, &this->error);
	if (this->error != 0) {
		std::cout << "Pardiso Clear Failed With Error " << this->error << ". Please Refine" << std::endl;
	}
}
