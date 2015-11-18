#include "../include/TridiagonalMatrices.hpp"

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "mkl.h"
#include "../include/SimulationData.hpp"
#include "../include/Potential.hpp"


TridiagonalMatrices::TridiagonalMatrices(SimulationData &simData, Potential &potentialData) {


	this->numMatrixElements = 3 * simData.getNumX();

	//Allocate memory for the matrices	
	try {
		this->aX = (MKL_Complex16*)mkl_malloc(this->numMatrixElements * sizeof(MKL_Complex16), 64);
		this->aY = (MKL_Complex16*)mkl_malloc(this->numMatrixElements * sizeof(MKL_Complex16), 64);
		this->bX = (MKL_Complex16*)mkl_malloc(this->numMatrixElements * sizeof(MKL_Complex16), 64);
		this->bY = (MKL_Complex16*)mkl_malloc(this->numMatrixElements * sizeof(MKL_Complex16), 64);

		if (this->aX == NULL) {
			throw -1;
		}
		if (this->aY == NULL) {
			throw -2;
		}
		if (this->bX == NULL) {
			throw -3;
		}
		if (this->bY == NULL) {
			throw -4;
		}
	}
	catch (int e) {
		if (e == -1) {
			std::cout << "EXCEPTION " << e << " A_X NOT ALLOCATED" << std::endl;
		}
		if (e == -2) {
			std::cout << "EXCEPTION " << e << " A_Y NOT ALLOCATED" << std::endl;
		}
		if (e == -3) {
			std::cout << "EXCEPTION " << e << "B_X NOT ALLOCATED" << std::endl;
		}
		if (e == -4) {
			std::cout << "EXCEPTION " << e << "B_Y NOT ALLOCATED" << std::endl;
		}
	}


	//Allocate memory for the column index arrays
	try {
		this->aColsX = (MKL_INT*)mkl_malloc(this->numMatrixElements * sizeof(MKL_INT), 64);
		this->aColsY = (MKL_INT*)mkl_malloc(this->numMatrixElements * sizeof(MKL_INT), 64);
		this->bColsX = (MKL_INT*)mkl_malloc(this->numMatrixElements * sizeof(MKL_INT), 64);
		this->bColsY = (MKL_INT*)mkl_malloc(this->numMatrixElements * sizeof(MKL_INT), 64);

		if (this->aColsX == NULL) {
			throw -1;
		}
		if (this->aColsY == NULL) {
			throw -2;
		}
		if (this->bColsX == NULL) {
			throw -3;
		}
		if (this->bColsY == NULL) {
			throw -4;
		}
	}
	catch (int e) {
		if (e == -1) {
			std::cout << "EXCEPTION " << e << " ACOLSX NOT ALLOCATED" << std::endl;
		}
		if (e == -2) {
			std::cout << "EXCEPTION " << e << " ACOLSY NOT ALLOCATED" << std::endl;
		}
		if (e == -3) {
			std::cout << "EXCEPTION " << e << "BCOLSX NOT ALLOCATED" << std::endl;
		}
		if (e == -4) {
			std::cout << "EXCEPTION " << e << "BCOLSY NOT ALLOCATED" << std::endl;
		}
	}

	//Allocate memory for the row index arrays
	try {
		this->aRowsX = (MKL_INT*)mkl_malloc((simData.getNumX() + 1) * sizeof(MKL_INT), 64);
		this->aRowsY = (MKL_INT*)mkl_malloc((simData.getNumY() + 1) * sizeof(MKL_INT), 64);
		this->bRowsX = (MKL_INT*)mkl_malloc((simData.getNumX() + 1) * sizeof(MKL_INT), 64);
		this->bRowsY = (MKL_INT*)mkl_malloc((simData.getNumY() + 1) * sizeof(MKL_INT), 64);

		if (this->aRowsX == NULL) {
			throw -1;
		}
		if (this->aRowsY == NULL) {
			throw -2;
		}
		if (this->bRowsX == NULL) {
			throw -3;
		}
		if (this->bRowsY == NULL) {
			throw -4;
		}
	}
	catch (int e) {
		if (e == -1) {
			std::cout << "EXCEPTION " << e << " AROWSX NOT ALLOCATED" << std::endl;
		}
		if (e == -2) {
			std::cout << "EXCEPTION " << e << " AROWSY NOT ALLOCATED" << std::endl;
		}
		if (e == -3) {
			std::cout << "EXCEPTION " << e << "BROWSX NOT ALLOCATED" << std::endl;
		}
		if (e == -4) {
			std::cout << "EXCEPTION " << e << "BROWSY NOT ALLOCATED" << std::endl;
		}
	}

	//Calculate required values. Currently assumes X and Y Dimensions are equal for derivativeCoefficient
	this->derivativeCoefficient = pow(simData.hbar, 2.0) / (2.0 * pow(simData.get_dx(), 2.0) * simData.mass);
	int index = 4;
	int columnIndex = 0;

	//Fill A,B and column array boundary values. 
	//First, second, second to last and last elements are assigned outside the loop as inside the loop, index - 1 segfaults for index == 0 and
	//index + 1 segfaults for index == N - 1
	//
	//simData.dt is set to -1j*dt initially in order to compute the ground state. Once computation is complete, reassignMatrixValues() is called
	//and simData.dt is set to be real again, and the matrices are repopulated.
	//This is implemented via the A and B matrices being purely real valued

	//Diagonal value
	this->aX[0].real = 1 - 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar; 
	this->aX[0].imag = 0;
	this->aY[0].real = 1 - 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar;
	this->aY[0].imag = 0; 

	this->bX[0].real = 1 + 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar;
	this->bX[0].imag = 0; 
	this->bY[0].real = 1 + 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar;
	this->bY[0].imag = 0; 

	//First off diagonal
	this->aX[1].imag = 0;
	this->aX[1].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->aY[1].imag = 0;
	this->aY[1].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	this->bX[1].imag = 0;
	this->bX[1].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->bY[1].imag = 0;
	this->bY[1].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	//Boundary value end of first row
	this->aX[2].imag = 0;
	this->aX[2].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->aY[2].imag = 0;
	this->aY[2].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	this->bX[2].imag = 0;
	this->bX[2].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->bY[2].imag = 0;
	this->bY[2].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;


	//Diagonal value
	this->aX[this->numMatrixElements - 1].imag = 0;
	this->aX[this->numMatrixElements - 1].real = 1 - 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[this->numMatrixElements - 1]) * simData.get_dt() / simData.hbar; 
	this->aY[this->numMatrixElements - 1].imag = 0;
	this->aY[this->numMatrixElements - 1].real = 1 - 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[this->numMatrixElements - 1]) * simData.get_dt() / simData.hbar;

	this->bX[this->numMatrixElements - 1].imag = 0;
	this->bX[this->numMatrixElements - 1].real = 1 + 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[this->numMatrixElements - 1]) * simData.get_dt() / simData.hbar;
	this->bY[this->numMatrixElements - 1].imag = 0;
	this->bY[this->numMatrixElements - 1].real = 1 + 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[this->numMatrixElements - 1]) * simData.get_dt() / simData.hbar;

	//First off diagonal
	this->aX[this->numMatrixElements - 2].imag = 0;
	this->aX[this->numMatrixElements - 2].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->aY[this->numMatrixElements - 2].imag = 0;
	this->aY[this->numMatrixElements - 2].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	this->bX[this->numMatrixElements - 2].imag = 0;
	this->bX[this->numMatrixElements - 2].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->bY[this->numMatrixElements - 2].imag = 0;
	this->bY[this->numMatrixElements - 2].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	//Boundary value end of first row
	this->aX[this->numMatrixElements - 3].imag = 0;
	this->aX[this->numMatrixElements - 3].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->aY[this->numMatrixElements - 3].imag = 0;
	this->aY[this->numMatrixElements - 3].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	this->bX[this->numMatrixElements - 3].imag = 0;
	this->bX[this->numMatrixElements - 3].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->bY[this->numMatrixElements - 3].imag = 0;
	this->bY[this->numMatrixElements - 3].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	//Fill First and last few column indices in A
	this->aColsX[0] = 0;
	this->aColsX[1] = 1;
	this->aColsX[2] = simData.getNumX() - 1;

	this->aColsY[0] = 0;
	this->aColsY[1] = 1;
	this->aColsY[2] = simData.getNumY() - 1;

	this->aColsX[this->numMatrixElements - 1] = simData.getNumX() - 1;
	this->aColsX[this->numMatrixElements - 2] = simData.getNumX() - 2;
	this->aColsX[this->numMatrixElements - 3] = simData.getNumX() - 3;

	this->aColsY[this->numMatrixElements - 1] = simData.getNumY() - 1;
	this->aColsY[this->numMatrixElements - 2] = simData.getNumY() - 2;
	this->aColsY[this->numMatrixElements - 3] = simData.getNumY() - 3;

	//Fill First and last few column indices in B
	this->bColsX[0] = 0;
	this->bColsX[1] = 1;
	this->bColsX[2] = simData.getNumX() - 1;

	this->bColsY[0] = 0;
	this->bColsY[1] = 1;
	this->bColsY[2] = simData.getNumY() - 1;

	this->bColsX[this->numMatrixElements - 1] = simData.getNumX() - 1;
	this->bColsX[this->numMatrixElements - 2] = simData.getNumX() - 2;
	this->bColsX[this->numMatrixElements - 3] = simData.getNumX() - 3;

	this->bColsY[this->numMatrixElements - 1] = simData.getNumY() - 1;
	this->bColsY[this->numMatrixElements - 2] = simData.getNumY() - 2;
	this->bColsY[this->numMatrixElements - 3] = simData.getNumY() - 3;

	//Iterate over remaining indices of arrays to fill them
	for(int i =1; i < simData.getNumX() - 1; ++i) {
		for(int j = 0; j < simData.getNumX(); ++j) {
			if (i == j) {

				try {
					if (index >= this->numMatrixElements) {
						throw -1;
					}
				}
				catch (int e) {
					std::cout << "ACCESSING OUTSIDE THE MATRIX ARRAYS" << std::endl;
				}

				//Diagonal elements assigned
				this->aX[index].imag = 0;
				this->aX[index].real = 1 - 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[i]) * simData.get_dt() / simData.hbar;
				this->aY[index].imag = 0;
				this->aY[index].real = 1 - 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[i]) * simData.get_dt() / simData.hbar;

				this->bX[index].imag = 0;
				this->bX[index].real = 1 + 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[i]) * simData.get_dt() / simData.hbar;
				this->bY[index].imag = 0;
				this->bY[index].real = 1 + 0.5 * (-2.0 * this->derivativeCoefficient + potentialData.harmonicTrap[i]) * simData.get_dt() / simData.hbar;

				//Upper diagonal elements assigned
				this->aX[index + 1].imag = 0;
				this->aX[index + 1].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
				this->aY[index + 1].imag = 0;
				this->aY[index + 1].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

				this->bX[index + 1].imag = 0;
				this->bX[index + 1].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
				this->bY[index + 1].imag = 0;
				this->bY[index + 1].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;				

				//Lower diagonal elements assigned
				this->aX[index - 1].imag = 0;
				this->aX[index - 1].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
				this->aY[index - 1].imag = 0;
				this->aY[index - 1].real = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

				this->bX[index - 1].imag = 0;
				this->bX[index - 1].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
				this->bY[index - 1].imag = 0;
				this->bY[index - 1].real = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;	

				//Assign column indices
				this->aColsX[index] = columnIndex + 1;
				this->aColsY[index] = columnIndex + 1;
				this->bColsX[index] = columnIndex + 1;
				this->bColsY[index] = columnIndex + 1;

				this->aColsX[index + 1] = columnIndex + 2;
				this->aColsY[index + 1] = columnIndex + 2;
				this->bColsX[index + 1] = columnIndex + 2;
				this->bColsY[index + 1] = columnIndex + 2;

				this->aColsX[index - 1] = columnIndex;
				this->aColsY[index - 1] = columnIndex;
				this->bColsX[index - 1] = columnIndex;
				this->bColsY[index - 1] = columnIndex;

				//Increment index and columIndex
				index += 3;
				columnIndex++;
			}
		}
	}





	//Fill arrays of row indices
	for (int i = 0; i < simData.getNumX() + 1; ++i) {
		try {
			if (index >= simData.getN()) {
				throw -1;
			}
		}
		catch (int e) {
			std::cout << "ACCESSING OUTSIDE THE ROWS ARRAYS" << std::endl;
		}
		this->aRowsX[i] = 3 * i;
		this->aRowsY[i] = 3 * i;
		this->bRowsX[i] = 3 * i;
		this->bRowsY[i] = 3 * i;

	}


};

//Destructor
TridiagonalMatrices::~TridiagonalMatrices() {
	mkl_free(aX);
	mkl_free(aY);
	mkl_free(bX);
	mkl_free(bY);

	mkl_free(aRowsX);
	mkl_free(aRowsY);
	mkl_free(bRowsX);
	mkl_free(bRowsY);

	mkl_free(aColsX);
	mkl_free(aColsY);
	mkl_free(bColsX);
	mkl_free(bColsY);

	std::cout << "Matrix Memory Freed" << std::endl;
};

//Reassign values for real time
void TridiagonalMatrices::reassignMatrixValues(SimulationData &simData, Potential &potentialData) {
	std::cout << "Reassigning Matrix Values For Real Time Evolution" << std::endl;

	int index = 4;
	//Diagonal value
	this->aX[0].real = 1;
	this->aX[0].imag = -0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar; 
	this->aY[0].real = 1;
	this->aY[0].imag = -0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar;

	this->bX[0].real = 1;
	this->bX[0].imag = 0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar;
	this->bY[0].real = 1;
	this->bY[0].imag = 0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar;

	//First off diagonal
	this->aX[1].real = 0;
	this->aX[1].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->aY[1].real = 0;
	this->aY[1].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	this->bX[1].real = 0;
	this->bX[1].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->bY[1].real = 0;
	this->bY[1].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	//Boundary value end of first row
	this->aX[2].real = 0;
	this->aX[2].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->aY[2].real = 0;
	this->aY[2].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	this->bX[2].real = 0;
	this->bX[2].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->bY[2].real = 0;
	this->bY[2].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;


	//Diagonal value
	this->aX[this->numMatrixElements - 1].real = 1;
	this->aX[this->numMatrixElements - 1].imag = -0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar; 
	this->aY[this->numMatrixElements - 1].real = 1;
	this->aY[this->numMatrixElements - 1].imag = -0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar;

	this->bX[this->numMatrixElements - 1].real = 1;
	this->bX[this->numMatrixElements - 1].imag = 0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar;
	this->bY[this->numMatrixElements - 1].real = 1;
	this->bY[this->numMatrixElements - 1].imag = 0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[0]) * simData.get_dt() / simData.hbar;

	//First off diagonal
	this->aX[this->numMatrixElements - 2].real = 0;
	this->aX[this->numMatrixElements - 2].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->aY[this->numMatrixElements - 2].real = 0;
	this->aY[this->numMatrixElements - 2].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	this->bX[this->numMatrixElements - 2].real = 0;
	this->bX[this->numMatrixElements - 2].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->bY[this->numMatrixElements - 2].real = 0;
	this->bY[this->numMatrixElements - 2].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	//Boundary value end of first row
	this->aX[this->numMatrixElements - 3].real = 0;
	this->aX[this->numMatrixElements - 3].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->aY[this->numMatrixElements - 3].real = 0;
	this->aY[this->numMatrixElements - 3].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

	this->bX[this->numMatrixElements - 3].real = 0;
	this->bX[this->numMatrixElements - 3].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
	this->bY[this->numMatrixElements - 3].real = 0;
	this->bY[this->numMatrixElements - 3].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;


	for(int i =1; i < simData.getNumX() - 1; ++i) {
		for(int j = 0; j < simData.getNumX(); ++j) {
			if (i == j) {

				try {
					if (index >= this->numMatrixElements) {
						throw -1;
					}
				}
				catch (int e) {
					std::cout << "ACCESSING OUTSIDE THE MATRIX ARRAYS" << std::endl;
				}

				//Diagonal elements assigned
				this->aX[index].real = 1;
				this->aX[index].imag = -0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[i]) * simData.get_dt() / simData.hbar;
				this->aY[index].real = 1;
				this->aY[index].imag = -0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[i]) * simData.get_dt() / simData.hbar;

				this->bX[index].real = 1;
				this->bX[index].imag = 0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[i]) * simData.get_dt() / simData.hbar;
				this->bY[index].real = 1;
				this->bY[index].imag = 0.5 * (-2.0 * this->derivativeCoefficient + 0*potentialData.harmonicTrap[i]) * simData.get_dt() / simData.hbar;

				//Upper diagonal elements assigned
				this->aX[index + 1].real = 0;
				this->aX[index + 1].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
				this->aY[index + 1].real = 0;
				this->aY[index + 1].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

				this->bX[index + 1].real = 0;
				this->bX[index + 1].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
				this->bY[index + 1].real = 0;
				this->bY[index + 1].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;				

				//Lower diagonal elements assigned
				this->aX[index - 1].real = 0;
				this->aX[index - 1].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
				this->aY[index - 1].real = 0;
				this->aY[index - 1].imag = -0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;

				this->bX[index - 1].real = 0;
				this->bX[index - 1].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;
				this->bY[index - 1].real = 0;
				this->bY[index - 1].imag = 0.5 * this->derivativeCoefficient * simData.get_dt() / simData.hbar;	



				//Increment index 
				index += 3;

			}
		}
	}
}

