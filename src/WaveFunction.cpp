#include "../include/WaveFunction.hpp"

#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "mkl.h"
#include "../include/SimulationData.hpp"

//Constructor
WaveFunction::WaveFunction(SimulationData &simData, double *harmonicTrap) {

	//Allocate Memory
	try {
		this->psi = (MKL_Complex16*)mkl_malloc(simData.getN() * sizeof(MKL_Complex16), 64);
		this->absPsi = (double*)mkl_malloc(simData.getN() * sizeof(double), 64);
		this->tempPsi = (MKL_Complex16*)mkl_malloc(simData.getN() * sizeof(MKL_Complex16), 64);
		if (psi == NULL) {
			throw -1;
		}
		if (absPsi == NULL) {
			throw -2;
		}
		if (tempPsi == NULL) {
			throw -3;
		}
	}
	catch (int e) {
		if (e == -1) {
			std::cout << "EXCEPTION " << e << " PSI NOT ALLOCATED CORRECTLY" << std::endl;
		}
		if (e == -2) {
			std::cout << "EXCEPTION " << e << " ABSPSI NOT ALLOCATED CORRECTLY" << std::endl;
		}
		if (e == -3) {
			std::cout << "EXCEPTION " << e << " TEMPPSI NOT ALLOCATED CORRECTLY" << std::endl;
		}
	}
	//Fill arrays
	int index;
	for (int i = 0; i < simData.getNumX(); ++i) {
		for (int j = 0; j < simData.getNumY(); ++j) {
			index = simData.getNumY() * i + j;
	/*	
			this->psi[index].real = exp(-1.0 * (pow(simData.x[i], 2.0) / pow(simData.sigma_x, 2.0) + pow(simData.y[j], 2.0) / pow(simData.sigma_y, 2.0)));
			this->psi[index].imag = 0;

			this->absPsi[index] = 0;
		
			

			if(pow((pow(simData.x[i],2.0) + pow(simData.y[j],2.0)),0.5) <= 5e-6)
			{
				this->psi[index].real = 1;
			}
			else
			{
				this->psi[index].real=0;	
			}
		*/

			if ((simData.chemicalPotential - harmonicTrap[index]) < 0) {
				psi[index].real = 0;
				psi[index].imag = 0;
			}
			else {
				psi[index].real = (simData.chemicalPotential - harmonicTrap[index]) / (simData.numAtoms * simData.U);
				psi[index].imag = 0;
			}

			
		}
	}
}

//Destructor
WaveFunction::~WaveFunction() {

	mkl_free(psi);
	mkl_free(absPsi);
	mkl_free(tempPsi);
	std::cout << "WaveFunction Memory Cleared" << std::endl;
}

void WaveFunction::getNorm(SimulationData &simData) {

	double psiSum = 0;
	double tempReal;
	double tempImag;
	for (int i = 0; i < simData.getN(); ++i) {
		psiSum += this->absPsi[i];
	}

	this->normPsi = sqrt(1.0 / (psiSum * simData.get_dx() * simData.get_dy()));
	
	for (int i = 0; i < simData.getN(); ++i) {
		tempReal = this->psi[i].real * this->normPsi;
		tempImag = this->psi[i].imag * this->normPsi;
		this->psi[i].real = tempReal;
		this->psi[i].imag = tempImag;
	}
}

void WaveFunction::getAbs(int N) {
	vzAbs(N, this->psi, this->absPsi);
	vdMul(N, this->absPsi, this->absPsi, this->absPsi);
}
