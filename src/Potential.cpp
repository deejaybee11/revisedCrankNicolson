#include "../include/Potential.hpp"

#include <stdlib.h>
#include <iostream>
#include <random>
#include <chrono>

#include "mkl.h"
#include "../include/SimulationData.hpp"

/**
 *This class constructs the arrays holding the potentials present in the GPE. The standard harmonic trap is present along with a
 *speckle potential used in simulations of Anderson localization. The speckle potential is modeled as an incident beam on a speckle
 *plate where a gaussian is assigned a random phase value at each point, and then passed through a bandpass filter, kFilter.
 *The resultant electric field then undergoes a Fourier transform and the absolute value/intensity is assigned to the array for the
 *potential. 
 */



Potential::Potential(SimulationData &simData) {
	//Allocate memory for potentials
	try {
		this->harmonicTrap = (double*)mkl_malloc(simData.getN() * sizeof(double), 64);
		this->specklePotential = (double*)mkl_malloc(simData.getN() * sizeof(double), 64);
		if (harmonicTrap == NULL) {
			throw -1;
		}
		if (specklePotential == NULL) {
			throw -2;
		}
	}
	catch (int e) {
		if (e == -1) {
			std::cout << "EXCEPTION " << e << " HARMONIC TRAP NOT ALLOCATED" << std::endl;
		}
		if ( e == -2) {
			std::cout << "EXCEPTION " << e << " SPECKLE POTENTIAL NOT ALLOCATED" << std::endl;
		}
	}

	//Declare and allocate memory for temporary arrays used in construction of the speckle potential
	double *phases = NULL;
	double *kFilter = NULL;
	double *vRandom = NULL;
	MKL_Complex16 *electricField = NULL;
	
	try {
		phases = (double*)mkl_malloc(simData.getN() * sizeof(double), 64);
		kFilter = (double*)mkl_malloc(simData.getN() * sizeof(double), 64);
		vRandom = (double*)mkl_malloc(simData.getN() * sizeof(double), 64);
		electricField = (MKL_Complex16*)mkl_malloc(simData.getN() * sizeof(MKL_Complex16), 64);
		if (phases == NULL) {
			throw -1;
		}
		if (kFilter == NULL) {
			throw -2;
		}
		if (electricField == NULL) {
			throw -3;
		}
		if (vRandom == NULL) {
			throw -4;
		}
	}
	catch (int e) {
		if (e == -1) {
			std::cout << "EXCEPTION " << e << " PHASES NOT ALLOCATED" << std::endl;
		}
		if (e == -2) {
			std::cout << "EXCEPTION " << e << " KFILTER NOT ALLOCATED" << std::endl;
		}
		if (e == -3) {
			std::cout << "EXCEPTION " << e << " ELECTRIC FIELD NOT ALLOCATED" << std::endl;
		}
		if (e == -4) {
			std::cout << "EXCEPTION " << e << " VRANDOM NOT ALLOCATED" << std::endl;
		}
	}


	//Construct the random seed for the speckle
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::exponential_distribution<double> distribtuion(1.0/200.0);
		
	//Calculate values relevant to speckle construction
	this->correlationLengthX = (1.0 / 1.5) * simData.healingLength;
	this->correlationLengthY = (1.0 / 1.5) * simData.healingLength;
	this->correlationLengthP_X = this->correlationLengthX * simData.getNumX() / simData.getLengthX();
	this->correlationLengthP_Y = this->correlationLengthY * simData.getNumY() / simData.getLengthY();

	this->correlationEnergy = pow(simData.hbar, 2.0) / (2.0 * simData.mass * pow(this->correlationLengthX, 2.0));
	this->specklePotentialStrength = 0.0325 * this->correlationEnergy;
	
	//Width of Gaussian for electricField
	double sig_x = 80, sig_y = 80;

	//Window in momentum space, acts as bandpass filter, ensures correlated numbers
	double kxCut = 1.0 / (2.0 * this->correlationLengthP_X);
	double kyCut = 1.0 / (2.0 * this->correlationLengthP_Y);

	//Time step is set based on the energy scale to ensure it is slowest process involved with speckle
	simData.set_dt(0.001 * simData.hbar / this->correlationEnergy);
		
	//Create exponentially correlated numbers
	for (int i = 0; i < simData.getN(); ++i) {
		phases[i] = fmod(distribtuion(generator), (2.0 * M_PI));
	}

	//Create k-space filter
	int index;
	double tempX, tempY;
	for (int i = 0; i < simData.getNumX(); ++i) {
		for (int j = 0; j < simData.getNumY(); ++j) {
			index = i * simData.getNumY() + j;

			tempX = simData.kX[i] * simData.getLengthX() / (double)simData.getN();
			tempY = simData.kY[i] * simData.getLengthY() / (double)simData.getN();

			if ((fabs(tempX < kxCut) && (fabs(tempY < kyCut)))) {
				kFilter[index] = 1;
			}
			else {
				kFilter[index] = 0;
			}
		}
	}

	//Construct electric field array
	double gaussian;
	for (int i = 0; i < simData.getNumX(); ++i) {
		for (int j = 0; j < simData.getNumY(); ++j) {
			index = i * simData.getNumY() + j;
			gaussian = exp(-1.0 * ((pow(simData.getNumX() / 2.0 + i, 2.0) / pow(2 * sig_x, 2.0)) + 
						(pow(simData.getNumY() / 2.0 + j, 2.0) / pow(2 * sig_y, 2.0))));
			electricField[index].real = gaussian * kFilter[index] * cos(phases[index]);
			electricField[index].imag = gaussian * kFilter[index] * sin(phases[index]);
		}
	}		

	//Fourier transform the electric field array using MKL DftiComputeForward
	//Handle for the fourier transform routine
	DFTI_DESCRIPTOR_HANDLE handle;
	//status variable
	MKL_LONG status = 0;
	//MKL_LONG N holds dimensions of array
	MKL_LONG N[2]; N[0] = simData.getNumX(); N[1] = simData.getNumY();
	//Create descriptor for type of transform
	status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, N);
	status = DftiCommitDescriptor(handle);
	//Compute the forward transform
	status = DftiComputeForward(handle, electricField, electricField);
	//Free descriptor from memory
	status = DftiFreeDescriptor(&handle);

	//absolute square of electric field taken
	for (int i = 0; i < simData.getN(); ++i) {
		vRandom[i] = pow((pow(electricField[i].real, 2.0) + pow(electricField[i].imag, 2.0)), 0.5);
	}

	//Compute the mean, and scale the potential from it
	double mean = 0;
	for (int i = 0; i < simData.getN(); ++i) {
		mean += vRandom[i] / simData.getN();
	}

	//Scale potential
	double temp = 0;
	for (int i = 0; i < simData.getN(); ++i) {
		temp = vRandom[i];
		vRandom[i] = this->specklePotentialStrength * (temp / abs(mean) - 1);
	}

	//Fill potential arrays with values
	for (int i = 0; i < simData.getNumX(); ++i) {
		for (int j = 0; j < simData.getNumY(); ++j) {
			index = i * simData.getNumY() + j;
			//harmonicTrapStrength takes on the value at point (i,j) of (m/2)*(omega_x^2 * x_i^2 + omega_y^2 * y_j^2) 
			//of a standard harmonic trap
			this->harmonicTrapStrength = (simData.mass / 2.0) * (pow(simData.omega_x, 2.0) * pow(simData.x[i], 2.0)
					+ pow(simData.omega_y, 2.0) * pow(simData.y[j], 2.0));
			this->harmonicTrap[index] = harmonicTrapStrength;

			//assign vRand to the speckle potential array
			this->specklePotential[index] = vRandom[index];

		}
	}

	//Free memory
	mkl_free(kFilter);
	mkl_free(electricField);
	mkl_free(phases);
	mkl_free(vRandom);
};


Potential::~Potential() {
	mkl_free(harmonicTrap);
	mkl_free(specklePotential);
	std::cout << "Potential Memory Cleared" << std::endl;
};



