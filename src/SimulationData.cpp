#include "../include/SimulationData.hpp"

#include <stdlib.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

//Constructor
SimulationData::SimulationData(int num_x, int num_y) {

	//Physical dimensions of grid
	this->num_x = num_x;
	this->num_y = num_y;
	this->N = num_x * num_y;
	this->length_x = 80e-6;
	this->length_y = 80e-6;

	//Iteration parameters
	this->dt = 0;
	this->numSteps = 30000;
	this->printSteps = 500;
	this->currStep = 0;
	this->fileCount = 0;

	//Time step is set based on the energy scale to ensure it is slowest process involved with speckle
	this->dt = 1e-8;
	std::cout << "DT = " << this->dt << std::endl;

	//Physical properties of Gaussian wavefunction
	this->sigma_x = 3e-6;
	this->sigma_y = 3e-6;

	//Harmonic trap properties
	this->omega_x = 2 * M_PI * 150;
	this->omega_y = 2 * M_PI * 150;

	//BEC parameters
	this->hbar = 1.054e-34;
	this->mass = 87 * 1.667e-27;
	this->scatteringLength = 5.186e-9;
	this->numAtoms = 30000.0;
	this->omega_bar = pow(this->omega_x * this->omega_x * this->omega_y, (1.0/3.0));
	this->U =(1.0 / pow(this->hbar / (this->omega_bar * this->mass), 0.5)) *  (4 * M_PI * pow(this->hbar, 2.0) * this->scatteringLength) / this->mass;
	this->chemicalPotential = (this->hbar * this->omega_bar / 2.0) * pow((15 * this->numAtoms * this->scatteringLength *
			       	(1.0 / pow(this->hbar / (this->omega_bar * this->mass), 0.5))), (2.0 / 5.0));
	this->healingLength = this->hbar / pow(2.0 * this->mass * this->chemicalPotential, 0.5);


	//Allocate memory
	try{

		this->x = (double*)mkl_malloc(this->num_x * sizeof(double), 64);
		if (this->x == NULL) { throw -1; };
		this->y = (double*)mkl_malloc(this->num_y * sizeof(double), 64);
		if (this->y == NULL) { throw -1; };
		this->kX = (double*)mkl_malloc(this->num_x * sizeof(double), 64);
		if (this->kX == NULL) { throw -1; };
		this->kY = (double*)mkl_malloc(this->num_y * sizeof(double), 64);
		if (this->kY == NULL) { throw -1; };

	}
	catch (int e){
		std::cout << "EXCEPTION CAUGHT" << std::endl;
		if(e == -1){
			std::cout << "EXCEPTION " << e << " MEMORY NOT ALLOCATED" << std::endl;
		}
	}

	//Fill X and Y arrays
	for (int i = 0; i < this->num_x; ++i) {
		this->x[i] = -1.0 * this->length_x / 2.0 + i * this->length_x / ((double)this->num_x);
	}
	for (int i = 0; i < this->num_y; ++i) {
		this->y[i] = -1.0 * this->length_y / 2.0 + i * this->length_y / ((double)this->num_y);
	}

	//Calculate step size
	this->dx = this->x[1] - this->x[0];	
	this->dy = this->y[1] - this->x[0];


	//Calculate increments for wavenumber arrays
	double aX = -1.0 * this->num_x / 2.0;
	double aY = -1.0 * this->num_y / 2.0;
	double bX = this->num_x / 2.0 - 1;
	double bY = this->num_y / 2.0 - 1;
	double stepX = (2 * M_PI / this->length_x) * ((bX - aX) / (this->num_x - 1.0));
	double stepY = (2 * M_PI / this->length_y) * ((bY - aY) / (this->num_y - 1.0));

	//Fill kX and kY arrays
	for (int i = 0; i < this->num_x; ++i) {
		this->kX[i] = (2 * M_PI / this->length_x) * aX + i * stepX;
	}
	for (int i = 0; i < this->num_y; ++i) {
		this->kY[i] = (2 * M_PI / this->length_y) * aY + i * stepY;
	}






};

//Destructor
SimulationData::~SimulationData() {
	mkl_free(x);
	mkl_free(y);
	mkl_free(kX);
	mkl_free(kY);
	std::cout << "SimulationData Memory Successfully cleared." << std::endl;
};

void SimulationData::set_dt(double setdt) {
     this->dt = setdt;
}     
