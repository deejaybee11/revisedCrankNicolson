#ifndef CRANKNICOLSON_SIMULATION_DATA_H
#define CRANKNICOLSON_SIMULATION_DATA_H

#include <stdlib.h>

#include "mkl.h"

/**
 * Simulation Data Class
 *
 * Stores data and arrays used in the simulation
 */


class SimulationData {
public:
	SimulationData(int num_x, int num_y);
	~SimulationData();

	double *x;						//Pointer to x array
	double *y;						//Pointer to y array

	double *kX;						//Pointer to X wavenumber array
	double *kY;						//Pointer to Y wavenumber array

	double sigma_x;						//Width of Gaussian in X
	double sigma_y;						//Width of Gaussian in Y

	double omega_x;						//Trap frequency in X	
	double omega_y;						//Trap frequency in Y
	double omega_bar;					//Geometric mean of trap frequencies

	double mass;						//Mass of atomic species (Rb87)
	double U;						//Nonlinear interaction strength
	double scatteringLength;				//S-wave scattering length of Rb87
	double numAtoms;					//Number of atoms in BEC

	double chemicalPotential;				//Chemical potential mu
	double healingLength;					//healing length of condensate

	double hbar;						//Reduced plancks constant
	int numSteps;
	int printSteps;
	int currStep;
	int fileCount;
	const char *R = "REAL";
	const char *I = "IMAG";

	int getNumX() { return this->num_x; };			//Retrieves number of points in X
	int getNumY() { return this->num_y; };			//Retrieves number of points in Y
	int getN() { return this->N; };				//Retrieves total number of points

	double get_dx() { return this->dx; };			//Retrieves grid spacing in X
	double get_dy() { return this->dy; };			//Retrieves grid spacing in Y

	double getLengthX() { return this->length_x; };		//Retrieves physical length of X
	double getLengthY() { return this->length_y; };		//Retrieves physical length of Y

	void set_dt(double setdt);				//Sets the size of the time step dt
	double get_dt() { return this->dt; };			//Returns value of dt

private:
	int num_x;						//Number of points in X direction
	int num_y;						//Number of points in Y direction
	int N;							//Total number of grid points

	double length_x;					//Length of X in meters
	double length_y;					//Length of Y in meters
		
	double dx;						//Grain size in meters
	double dy;						//Grain size in meters

	double dt;						//Time step size. Time step is set in Potential class as a fraction of the correlation frequency
}; 





#endif    //    CRANKNICOLSON_SIMULATION_DATA_H
