#include "fitsio.h"

#include <stdlib.h>
#include <iostream>

#include "../include/SimulationData.hpp"
#include "mkl.h"


void saveFITS(double *array, const char *filename, SimulationData &simData) {

	//Initialise a pointer an allocate memory for saving to FITS
	double *saveArray;
	try {
		saveArray = (double*)mkl_malloc(simData.getN() * sizeof(double), 64);
		if (saveArray == NULL) { throw -1; };
	}
	catch (int e) {
		if (e == -1) {
			std::cout << "EXCEPTION " << e << " MEMORY NOT ALLOCATED IN saveFITS" << std::endl;
		}
	}

	//Pointer to fits file
	fitsfile *fptr;

	//Status reports success of save routine
	int status = 0;

	//Values for number of axes, number of pixels
	long fpixel = 1, naxis = 2, nelements;
	long naxes[2] = {simData.getNumY(), simData.getNumX()};

	//Fills saveArray with input array data
	for (int i = 0; i < simData.getNumX()*simData.getNumY(); ++i) {
		saveArray[i] = array[i];
	}

	//Create fits file
	try {
		fits_create_file(&fptr, filename, &status);
		if (status != 0) {
			throw -1;
		}
		fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
		if (status != 0) {
			throw -2;
		}
		nelements = naxes[0] * naxes[1];
		fits_write_img(fptr, TDOUBLE, fpixel, nelements, saveArray, &status);
		if (status != 0) {
			throw -3;
		}
		fits_close_file(fptr, &status);
		if (status != 0) {
			throw -4;
		}
		fits_report_error(stderr, status);
	}
	catch (int e) {
		std::cout << "EXCEPTION " << e << " CAUGHT. STATUS = " << status << std::endl;
	}

	mkl_free(saveArray);
}


