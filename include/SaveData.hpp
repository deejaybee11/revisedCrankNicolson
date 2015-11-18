#ifndef _SAVE_FITS_DATA_H
#define _SAVE_FITS_DATA_H

#include "fitsio.h"

#include <stdlib.h>

#include "mkl.h"

//Saves a 2-dimensional array of doubles as a FITS file (Flexible Image Transport System) as developed by NASA.
//INPUT array MUST be of type double. FITS files can be opened in python with the astropy.io.fits library
void saveFITS(double *array, const char *filename, SimulationData &simData);


#endif    //    _SAVE_FITS_DATA_H

