/*
 * refractiveMap.cpp
 *
 *  Created on: Apr 23, 2012
 *      Author: StaleyJW
 */


#include "refractiveMap.h"
#include "vector3D.h"
#include <boost/lexical_cast.hpp>





RefractiveMap::RefractiveMap()
: X_PML_OFFSET(25),
  Y_PML_OFFSET(10),
  Z_PML_OFFSET(10)
{
    refractive_x = refractive_y = refractive_z = NULL;
}


RefractiveMap::~RefractiveMap()
{

    /// Clean up memory.
    if (refractive_map)
    {
        delete refractive_map;
        refractive_map = NULL;
    }
}


//void RefractiveMap::loadRefractiveMap(void)
//{
//
//	// Ensure memory has been allocated for the pressure values that
//	// will be read in from file.  That is, initCommon() has already
//	// been called.
//	assert(refractive_grid != NULL);
//	// Verify that the calculations needed for the refractive index use values
//	// that have been set in the medium's call to this method.
//	assert(density > 0.0 &&
//			speed_of_sound > 0 &&
//			pezio_optical_coeff > 0);
//
//	// Open the file that contains the pressure values.
//	pressure_file_stream.open(pressure_file.c_str());
//
//	if (!pressure_file_stream) {
//		cout << "!!! Error opening pressure map file !!!\n";
//		exit(1);
//	}
//	else {
//		cout << "Pressure map " << pressure_file.c_str() << "...  opened successfully\n";
//		cout << "Loading PRESSURE values for calculating REFRACTIVE_INDEX values...\n";
//	}
//
//
//	double pressure_val = 0.0;
//		double M = 0.0;  // Coefficient of modulation.
//
//		for (array_index a = 0; a < Nx && pressure_file_stream.good(); a++)
//			for (array_index b = 0; b < Nz; b++)
//				for (array_index c = 0; c < Ny; c++)
//				{
//					pressure_file_stream >> pressure_val;
//					M = 2.0 * pezio_optical_coeff * pressure_val / (density * speed_of_sound * speed_of_sound);
//					(*refractive_grid)[a][b][c] = n_background * (1 + 0.5 * M);
//	#ifdef DEBUG
//					cout << (*refractive_grid)[a][b][c] << endl;
//	#endif
//				}
//
//	pressure_file_stream.close();
//}



//void RefractiveMap::loadRefractiveMap(const std::string &filename)
//{
//
//	// Ensure memory has been allocated for the pressure values that
//	// will be read in from file.  That is, initCommon() has already
//	// been called.
//	assert(refractive_grid != NULL);
//	// Verify that the calculations needed for the refractive index use values
//	// that have been set in the medium's call to this method.
//	assert(density > 0.0 &&
//			speed_of_sound > 0 &&
//			pezio_optical_coeff > 0);
//
//
//	std::string file_to_open = filename;
//	pressure_file_stream.open(file_to_open.c_str());
//
//	// Check for successful opening of the file.
//	if (!pressure_file_stream)
//	{
//		cout << "!!! Error opening pressure map file " << file_to_open.c_str() << "!!!\n";
//		exit(1);
//	}
//	else
//	{
//		cout << "Pressure map " << file_to_open.c_str() << " opened successfully. ";
//		cout << "Loading pressure values...\n";
//	}
//
//
//	//#define DEBUG
//
//
//	double pressure_val = 0.0;
//	double M = 0.0;  // Coefficient of modulation.
//
//	for (array_index a = 0; a < Nx && pressure_file_stream.good(); a++)
//		for (array_index b = 0; b < Nz; b++)
//			for (array_index c = 0; c < Ny; c++)
//			{
//				pressure_file_stream >> pressure_val;
//				M = 2.0 * pezio_optical_coeff * pressure_val / (density * speed_of_sound * speed_of_sound);
//				(*refractive_grid)[a][b][c] = n_background * (1 + 0.5 * M);
//#ifdef DEBUG
//				cout << (*refractive_grid)[a][b][c] << endl;
//#endif
//			}
//
//	pressure_file_stream.close();
//}


//void RefractiveMap::loadRefractiveMap(const std::string &filename, const int timeStep)
//{
//
//	// Ensure memory has been allocated for the pressure values that
//	// will be read in from file.  That is, initCommon() has already
//	// been called.
//	assert(refractive_grid != NULL);
//	// Verify that the calculations needed for the refractive index use values
//	// that have been set in the medium's call to this method.
//	assert(density > 0.0 &&
//			speed_of_sound > 0 &&
//			pezio_optical_coeff > 0);
//
//
//	// Concatonate the values passed in to form a filename to read in.
//	std::string file_to_open = filename + boost::lexical_cast<std::string>(timeStep);
//	pressure_file_stream.open(file_to_open.c_str());
//
//
//	// Check for successful opening of the file.
//	if (!pressure_file_stream)
//	{
//		cout << "!!! Error opening pressure map file " << file_to_open.c_str() << "!!!\n";
//		exit(1);
//	}
//	else
//	{
//		cout << "Pressure map " << file_to_open.c_str() << " opened successfully. ";
//		cout << "Loading pressure values...\n";
//	}
//
//
//	double pressure_val = 0.0;
//	double M = 0.0;  // Coefficient of modulation.
//
//	for (array_index a = 0; a < Nx && pressure_file_stream.good(); a++)
//		for (array_index b = 0; b < Nz; b++)
//			for (array_index c = 0; c < Ny; c++)
//			{
//				pressure_file_stream >> pressure_val;
//				M = 2.0 * pezio_optical_coeff * pressure_val / (density * speed_of_sound * speed_of_sound);
//				(*refractive_grid)[a][b][c] = n_background * (1 + 0.5 * M);
//#ifdef DEBUG
//				cout << (*refractive_grid)[a][b][c] << endl;
//#endif
//			}
//
//	pressure_file_stream.close();
//}




// Get the refractive index by converting the caartesian coords to voxel indeces.
//double RefractiveMap::getRefractiveIndex(const double x, const double y, const double z)
//{
//
//    cout << "ERROR: Stub function RefractiveMap::getRefractiveIndex(const double x, const double y, const double z)\n";
//    
//
//	// Indices into the displacement grids.
//	int _x, _y, _z;
//
//	boost::mutex::scoped_lock lock(m_refractive_mutex);
//	{
//		// Indices into the displacement grids.
//		_x = x/dx - (x/dx)/Nx;
//		_y = y/dy - (y/dy)/Ny;
//		_z = z/dz - (z/dz)/Nz;
//	}
//
//
//	return getRefractiveIndexFromGrid(_x, _y, _z);
//}


// Returns the individual axis displacement value from their location in the grid.
float
RefractiveMap::getRefractiveIndexFromGrid(const int x_photon, const int y_photon, const int z_photon)
{
    /// NEW: Used to access k-Wave data structure.
    return Get_refractive_index_TRealMatrix(x_photon + X_PML_OFFSET,
                                            y_photon + Y_PML_OFFSET,
                                            z_photon + Z_PML_OFFSET);
}


/// Return the refractive index from the TRealMatrix.
float
RefractiveMap::Get_refractive_index_TRealMatrix(const int x_photon, const int y_photon, const int z_photon)
{
    assert(refractive_map != NULL);
    
    /// NOTE: 'GetElementFrom3D' is a method of class 'TRealMatrix' provided by kWave.
    ///       To have the ultrasound z-axis orthogonal to the light z-axis we need to take
    ///       into account the way in which kWave stores data.  3-D grids are accessed as,
    ///       GetElementFrom3D(x-axis, y-axis, z-axis).  To accomodate this, Monte-carlo
    ///       x-coordinate is used to obtain kWave's z-coordinate data.  This means the
    ///       x-y plane in Monte-carlo forms the y-z plane in kWave.  Simply put, a point
    ///       in x-y from Monte-carlo retrieves data from z in kWave.
    ///       Below is Monte-carlo coordinates used to access kWave data so that they propagate
    ///       orthogonally.
    
    /// ON FURTHER INSPECTION I BELIEVE THE ABOVE IS NOT TRUE.  I THINK IT IS ALREADY MADE BY LOOKING
    /// AT THE TRANSDUCER PLOT IN KWAVE.  VERIFY!
	///                 Verified by Jiri (k-Wave developer)
    return refractive_map->GetElementFrom3D(x_photon, y_photon, z_photon);

}









