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
    refractive_total = refractive_x = refractive_y = refractive_z = NULL;
}


RefractiveMap::~RefractiveMap()
{

    /// Clean up memory.
    if (refractive_total)
    {
        delete refractive_total;
        refractive_total = NULL;
    }
    
    if (refractive_x)
    {
        delete refractive_x;
        refractive_x = NULL;
    }
    
    if (refractive_y)
    {
        delete refractive_y;
        refractive_y = NULL;
    }
    
    if (refractive_z)
    {
        delete refractive_z;
        refractive_z = NULL;
    }
}



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






float
RefractiveMap::getRefractiveIndexFromGradientGrid(const char axis, const int x_photon, const int y_photon, const int z_photon)
{
    /// Refractive index changes induced from pressure propagating along x-axis.
    if (axis == 'x')
    {
        assert(refractive_x != NULL);
        return refractive_x->GetElementFrom3D(x_photon + X_PML_OFFSET,
                                       y_photon + Y_PML_OFFSET,
                                       z_photon + Z_PML_OFFSET);
    }

    /// Refractive index changes induced from pressure propagating along y-axis.
    if (axis == 'y')
    {
        assert(refractive_y != NULL);
        return refractive_y->GetElementFrom3D(x_photon + X_PML_OFFSET,
                                       y_photon + Y_PML_OFFSET,
                                       z_photon + Z_PML_OFFSET);
    }
    
    /// Refractive index changes induced from pressure propagating along z-axis.    
    if (axis == 'z')
    {
        assert(refractive_z != NULL);
        return refractive_z->GetElementFrom3D(x_photon + X_PML_OFFSET,
                                       y_photon + Y_PML_OFFSET,
                                       z_photon + Z_PML_OFFSET);
    }
    
    /// Error
    return -1.0;
}



// Returns the spatially located refractive index value from the supplied location in the grid.
/// x_photon, y_photon and z_photon are the voxel coordinates obtained by translating the spatial coordinate
/// of the current scattering event.
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
    assert(refractive_total != NULL);
    
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
    ///
	///              !!!!   Verified by Jiri (k-Wave developer) !!!!
    return refractive_total->GetElementFrom3D(x_photon, y_photon, z_photon);

}


/// Invert the phase of the refractive index data 180 degrees by multiplying through the matrix by -1.
void
RefractiveMap::Invert_phase(void)
{
    float * raw_data        = refractive_total->GetRawData();
    const size_t sensor_size  = refractive_total->GetTotalElementCount();

    for (size_t i = 0; i < sensor_size; i++)
    {
        /// Perform the inversion.
        raw_data[i] *= -1.0f;
    }
}








