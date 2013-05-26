/*
 * refractiveMap.cpp
 *
 *  Created on: Apr 23, 2012
 *      Author: StaleyJW
 */


#include "refractiveMap.h"
#include "vector3D.h"
#include <boost/lexical_cast.hpp>




RefractiveMap::RefractiveMap(TRealMatrix *  pressure,
                             TRealMatrix *  rhox,
                             TRealMatrix *  rhoy,
                             TRealMatrix *  rhoz,
                             TRealMatrix *  rho0,
                             TRealMatrix *  c2,
                             const double   pezio_optical_coeff,
                             const double   n_background)
: X_PML_OFFSET(25),
  Y_PML_OFFSET(10),
  Z_PML_OFFSET(10)
{

    refractive_map = new TRealMatrix(pressure->GetDimensionSizes());

    this->pezio_optical_coeff   = pezio_optical_coeff;
    this->n_background          = n_background;

    //initCommon(voxel_dims);
    
    Update_refractive_map(pressure,
                          rhox,
                          rhoy,
                          rhoz,
                          rho0,
                          c2);
}


RefractiveMap::RefractiveMap(TRealMatrix *  pressure,
                             TRealMatrix *  rhox,
                             TRealMatrix *  rhoy,
                             TRealMatrix *  rhoz,
                             TRealMatrix *  rho0,
                             TRealMatrix *  c2,
                             size_t x_pml_offset,
                             size_t y_pml_offset,
                             size_t z_pml_offset,
                             const double   pezio_optical_coeff,
                             const double   n_background)
{
    
    
    refractive_map = new TRealMatrix(pressure->GetDimensionSizes());
    
    this->pezio_optical_coeff   = pezio_optical_coeff;
    this->n_background          = n_background;
    
    this->X_PML_OFFSET = x_pml_offset;
    this->Y_PML_OFFSET = y_pml_offset;
    this->Z_PML_OFFSET = z_pml_offset;
    
    Update_refractive_map(pressure,
                          rhox,
                          rhoy,
                          rhoz,
                          rho0,
                          c2);
    
}


//RefractiveMap::RefractiveMap(const std::string &filename, const int Nx, const int Nz, const int Ny, const int grid_size)
//{
//	// Assign the number of grid points (pixels in k-wave) used in the simulation.
//	this->Nx = Nx;
//	this->Ny = Ny;
//	this->Nz = Nz;
//
//	x_bound = y_bound = z_bound = grid_size;  // (meters)
//
//
//	// Initialize the data structures and values for the pressure map.
//	initCommon();
//
//
//	// FIXME: SHOULD USE BOOST FILESYSTEM TO GET initial_path AND DECIDE
//	//        WHAT FILE TO LOAD BASED ON CURRENT WORKING DIRECTORY.
//	// Assign the name to the pressure file.
//	pressure_file = filename;
//
//	// Load the pressure map values from disk file into pressure map array.
//	// NOTE: Order is important, this should be called after initCommon().
//	loadRefractiveMap();
//}


//RefractiveMap::RefractiveMap(const int Nx, const int Nz, const int Ny, const int grid_size)
//{
//	// Assign the number of grid points (pixels in k-wave) used in the simulation.
//	this->Nx = Nx;
//	this->Ny = Ny;
//	this->Nz = Nz;
//
//	// Sets the bounds of the pressure map grid.  Assumes uniform grid in each dimension.
//	x_bound = y_bound = z_bound = grid_size;  // (meters)
//
//	// Initialize the data structures and values for the pressure map.
//	initCommon();
//}

void RefractiveMap::initCommon(VoxelAttributes vox_attr)
{
//    voxel_dims.Nx = vox_attr.Nx;
//    voxel_dims.Ny = vox_attr.Ny;
//    voxel_dims.Nz = vox_attr.Nz;
//    
//    voxel_dims.dx = vox_attr.dx;
//    voxel_dims.dy = vox_attr.dy;
//    voxel_dims.dz = vox_attr.dz;
//    
//    
    
    
//	// Make sure the grid size (voxels in each axis) has been defined.
//	assert(Nx != 0 &&
//			Ny != 0 &&
//			Nz != 0);
//
//	dx = (double)x_bound / (double)Nx; // (meters)
//	dy = (double)y_bound / (double)Ny; // (meters)
//	dz = (double)z_bound / (double)Nz; // (meters)
//
//	// Set defaults for the refractive index grid attributes.
//	density = speed_of_sound = pezio_optical_coeff = n_background = 0.0;
//
//	refractive_grid = new three_dim_array (boost::extents[Nx][Nz][Ny]);
    
    
    
}


RefractiveMap::~RefractiveMap()
{
//	if (refractive_grid)
//	{
//		delete refractive_grid;
//		refractive_grid = NULL;
//	}
    
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


// Returns a Vector3d object holding values for displacements in all axes.
double RefractiveMap::getRefractiveIndex(const Vector3d &photonLocation)
{
    
    cout << "ERROR: Stub function RefractiveMap::getRefractiveIndex(const Vector3d &photonLocation)\n";
//	int _x, _y, _z;
//
//	boost::mutex::scoped_lock lock(m_refractive_mutex);
//	{
//		// Indices into the displacement grids.
//		_x = photonLocation.location.x/dx - (photonLocation.location.x/dx)/Nx;
//		_y = photonLocation.location.y/dy - (photonLocation.location.y/dy)/Ny;
//		_z = photonLocation.location.z/dz - (photonLocation.location.z/dz)/Nz;
//
//#ifdef DEBUG
//	// Sanity check.
//	assert(((_x < Nx && _x >= 0) &&
//			(_y < Ny && _y >= 0) &&
//			(_z < Nz && _z >= 0)) ||
//			assert_msg("_x=" << _x << " _y=" << _y << " _z=" << _z << "\n"
//					<< photonLocation.location.x << " "
//					<< photonLocation.location.y << " "
//					<< photonLocation.location.z));
//#endif
//	}
//
//
//	return getRefractiveIndexFromGrid(_x, _y, _z);
}

// Get the refractive index by converting the caartesian coords to voxel indeces.
double RefractiveMap::getRefractiveIndex(const double x, const double y, const double z)
{

    cout << "ERROR: Stub function RefractiveMap::getRefractiveIndex(const double x, const double y, const double z)\n";
    
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
}


// Returns the individual axis displacement value from their location in the grid.
float
RefractiveMap::getRefractiveIndexFromGrid(const int x_photon, const int y_photon, const int z_photon)
{
    /// OLD: Used to access boost 3D array.
	//return (*refractive_grid)[(array_index)a][(array_index)b][(array_index)c];
    
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
	/// Verified by Jiri (k-Wave developer)
    return refractive_map->GetElementFrom3D(x_photon, y_photon, z_photon);

}



/// Updates the current values in 'refractive_map' when a new pressure matrix is obtained
/// from a time step of the k-Wave simulation.
void
RefractiveMap::Update_refractive_map(TRealMatrix *pressure,
                                     TRealMatrix *rhox,
                                     TRealMatrix *rhoy,
                                     TRealMatrix *rhoz,
                                     TRealMatrix *rho0,
                                     TRealMatrix *c2)
{
    assert(pressure != NULL);
    assert(refractive_map != NULL);
    assert(pressure->GetTotalElementCount() == refractive_map->GetTotalElementCount());
    
    
    float M = 0.0;
    float density = 0.0;
    
    float * p_data      = pressure->GetRawData();
    float * rhox_data   = rhox->GetRawData();
    float * rhoy_data   = rhoy->GetRawData();
    float * rhoz_data   = rhoz->GetRawData();
    float * rho0_data   = rho0->GetRawData();
    float * c2_data     = c2->GetRawData();
    float * n_data      = refractive_map->GetRawData();
    
    
    /// Update the values of refraction in the TRealMatrix, based on the pressure passed in, at this timestep.
    for (int i = 0; i < refractive_map->GetTotalElementCount(); i++)
    {
        /// Calculate the background density with the addition of the pressure induced variations.
        /// NOTE:
        /// The density passed in to this function is density obtained from k-Wave.  In the description
        /// of how this density is calculated, described in k-wave_user_manual_1.0.1.pdf, it is not fully
        /// accurate do to the removal of the -u*grad(rho0) term in the mass conservation equation (Eq. 2.4).
        /// Three options exist:
        /// 1) Verify that the error is not significant and live with it.
        /// 2) Implement the term in k-Wave and pass it in here as an addition to rho (need each axial component).
        /// 3) Use a 1st order approximation from (p=c0^2*rho).
        ///density = rhox_data[i] + rhoy_data[i] + rhoz_data[i];       /// Density with error.
        ///density = p_data[i] / c2_data[i];                           /// 1st order approxmation.
        
        

        ///  Calculate the modulation coefficient as described by Skadazac and Wang.
	/// --------------------- THIS IS WRONG FOR PRESSURES I'M USING -------------------
        ///M = 2.0 * pezio_optical_coeff * (p_data[i] / (density * c2_data[i]));
        /// Update the refractive index value based the pressure induced changes.
        ///n_data[i] = n_background * (1 + 0.5 * M);
        


	/// Spatially varying pressure induced density changes
	density = sqrt(rhox_data[i]*rhox_data[i] +
		       rhoy_data[i]*rhoy_data[i] + 
		       rhoz_data[i]*rhoz_data[i]);
	/// "Optical Measurement of Ultrasonic Poynting and Velocity Vector Fields".  (Pitts, 2002)	
	/// Below uses the elasto-optical coefficient 
	n_data[i] = n_background + (n_background*n_background - 1)/(2*rho0_data[i]*n_background) *(density - rho0_data[i]);

    }

}










