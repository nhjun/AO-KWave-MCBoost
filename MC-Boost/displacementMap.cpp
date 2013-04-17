/*
 * displacementMap.cpp
 *
 *  Created on: 3 aug. 2011
 *      Author: StaleyJW
 */

#include "debug.h"
#include "vector3D.h"
#include "displacementMap.h"
#include <boost/lexical_cast.hpp>


// It's an error to create a DisplacementMap object without specifying attributes,
// therefore the default constructor should never be called.
DisplacementMap::DisplacementMap()
{
	cout << "Error: Default DisplacementMap() constructor called\n";
}


DisplacementMap::DisplacementMap(const std::string &filename, const int Nx, const int Nz, const int Ny, const int grid_size)
{
	// Assign the number of grid points (pixels in k-wave) used in the simulation.
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;

	x_bound = y_bound = z_bound = grid_size;  // (meters)


	initCommon();
}


DisplacementMap::DisplacementMap(const int Nx, const int Nz, const int Ny, const int grid_size)
{
	// Assign the number of grid points (pixels in k-wave) used in the simulation.
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;

	x_bound = y_bound = z_bound = grid_size;  // (meters)

	initCommon();
}



/// Displacement map containing the velocities in each axial direction for this time step.
DisplacementMap::DisplacementMap(TRealMatrix * velocity_x,
                                 TRealMatrix * velocity_y,
                                 TRealMatrix * velocity_z,
                                 const float US_freq,
                                 const float dt)
{
    /// Allocate displacement matrices for holding particle displacements used with MC sim.
    displacement_map_x = new TRealMatrix(velocity_x->GetDimensionSizes());
    displacement_map_y = new TRealMatrix(velocity_y->GetDimensionSizes());
    displacement_map_z = new TRealMatrix(velocity_z->GetDimensionSizes());
    
    /// Sanity check.
    assert (displacement_map_x != NULL);
    assert (displacement_map_y != NULL);
    assert (displacement_map_z != NULL);
    
    /// Displacements are updated per k-Wave time-step using a simple finite-difference scheme
    /// d = d + dt*u;  where u is the velocity.
    /// Since they are updated per time-step we need to zero them out.
    displacement_map_x->ZeroMatrix();
    displacement_map_y->ZeroMatrix();
    displacement_map_z->ZeroMatrix();
    
    
    
    velocity_map_x  = velocity_x;
    velocity_map_y  = velocity_y;
    velocity_map_z  = velocity_z;
    this->US_freq   = US_freq;
    this->dt        = dt;
    
    
    
    Update_displacement_map(velocity_x,
                            velocity_y,
                            velocity_z);
    
}

void DisplacementMap::initCommon()
{
//	assert(Nx != 0 &&
//			Ny != 0 &&
//			Nz != 0);
//
//	dx = (double)x_bound / (double)Nx; // (meters)
//	dy = (double)y_bound / (double)Ny; // (meters)
//	dz = (double)z_bound / (double)Nz; // (meters)
//	displacement_gridX = new three_dim_array (boost::extents[Nx][Nz][Ny]);
//	displacement_gridY = new three_dim_array (boost::extents[Nx][Nz][Ny]);
//	displacement_gridZ = new three_dim_array (boost::extents[Nx][Nz][Ny]);
}


DisplacementMap::~DisplacementMap()
{
	if (displacement_gridX)
		delete displacement_gridX;

	if (displacement_gridY)
		delete displacement_gridY;

	if (displacement_gridZ)
		delete displacement_gridZ;
    
    
    
    if (displacement_map_x)
    {
        delete displacement_map_x;
        displacement_map_x = NULL;
    }
    
    if (displacement_map_y)
    {
        delete displacement_map_y;
        displacement_map_y = NULL;
    }
    
    if (displacement_map_z)
    {
        delete displacement_map_z;
        displacement_map_z = NULL;
    }
    
}



/// XXX:
/// Old approach, no longer used.  No longer using offline calculations for particle displacements.  Currently data is
/// obtained from k-Wave (c++) at runtime.
// Loads a text files containing discrete displacement values at a given time step, in all dimensions (i.e. x, y, z),
// that were obtained from kWave simulation post-processed data.
void DisplacementMap::loadDisplacementMaps(const std::string &filename, const int timeStep)
{



//	// Assure memory has been allocated for the pressure values that
//	// will be read in from file.  That is, initCommon() has already
//	// been called.
//	assert(displacement_gridX != NULL);
//	assert(displacement_gridY != NULL);
//	assert(displacement_gridZ != NULL);

//	// A pointer to the string that is set for opening the displacement file.
//	std::string file_to_open;

//	// A pointer to one of the displacement grid arrays.  Set depending on
//	// which grid should be filled below.
//	three_dim_array *p_displacement_grid = NULL;

//	// Load each data structure with their respective displacement data.
//	for (int i = 0; i < 3; i++)
//	{


//		// Open the file that contains the pressure values for the specific
//		// dimension based on the loop index.
//		//
//		if (i == 0)
//		{  // X-displacement file.

//			// Clear the string.
//			file_to_open.clear();

//			// Concatonate the values passed in to form a filename to read in.
//			file_to_open = filename + "X-" + boost::lexical_cast<std::string>(timeStep) + ".txt";
//			disp_file_stream.open(file_to_open.c_str());

//			// The appropriate displacement grid is assigned to be filled below.
//			p_displacement_grid = displacement_gridX;
//		}
//		else if (i == 1)
//		{  // Y-displacement file.

//			// Clear the string.
//			file_to_open.clear();

//			// Concatonate the values passed in to form a filename to read in.
//			file_to_open = filename + "Y-" + boost::lexical_cast<std::string>(timeStep) + ".txt";
//			disp_file_stream.open(file_to_open.c_str());

//			// The appropriate displacement grid is assigned to be filled below.
//			p_displacement_grid = displacement_gridY;
//		}
//		else
//		{  // Z-displacement file.
//			file_to_open.clear();

//			// Concatonate the values passed in to form a filename to read in.
//			file_to_open = filename + "Z-" + boost::lexical_cast<std::string>(timeStep) + ".txt";
//			disp_file_stream.open(file_to_open.c_str());

//			// The appropriate displacement grid is assigned to be filled below.
//			p_displacement_grid = displacement_gridZ;
//		}


//		// Check for successful opening of the file.
//		if (!disp_file_stream)
//		{
//			cout << "!!! Error opening displacement map file " << file_to_open.c_str() << "!!!\n";
//			exit(1);
//		}
//		else
//		{
//			cout << "Displacement map " << file_to_open.c_str() << " opened successfully. ";
//			cout << "Loading displacement values...\n";
//		}


//		double data = 0.0;
//		// Read in data to the proper displacement array.
//		for (array_index a = 0; a < Nx && disp_file_stream.good(); a++)
//		{
//			for (array_index b = 0; b < Nz; b++)
//			{
//				for (array_index c = 0; c < Ny; c++)
//				{
//					disp_file_stream >> data;
//					(*p_displacement_grid)[a][b][c] = data;
//					//cout << (*p_displacement_grid)[a][b][c] << endl;
//				}
//			}
//		}


//		disp_file_stream.close();
//	}
}


/// XXX:
/// Old approach, no longer used.  No longer using offline calculations for particle displacements.  Currently data is
/// obtained from k-Wave (c++) at runtime.
void DisplacementMap::loadPressureAndCalculateDisplacements(const std::string &filename, const int timeStep,
															const double density,
															const double speed_of_sound,
															const double pezio_optic_coeff,
															const double background_refractive_index)
{

//	// Assure memory has been allocated for the pressure values that
//		// will be read in from file.  That is, initCommon() has already
//		// been called.
//		assert(displacement_gridX != NULL);
//		assert(displacement_gridY != NULL);
//		assert(displacement_gridZ != NULL);
//
//		// A pointer to the string that is set for opening the displacement file.
//		std::string file_to_open;
//
//		// A pointer to one of the displacement grid arrays.  Set depending on
//		// which grid should be filled below.
//		three_dim_array *p_displacement_grid = NULL;
//
//		// Load each data structure with their respective displacement data.
//		for (int i = 0; i < 3; i++)
//		{
//
//
//			// Open the file that contains the pressure values for the specific
//			// dimension based on the loop index.
//			//
//			if (i == 0)
//			{  // X-displacement file.
//
//				// Clear the string.
//				file_to_open.clear();
//
//				// Concatonate the values passed in to form a filename to read in.
//				file_to_open = filename + "X-" + boost::lexical_cast<std::string>(timeStep) + ".txt";
//				disp_file_stream.open(file_to_open.c_str());
//
//				// The appropriate displacement grid is assigned to be filled below.
//				p_displacement_grid = displacement_gridX;
//			}
//			else if (i == 1)
//			{  // Y-displacement file.
//
//				// Clear the string.
//				file_to_open.clear();
//
//				// Concatonate the values passed in to form a filename to read in.
//				file_to_open = filename + "Y-" + boost::lexical_cast<std::string>(timeStep) + ".txt";
//				disp_file_stream.open(file_to_open.c_str());
//
//				// The appropriate displacement grid is assigned to be filled below.
//				p_displacement_grid = displacement_gridY;
//			}
//			else
//			{  // Z-displacement file.
//				file_to_open.clear();
//
//				// Concatonate the values passed in to form a filename to read in.
//				file_to_open = filename + "Z-" + boost::lexical_cast<std::string>(timeStep) + ".txt";
//				disp_file_stream.open(file_to_open.c_str());
//
//				// The appropriate displacement grid is assigned to be filled below.
//				p_displacement_grid = displacement_gridZ;
//			}
//
//
//			// Check for successful opening of the file.
//			if (!disp_file_stream)
//			{
//				cout << "!!! Error opening displacement map file " << file_to_open.c_str() << "!!!\n";
//				exit(1);
//			}
//			else
//			{
//				cout << "Displacement map " << file_to_open.c_str() << " opened successfully. ";
//				cout << "Loading displacement values...\n";
//			}
//
//
//			//double pressure_val = 0.0;
//			// Read in data to the proper displacement array.
//			for (array_index a = 0; a < Nx && disp_file_stream.good(); a++)
//			{
//				for (array_index b = 0; b < Nz; b++)
//				{
//					for (array_index c = 0; c < Ny; c++)
//					{
//						
//                        cout << "DisplacementMap::loadPressureAndCalculateDisplacements NOT IMPLEMENTED!!!\n";
//						(*p_displacement_grid)[a][b][c] = -1;
//						//cout << (*p_displacement_grid)[a][b][c] << endl;
//					}
//				}
//			}
//
//
//			disp_file_stream.close();
//		}
}




// Returns a Vector3d object holding values for displacements in all axes.
boost::shared_ptr<Vector3d>
DisplacementMap::getDisplacements(const Vector3d &photonLocation)
{

//	boost::mutex::scoped_lock lock(m_displacement_mutex);

//	boost::shared_ptr<Vector3d> result(new Vector3d);


//	// Indices into the displacement grids.
//	int _x = photonLocation.location.x/dx - (photonLocation.location.x/dx)/Nx;
//	int _y = photonLocation.location.y/dy - (photonLocation.location.y/dy)/Ny;
//	int _z = photonLocation.location.z/dz - (photonLocation.location.z/dz)/Nz;


//	// Sanity check.
//	assert(((_x < Nx && _x >= 0) &&
//			(_y < Ny && _y >= 0) &&
//			(_z < Nz && _z >= 0)) ||
//			assert_msg("_x=" << _x << " _y=" << _y << " _z=" << _z << "\n"
//					<< photonLocation.location.x << " "
//					<< photonLocation.location.y << " "
//					<< photonLocation.location.z));




//	result->location.x = getDisplacementFromGridX(_x, _y, _z);
//	result->location.y = getDisplacementFromGridY(_x, _y, _z);
//	result->location.z = getDisplacementFromGridZ(_x, _y, _z);

//	return result;

}



boost::shared_ptr<Vector3d>
DisplacementMap::getDisplacements(const double x, const double y, const double z)
{


//	// Indices into the displacement grids.

//	boost::mutex::scoped_lock lock(m_displacement_mutex);

//	// Indices into the displacement grids.
//	int _x = x/dx - (x/dx)/Nx;
//	int _y = y/dy - (y/dy)/Ny;
//	int _z = z/dz - (z/dz)/Nz;

//	boost::shared_ptr<Vector3d> result (new Vector3d);
//	result->location.x = getDisplacementFromGridX(_x, _y, _z);
//	result->location.y = getDisplacementFromGridY(_x, _y, _z);
//	result->location.z = getDisplacementFromGridZ(_x, _y, _z);

//	return result;
}



void DisplacementMap::Update_displacement_map(TRealMatrix * velocity_x,
                                              TRealMatrix * velocity_y,
                                              TRealMatrix * velocity_z)
{
    assert (velocity_x != NULL);
    assert (velocity_y != NULL);
    assert (velocity_z != NULL);

    float *ux = velocity_x->GetRawData();
    float *uy = velocity_y->GetRawData();
    float *uz = velocity_z->GetRawData();
    
    float *d_map_x = displacement_map_x->GetRawData();
    float *d_map_y = displacement_map_y->GetRawData();
    float *d_map_z = displacement_map_z->GetRawData(); 
    
    size_t totalElements = velocity_x->GetTotalElementCount();
    
    for (size_t i = 0; i < totalElements; i++)
    {
        d_map_x[i] += ux[i]*dt;
        d_map_y[i] += uy[i]*dt;
        d_map_z[i] += uz[i]*dt;
    }
}


// Returns the individual axis displacement value from their location in the grid.
double DisplacementMap::getDisplacementFromGridX(const int x, const int y, const int z)
{
	//return (*displacement_gridX)[(array_index)a][(array_index)b][(array_index)c];
    
    return Get_displacement_X_TRealMatrix(x, y, z);
}

double DisplacementMap::getDisplacementFromGridY(const int x, const int y, const int z)
{
	//return (*displacement_gridY)[(array_index)a][(array_index)b][(array_index)c];
    
    return Get_displacement_Y_TRealMatrix(x, y, z);
}

double DisplacementMap::getDisplacementFromGridZ(const int x, const int y, const int z)
{
	//return (*displacement_gridZ)[(array_index)a][(array_index)b][(array_index)c];
    
    return Get_displacement_Z_TRealMatrix(x, y, z);
}
