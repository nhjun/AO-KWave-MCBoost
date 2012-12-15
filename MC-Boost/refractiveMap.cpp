/*
 * refractiveMap.cpp
 *
 *  Created on: Apr 23, 2012
 *      Author: StaleyJW
 */


#include "refractiveMap.h"
#include "vector3D.h"
#include <boost/lexical_cast.hpp>



// It's an error to create a PressureMap object without specifying attributes,
// therefore the default constructor should never be called.
RefractiveMap::RefractiveMap()
{
	cout << "Error: PressureMap() default constructor called\n";
}


RefractiveMap::RefractiveMap(const std::string &filename, const int Nx, const int Nz, const int Ny, const int grid_size)
{
	// Assign the number of grid points (pixels in k-wave) used in the simulation.
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;

	x_bound = y_bound = z_bound = grid_size;  // (meters)


	// Initialize the data structures and values for the pressure map.
	initCommon();


	// FIXME: SHOULD USE BOOST FILESYSTEM TO GET initial_path AND DECIDE
	//        WHAT FILE TO LOAD BASED ON CURRENT WORKING DIRECTORY.
	// Assign the name to the pressure file.
	pressure_file = filename;

	// Load the pressure map values from disk file into pressure map array.
	// NOTE: Order is important, this should be called after initCommon().
	loadRefractiveMap();
}


RefractiveMap::RefractiveMap(const int Nx, const int Nz, const int Ny, const int grid_size)
{
	// Assign the number of grid points (pixels in k-wave) used in the simulation.
	this->Nx = Nx;
	this->Ny = Ny;
	this->Nz = Nz;

	// Sets the bounds of the pressure map grid.  Assumes uniform grid in each dimension.
	x_bound = y_bound = z_bound = grid_size;  // (meters)

	// Initialize the data structures and values for the pressure map.
	initCommon();
}

void RefractiveMap::initCommon(void)
{
	// Make sure the grid size (voxels in each axis) has been defined.
	assert(Nx != 0 &&
			Ny != 0 &&
			Nz != 0);

	dx = (double)x_bound / (double)Nx; // (meters)
	dy = (double)y_bound / (double)Ny; // (meters)
	dz = (double)z_bound / (double)Nz; // (meters)

	// Set defaults for the refractive index grid attributes.
	density = speed_of_sound = pezio_optical_coeff = n_background = 0.0;

	refractive_grid = new three_dim_array (boost::extents[Nx][Nz][Ny]);
}


RefractiveMap::~RefractiveMap()
{
	if (refractive_grid)
	{
		delete refractive_grid;
		refractive_grid = NULL;
	}
}


void RefractiveMap::loadRefractiveMap(void)
{

	// Ensure memory has been allocated for the pressure values that
	// will be read in from file.  That is, initCommon() has already
	// been called.
	assert(refractive_grid != NULL);
	// Verify that the calculations needed for the refractive index use values
	// that have been set in the medium's call to this method.
	assert(density > 0.0 &&
			speed_of_sound > 0 &&
			pezio_optical_coeff > 0);

	// Open the file that contains the pressure values.
	pressure_file_stream.open(pressure_file.c_str());

	if (!pressure_file_stream) {
		cout << "!!! Error opening pressure map file !!!\n";
		exit(1);
	}
	else {
		cout << "Pressure map " << pressure_file.c_str() << "...  opened successfully\n";
		cout << "Loading PRESSURE values for calculating REFRACTIVE_INDEX values...\n";
	}


	double pressure_val = 0.0;
		double M = 0.0;  // Coefficient of modulation.

		for (array_index a = 0; a < Nx && pressure_file_stream.good(); a++)
			for (array_index b = 0; b < Nz; b++)
				for (array_index c = 0; c < Ny; c++)
				{
					pressure_file_stream >> pressure_val;
					M = 2.0 * pezio_optical_coeff * pressure_val / (density * speed_of_sound * speed_of_sound);
					(*refractive_grid)[a][b][c] = n_background * (1 + 0.5 * M);
	#ifdef DEBUG
					cout << (*refractive_grid)[a][b][c] << endl;
	#endif
				}

	pressure_file_stream.close();
}



void RefractiveMap::loadRefractiveMap(const std::string &filename)
{

	// Ensure memory has been allocated for the pressure values that
	// will be read in from file.  That is, initCommon() has already
	// been called.
	assert(refractive_grid != NULL);
	// Verify that the calculations needed for the refractive index use values
	// that have been set in the medium's call to this method.
	assert(density > 0.0 &&
			speed_of_sound > 0 &&
			pezio_optical_coeff > 0);


	std::string file_to_open = filename;
	pressure_file_stream.open(file_to_open.c_str());

	// Check for successful opening of the file.
	if (!pressure_file_stream)
	{
		cout << "!!! Error opening pressure map file " << file_to_open.c_str() << "!!!\n";
		exit(1);
	}
	else
	{
		cout << "Pressure map " << file_to_open.c_str() << " opened successfully. ";
		cout << "Loading pressure values...\n";
	}


	//#define DEBUG


	double pressure_val = 0.0;
	double M = 0.0;  // Coefficient of modulation.

	for (array_index a = 0; a < Nx && pressure_file_stream.good(); a++)
		for (array_index b = 0; b < Nz; b++)
			for (array_index c = 0; c < Ny; c++)
			{
				pressure_file_stream >> pressure_val;
				M = 2.0 * pezio_optical_coeff * pressure_val / (density * speed_of_sound * speed_of_sound);
				(*refractive_grid)[a][b][c] = n_background * (1 + 0.5 * M);
#ifdef DEBUG
				cout << (*refractive_grid)[a][b][c] << endl;
#endif
			}

	pressure_file_stream.close();
}


void RefractiveMap::loadRefractiveMap(const std::string &filename, const int timeStep)
{

	// Ensure memory has been allocated for the pressure values that
	// will be read in from file.  That is, initCommon() has already
	// been called.
	assert(refractive_grid != NULL);
	// Verify that the calculations needed for the refractive index use values
	// that have been set in the medium's call to this method.
	assert(density > 0.0 &&
			speed_of_sound > 0 &&
			pezio_optical_coeff > 0);


	// Concatonate the values passed in to form a filename to read in.
	std::string file_to_open = filename + boost::lexical_cast<std::string>(timeStep);
	pressure_file_stream.open(file_to_open.c_str());


	// Check for successful opening of the file.
	if (!pressure_file_stream)
	{
		cout << "!!! Error opening pressure map file " << file_to_open.c_str() << "!!!\n";
		exit(1);
	}
	else
	{
		cout << "Pressure map " << file_to_open.c_str() << " opened successfully. ";
		cout << "Loading pressure values...\n";
	}


	double pressure_val = 0.0;
	double M = 0.0;  // Coefficient of modulation.

	for (array_index a = 0; a < Nx && pressure_file_stream.good(); a++)
		for (array_index b = 0; b < Nz; b++)
			for (array_index c = 0; c < Ny; c++)
			{
				pressure_file_stream >> pressure_val;
				M = 2.0 * pezio_optical_coeff * pressure_val / (density * speed_of_sound * speed_of_sound);
				(*refractive_grid)[a][b][c] = n_background * (1 + 0.5 * M);
#ifdef DEBUG
				cout << (*refractive_grid)[a][b][c] << endl;
#endif
			}

	pressure_file_stream.close();
}


// Returns a Vector3d object holding values for displacements in all axes.
double RefractiveMap::getRefractiveIndex(const Vector3d &photonLocation)
{
	int _x, _y, _z;

	boost::mutex::scoped_lock lock(m_refractive_mutex);
	{
		// Indices into the displacement grids.
		_x = photonLocation.location.x/dx - (photonLocation.location.x/dx)/Nx;
		_y = photonLocation.location.y/dy - (photonLocation.location.y/dy)/Ny;
		_z = photonLocation.location.z/dz - (photonLocation.location.z/dz)/Nz;

#ifdef DEBUG
	// Sanity check.
	assert(((_x < Nx && _x >= 0) &&
			(_y < Ny && _y >= 0) &&
			(_z < Nz && _z >= 0)) ||
			assert_msg("_x=" << _x << " _y=" << _y << " _z=" << _z << "\n"
					<< photonLocation.location.x << " "
					<< photonLocation.location.y << " "
					<< photonLocation.location.z));
#endif
	}


	return getRefractiveIndexFromGrid(_x, _y, _z);
}

// Get the refractive index by converting the caartesian coords to voxel indeces.
double RefractiveMap::getRefractiveIndex(const double x, const double y, const double z)
{


	// Indices into the displacement grids.
	int _x, _y, _z;

	boost::mutex::scoped_lock lock(m_refractive_mutex);
	{
		// Indices into the displacement grids.
		_x = x/dx - (x/dx)/Nx;
		_y = y/dy - (y/dy)/Ny;
		_z = z/dz - (z/dz)/Nz;
	}


	return getRefractiveIndexFromGrid(_x, _y, _z);
}


// Returns the individual axis displacement value from their location in the grid.
double RefractiveMap::getRefractiveIndexFromGrid(const int a, const int b, const int c)
{
	return (*refractive_grid)[(array_index)a][(array_index)b][(array_index)c];
}



