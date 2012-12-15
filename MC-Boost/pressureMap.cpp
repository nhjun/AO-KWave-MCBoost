/*
 * loadPressureMap.cpp
 *
 *  Created on: 18 mei 2011
 *      Author: StaleyJW
 */

#include "pressureMap.h"
#include "vector3D.h"
#include <boost/lexical_cast.hpp>



// It's an error to create a PressureMap object without specifying attributes,
// therefore the default constructor should never be called.
PressureMap::PressureMap()
{
	cout << "Error: PressureMap() default constructor called\n";
}


PressureMap::PressureMap(const std::string &filename, const int Nx, const int Nz, const int Ny, const double grid_size)
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
	loadPressureMap();
}


PressureMap::PressureMap(const int Nx, const int Nz, const int Ny, const double grid_size)
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


void PressureMap::initCommon(void)
{
	// Make sure the grid size (voxels in each axis) has been defined.
	assert(Nx != 0 &&
			Ny != 0 &&
			Nz != 0);

	dx = (double)x_bound / (double)Nx; // (meters)
	dy = (double)y_bound / (double)Ny; // (meters)
	dz = (double)z_bound / (double)Nz; // (meters)
	pressure_grid = new three_dim_array (boost::extents[Nx][Nz][Ny]);
}


PressureMap::~PressureMap()
{
	if (pressure_grid)
	{
		delete pressure_grid;
		pressure_grid = NULL;
	}
}


void PressureMap::loadPressureMap(void)
{

	// Assure memory has been allocated for the pressure values that
	// will be read in from file.  That is, initCommon() has already
	// been called.
	assert(pressure_grid != NULL);

	// Open the file that contains the pressure values.
	pressure_file_stream.open(pressure_file.c_str());

	if (!pressure_file_stream) {
		cout << "!!! Error opening pressure map file !!!\n";
		exit(1);
	}
	else {
		cout << "Pressure map " << pressure_file.c_str() << "...  opened successfully\n";
		cout << "Loading pressure values...\n";
	}


	double data = 0.0;

	for (array_index a = 0; a < Nx && pressure_file_stream.good(); a++)
		for (array_index b = 0; b < Nz; b++)
			for (array_index c = 0; c < Ny; c++)
			{
				pressure_file_stream >> data;
				(*pressure_grid)[a][b][c] = data;
#ifdef DEBUG
				cout << (*pressure_grid)[a][b][c] << endl;
#endif
			}


	pressure_file_stream.close();
}


void PressureMap::loadPressureMap(const std::string &filename)
{
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


	double data = 0.0;


	for (array_index a = 0; a < Nx && pressure_file_stream.good(); a++)
		for (array_index b = 0; b < Nz; b++)
			for (array_index c = 0; c < Ny; c++)
			{
				pressure_file_stream >> data;
				(*pressure_grid)[a][b][c] = data;
				if (data == 643110)
					cout << "a,b,c = " << a << "," << b << "," << c << endl;
#ifdef DEBUG
				cout << (*pressure_grid)[a][b][c] << endl;
#endif
			}

	pressure_file_stream.close();
}


void PressureMap::loadPressureMap(const std::string &filename, const int timeStep)
{

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


	double data = 0.0;


	for (array_index a = 0; a < Nx && pressure_file_stream.good(); a++)
		for (array_index b = 0; b < Nz; b++)
			for (array_index c = 0; c < Ny; c++)
			{
				pressure_file_stream >> data;
				(*pressure_grid)[a][b][c] = data;
#ifdef DEBUG
				cout << (*pressure_grid)[a][b][c] << endl;
#endif
			}

	pressure_file_stream.close();
}


// FIXME: ENSURE INDICES ARE WITHIN THE DIMENSIONS OF THE GRID.
double PressureMap::getPressureFromGrid(int a, int b, int c)
{
	//	array_index x = x_location;
	//	array_index z = z_location;
	//	array_index y = y_location;
	//	return (*pressure_grid)[x][z][y];
	//	cout << "PressureMap::getPressureFromGrid\n";
	//	cout << "a=" << a << ", b=" << b << ", c =" << c  << endl;
	return (*pressure_grid)[(array_index)a][(array_index)b][(array_index)c];
}

// Returns the pressure from the grid based on supplied coordinates.
double PressureMap::getPressure(double a, double b, double c)
{

	//	int _x = floor(a/dx);
	//	int _y = floor(b/dz);
	//	int _z = floor(c/dy);
    int _x, _y, _z;
    
    boost::mutex::scoped_lock lock(m_pressure_grid_mutex);
    {
        _x = _y = _z = 0;
        
        _x = a/dx - (a/dx)/Nx;
        _y = b/dy - (b/dy)/Ny;
        _z = c/dz - (c/dz)/Nz;
    }
    
#ifdef DEBUG
	// Sanity check.
	assert((_x <= Nx && _x >= 0) &&
			(_y <= Ny && _y >= 0) &&
			(_z <= Nz && _z >= 0));
#endif

	//	cout << "PressureMap::getPressureCartCords\n";
	//cout << "a=" << _x << ", b=" << _z << ", c=" << _y << endl;
	return getPressureFromGrid(_x, _y, _z);

}


// Return pressure based on the location of the photon within the grid.
double PressureMap::getPressure(const Vector3d &photonLocation)
{
    int _x, _y, _z;
	// Indices into the displacement grids.
    boost::mutex::scoped_lock lock(m_pressure_grid_mutex);
    {
        _x = _y = _z = 0;
        
        _x = photonLocation.location.x/dx - (photonLocation.location.x/dx)/Nx;
        _y = photonLocation.location.y/dy - (photonLocation.location.y/dy)/Ny;
        _z = photonLocation.location.z/dz - (photonLocation.location.z/dz)/Nz;
    }
        
	// Sanity check.
	assert((_x <= Nx && _x >= 0) &&
			(_y <= Ny && _y >= 0) &&
			(_z <= Nz && _z >= 0));

	//	cout << "PressureMap::getPressureCartCords\n";
	//cout << "a=" << _x << ", b=" << _z << ", c=" << _y << endl;
	return getPressureFromGrid(_x, _y, _z);
}



