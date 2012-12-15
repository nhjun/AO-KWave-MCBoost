/*
 * loadPressureMap.h
 *
 *  Created on: 18 mei 2011
 *      Author: StaleyJW
 */

#ifndef PRESSUREMAP_H_
#define PRESSUREMAP_H_

#include <boost/multi_array.hpp>
#include <boost/thread/mutex.hpp>
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::fstream;
#include <string>




// Boost array that holds the pressure values after being
// loaded in from file.
typedef boost::multi_array<double, 3> three_dim_array;
typedef three_dim_array::index array_index;


// Forward declaration of Vector3d class.
class Vector3d;



class PressureMap
{
public:

	PressureMap(const std::string &filename, const int Nx, const int Nz, const int Ny, const double grid_size);
	PressureMap(const int Nx, const int Nz, const int Ny, const double grid_size);
	~PressureMap();

	void 	loadPressureMap(void);
	void	loadPressureMap(const std::string &filename, const int time_step);
	void 	loadPressureMap(const std::string &filename);

    // Returns the pressure from translating photon coordinates into pressure grid
    // indices.
	double 	getPressureFromGrid(int x, int z, int y);
    
    // Returns the pressure from cartesian coordinates of the photon.
	double 	getPressure(double x, double z, double y);
    // Returns the pressure from cartesian coordinates in a position vector.
    double  getPressure(const Vector3d &location);
    
    int 	getNumVoxelsXaxis(void) {return Nx;}
    int		getNumVoxelsYaxis(void) {return Ny;}
    int		getNumVoxelsZaxis(void) {return Nz;}
    
    double  getDx(void) {return dx;}
    double 	getDy(void) {return dy;}
    double  getDz(void) {return dz;}

    
    



private:
	// Ensure the default constructor is private.
	PressureMap();

	// Common init function for constructors of the class.
	void initCommon(void);

	// The bounds of the pressure grid [cm].
	double x_bound, y_bound, z_bound;

	// The number of voxels in the x, y, and z directions. (meters)
	int Nx, Nz, Ny;

	// The voxel size. (meters)
	double dx, dz, dy;


	// Input stream
	std::ifstream pressure_file_stream;

	// Holds the name of the text file that contains the pressure values.
	std::string pressure_file;

	// Holds the pressure values obtained from k-Wave in a 3-dimensional grid
	// allowing us to index into the grid based on the coordinates of the phton
	// and retrieve the localized pressure.
	three_dim_array * pressure_grid;
    
    // Mutex to serialize access to the pressure grid
	boost::mutex m_pressure_grid_mutex;
    

};


#endif /* PRESSUREMAP_H_ */
