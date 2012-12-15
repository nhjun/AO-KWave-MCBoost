/*
 * refractiveMap.h
 *
 *  Created on: Apr 23, 2012
 *      Author: StaleyJW
 */

#ifndef REFRACTIVEMAP_H_
#define REFRACTIVEMAP_H_

#include <boost/thread/mutex.hpp>
#include <boost/multi_array.hpp>
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <string>


// Boost array that holds the displacement values after being
// loaded in from file.
typedef boost::multi_array<double, 3> three_dim_array;
typedef three_dim_array::index array_index;

// Forward declaration of Vector3d class.
class Vector3d;


class RefractiveMap
{
public:
	RefractiveMap(const std::string &filename, const int Nx, const int Nz, const int Ny, const int grid_size);
	RefractiveMap(const int Nx, const int Nz, const int Ny, const int grid_size);

	~RefractiveMap();

	void 	loadRefractiveMap(void);
	void	loadRefractiveMap(const std::string &filename, const int time_step);
	void 	loadRefractiveMap(const std::string &filename);

	// Returns the pressure from translating photon coordinates into pressure grid
	// indices.
	double 	getRefractiveIndexFromGrid(int x, int z, int y);

	// Returns the pressure from cartesian coordinates of the photon.
	double 	getRefractiveIndex(double x, double z, double y);
	// Returns the pressure from cartesian coordinates in a position vector.
	double  getRefractiveIndex(const Vector3d &location);

	int 	getNumVoxelsXaxis(void) {return Nx;}
	int		getNumVoxelsYaxis(void) {return Ny;}
	int		getNumVoxelsZaxis(void) {return Nz;}

	double  getDx(void) {return dx;}
	double 	getDy(void) {return dy;}
	double  getDz(void) {return dz;}


	void	setDensity(const double density) {this->density = density;}
	void	setSOS(const double SOS) {this->speed_of_sound = SOS;}
	void	setBackgroundRefractiveIndex(const double n_bg) {this->n_background = n_bg;}
	void	setPezioOpticalCoeff(const double pzo) {this->pezio_optical_coeff = pzo;}

private:
	// Ensure the default constructor is private.
	RefractiveMap();

	// Common init function for constructors of the class.
	void initCommon(void);

	// The bounds of the pressure grid [cm].
	double x_bound, y_bound, z_bound;

	// The number of voxels in the x, y, and z directions. (meters)
	int Nx, Nz, Ny;

	// The voxel size. (meters)
	double dx, dz, dy;

	// The density, speed of sound and the pezio-optical coefficient of the medium simulated in k-Wave, which are
	// used to calculate the refractive indices throughout the medium.
	double density;
	double speed_of_sound;
	double pezio_optical_coeff;
	double n_background;

	// Mutex to serialize access to the displacement arrays.
	boost::mutex m_refractive_mutex;

	// Input stream
	std::ifstream pressure_file_stream;

	// Holds the name of the text file that contains the pressure values.
	std::string pressure_file;

	// Holds the pressure values obtained from k-Wave in a 3-dimensional grid
	// allowing us to index into the grid based on the coordinates of the phton
	// and retrieve the localized pressure.
	three_dim_array * refractive_grid;
};



#endif /* REFRACTIVEMAP_H_ */
