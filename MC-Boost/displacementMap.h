/*
 * displacementMap.h
 *
 *  Created on: 3 aug. 2011
 *      Author: StaleyJW
 */

#ifndef DISPLACEMENTMAP_H_
#define DISPLACEMENTMAP_H_


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


class DisplacementMap
{
public:
	DisplacementMap(const std::string &filename, const int Nx, const int Nz, const int Ny, const int grid_size);
    DisplacementMap(const int Nx, const int Nz, const int Ny, const int grid_size);

	~DisplacementMap();

	// Loads a text file containing discrete displacement values at a given time step
	// that were obtained from kWave simulation post-processed data.
	void	loadDisplacementMaps(const std::string &filename, const int timeStep);
	void	loadPressureAndCalculateDisplacements(const std::string &filename, const int dt,
			  	  	  	  	  	  	  	  	  	  const double density,
			  	  	  	  	  	  	  	  	  	  const double speed_of_sound,
			  	  	  	  	  	  	  	  	  	  const double pezio_optical_coeff,
			  	  	  	  	  	  	  	  	  	  const double background_refractive_index);


	// Returns a Vector3d object holding values for displacements in all axes.
	// That is the returned Vector3d objects holds the values the coordinates of
	// the photon should be displaced accordingly.
    boost::shared_ptr<Vector3d>	getDisplacements(const Vector3d &photonLocation);
    boost::shared_ptr<Vector3d> getDisplacements(const double x, const double y, const double z);
    
    // Returns the individual axis displacements.
    double  getDisplacementFromGridX(const int a, const int b, const int c);
    double  getDisplacementFromGridY(const int a, const int b, const int c);
    double  getDisplacementFromGridZ(const int a, const int b, const int c);


    int 	getNumVoxelsXaxis(void) {return Nx;}
    int		getNumVoxelsYaxis(void) {return Ny;}
    int		getNumVoxelsZaxis(void) {return Nz;}

    double  getDx(void) {return dx;}
    double 	getDy(void) {return dy;}
    double  getDz(void) {return dz;}

private:
	// Ensure the default constructor can never be called.
	DisplacementMap();

	// Common initialization function.
	void	initCommon();

	// Input stream.
	std::ifstream disp_file_stream;

	// The bounds of the pressure grid. (meters)
	int x_bound, y_bound, z_bound;

	// The number of voxels in the x, y, and z directions. (meters)
	int Nx, Nz, Ny;

	// The voxel size. (meters)
	double dx, dz, dy;

	// Mutex to serialize access to the displacement arrays.
	boost::mutex m_displacement_mutex;

	// Holds the displacement values obtained from k-Wave in a 3-dimensional grid
	// allowing indexing into the grid based on the coordinates of the photon
	// and retrieve the localized displacement.
	// NOTE:
	// - Displacement happens on each dimension separate from the other, based on the
	//   speed of sound in each direction.  Therefore post-processing of velocity data
	//   leaves displacement values in each direction, therefore we need 3 grids.
	three_dim_array * displacement_gridX;
	three_dim_array * displacement_gridY;
	three_dim_array * displacement_gridZ;



};

#endif /* DISPLACEMENTMAP_H_ */
