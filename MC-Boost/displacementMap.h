/*
 * displacementMap.h
 *
 *  Created on: 3 aug. 2011
 *      Author: StaleyJW
 */

#ifndef DISPLACEMENTMAP_H_
#define DISPLACEMENTMAP_H_



#include <MatrixClasses/RealMatrix.h>


#include <boost/thread/mutex.hpp>
#include <boost/multi_array.hpp>
#include "vectorMath.h"
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

    
    
    /// Constructor for forming a displacement map object using k-Wave C++ version.
    DisplacementMap(TRealMatrix * velocity_x,
                    TRealMatrix * velocity_y,
                    TRealMatrix * velocity_z,
                    const float US_freq,
                    const float dt);
    
    
    /// Update the pointers to the new velocities.
    void    Update_displacement_map(TRealMatrix * velocity_x,
                                    TRealMatrix * velocity_y,
                                    TRealMatrix * velocity_z);
    
    
    /// Return the displacement from the TRealMatrix based on axis.
    double   Get_displacement_X_TRealMatrix(const int x_voxel_index, const int y_voxel_index, const int z_voxel_index)
            {
                return (displacement_map_x->GetElementFrom3D(x_voxel_index, y_voxel_index, z_voxel_index));
            }
    double   Get_displacement_Y_TRealMatrix(const int x_voxel_index, const int y_voxel_index, const int z_voxel_index)
            {
                return (displacement_map_y->GetElementFrom3D(x_voxel_index, y_voxel_index, z_voxel_index));
            }
    double   Get_displacement_Z_TRealMatrix(const int x_voxel_index, const int y_voxel_index, const int z_voxel_index)
            {
                return (displacement_map_z->GetElementFrom3D(x_voxel_index, y_voxel_index, z_voxel_index));
            }
    

    
    
    
                    
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
    double  getDisplacementFromGridX(const int x_voxel_index, const int y_voxel_index, const int z_voxel_index);
    double  getDisplacementFromGridY(const int x_voxel_index, const int y_voxel_index, const int z_voxel_index);
    double  getDisplacementFromGridZ(const int x_voxel_index, const int y_voxel_index, const int z_voxel_index);


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
    
    
    /// FIXME:
    /// These are the velocities maps which are used to calculate the instantaneous displacement for a monochromatic
    /// ultrasound wave.  It is an approximation and eventually needs to be changed to something more accurate.
    TRealMatrix * velocity_map_x;
    TRealMatrix * velocity_map_y;
    TRealMatrix * velocity_map_z;
    
    
    /// Hold the displacement values calculated from velocities obtained by k-Wave.  Updated per k-Wave time step.
    TRealMatrix * displacement_map_x;
    TRealMatrix * displacement_map_y;
    TRealMatrix * displacement_map_z;
    
    /// The central frequency of the ultrasound.
    float US_freq;
    
    /// The time-step used in the k-Wave simulation.
    float dt;
    
    /// The size of the sensor mask used in the kWave simulation.  That is, the number of voxels
    /// in the 3D grid used in the simulation.
    double  sensor_mask_size;


};

#endif /* DISPLACEMENTMAP_H_ */
