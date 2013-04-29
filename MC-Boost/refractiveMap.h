/*
 * refractiveMap.h
 *
 *  Created on: Apr 23, 2012
 *      Author: StaleyJW
 */

#ifndef REFRACTIVEMAP_H_
#define REFRACTIVEMAP_H_



#include <MatrixClasses/RealMatrix.h>



#include <MC-Boost/voxel_struct.h>

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
//    
//	RefractiveMap(const std::string &filename, const int Nx, const int Nz, const int Ny, const int grid_size);
//	RefractiveMap(const int Nx, const int Nz, const int Ny, const int grid_size);

    RefractiveMap(TRealMatrix *  pressure,
                  TRealMatrix *  rhox,
                  TRealMatrix *  rhoy,
                  TRealMatrix *  rhoz,
                  TRealMatrix *  rho0,
                  TRealMatrix *  c2,
                  const double   pezio_optical_coeff,
                  const double   n_background);
    
    
    RefractiveMap(TRealMatrix *  pressure,
                  TRealMatrix *  rhox,
                  TRealMatrix *  rhoy,
                  TRealMatrix *  rhoz,
                  TRealMatrix *  rho0,
                  TRealMatrix *  c2,
                  size_t x_pml_offset,
                  size_t y_pml_offset,
                  size_t z_pml_offset,
                  const double   pezio_optical_coeff,
                  const double   n_background);
    
	
    ~RefractiveMap();
    
    
    // Common init function for constructors of the class.
	void initCommon(VoxelAttributes voxel_dims);
    
    
    /// Updates the current values in 'refractive_map' when a new pressure matrix is obtained
    /// from a time step of the k-Wave simulation.
    void    Update_refractive_map(TRealMatrix * pressure,
                                  TRealMatrix * rhox,
                                  TRealMatrix * rhoy,
                                  TRealMatrix * rhoz,
                                  TRealMatrix * rho0,
                                  TRealMatrix * c2);
    
	/// Returns the refractive index from the grid based on localized pressure
    ///
	float 	getRefractiveIndexFromGrid(const int x, const int z, const int y);

    /// Return the refractive index from the TRealMatrix.
    float  Get_refractive_index_TRealMatrix(const int x, const int y, const int z);
    
    
	// Returns the pressure from cartesian coordinates of the photon.
	double 	getRefractiveIndex(double x, double z, double y);
	// Returns the pressure from cartesian coordinates in a position vector.
	double  getRefractiveIndex(const Vector3d &location);


	//void	setDensity(const double density) {this->density = density;}
	//void	setSOS(const double SOS) {this->speed_of_sound = SOS;}
	void	setBackgroundRefractiveIndex(const double n_bg) {this->n_background = n_bg;}
	void	setPezioOpticalCoeff(const double pzo) {this->pezio_optical_coeff = pzo;}

private:
	
	

	// The bounds of the grid [meters].
	double x_bound, y_bound, z_bound;

	/// The number of voxels in the x, y, and z directions. (meters)
	/// The voxel size. (meters)
	VoxelAttributes voxel_dims;

	// The density, speed of sound and the pezio-optical coefficient of the medium simulated in k-Wave, which are
	// used to calculate the refractive indices throughout the medium.
	//double speed_of_sound;
	double pezio_optical_coeff;
	double n_background;

	// Mutex to serialize access to the displacement arrays.
	boost::mutex m_refractive_mutex;

	// Input stream
	std::ifstream pressure_file_stream;

	// Holds the name of the text file that contains the pressure values.
	std::string pressure_file;


    /// For correct indexing into the refractive grid (taking into account the pml and size of MC grid).
    size_t X_PML_OFFSET;
    size_t Y_PML_OFFSET;
    size_t Z_PML_OFFSET;
    
    /// NEW WAY:
    /// This holds the refractive map at a single time step.  That is the
    /// refractive map as it is obtained from 'KSpaceSolver' by converting pressure values
    /// to refractive index values, without any offline processing.
    TRealMatrix * refractive_map;
    
    /// The size of the sensor mask used in the kWave simulation.  That is, the number of voxels
    /// in the 3D grid used in the simulation.
    double  sensor_mask_size;
};



#endif /* REFRACTIVEMAP_H_ */
