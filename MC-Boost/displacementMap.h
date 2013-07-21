/*
 * displacementMap.h
 *
 *  Created on: 3 aug. 2011
 *      Author: StaleyJW
 */

#ifndef DISPLACEMENTMAP_H_
#define DISPLACEMENTMAP_H_



#include <MatrixClasses/RealMatrix.h>

#include <MC-Boost/vectorMath.h>

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
	DisplacementMap();

    ~DisplacementMap();
    

    /// Assign the displacements along each respective acess from data provided by KSpaceSolver
    void    Assign_displacement_map(TRealMatrix * disp_x,
                                    TRealMatrix * disp_y,
                                    TRealMatrix * disp_z)
    {
        displacement_map_x = disp_x;
        displacement_map_y = disp_y;
        displacement_map_z = disp_z;
    }
    
    
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
    
    
    /// Is simulation of the refractive gradient enabled.
    bool    IsSim_displacement(void)
            {
                return ((displacement_map_x != NULL) || (displacement_map_y != NULL) || (displacement_map_z != NULL));
            }


    
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
	// Input stream.
	std::ifstream disp_file_stream;

	// Mutex to serialize access to the displacement arrays.
	boost::mutex m_displacement_mutex;

	
    int Nx, Ny, Nz;
    double dx, dy, dz;
    
    /// To keep the displacement computational grid as the same size of the k-Wave velocity grid,
    /// but taking into account the PML, we have to offset into the displacement grid to match
    /// the monte-carlo grid.  These take care of that.
    size_t X_PML_OFFSET;
    size_t Y_PML_OFFSET;
    size_t Z_PML_OFFSET;
    
    
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
