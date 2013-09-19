/*
 * displacementMap.cpp
 *
 *  Created on: 3 aug. 2011
 *      Author: StaleyJW
 */

#include "debug.h"
#include <MC-Boost/vector3D.h>
#include "displacementMap.h"
#include <boost/lexical_cast.hpp>


// It's an error to create a DisplacementMap object without specifying attributes,
// therefore the default constructor should never be called.
DisplacementMap::DisplacementMap()
: X_PML_OFFSET(25), // defaults
  Y_PML_OFFSET(10),
  Z_PML_OFFSET(10)
{
    displacement_map_x = displacement_map_y = displacement_map_z = NULL;
}



DisplacementMap::~DisplacementMap()
{
//	if (displacement_gridX)
//		delete displacement_gridX;
//
//	if (displacement_gridY)
//		delete displacement_gridY;
//
//	if (displacement_gridZ)
//		delete displacement_gridZ;
//    
//    
//    
//    if (displacement_map_x)
//    {
//        delete displacement_map_x;
//        displacement_map_x = NULL;
//    }
//    
//    if (displacement_map_y)
//    {
//        delete displacement_map_y;
//        displacement_map_y = NULL;
//    }
//    
//    if (displacement_map_z)
//    {
//        delete displacement_map_z;
//        displacement_map_z = NULL;
//    }
//    
}



/// XXX:
/// Old approach, no longer used.  No longer using offline calculations for particle displacements.  Currently data is
/// obtained from k-Wave (c++) at runtime.
// Loads a text files containing discrete displacement values at a given time step, in all dimensions (i.e. x, y, z),
// that were obtained from kWave simulation post-processed data.
//void DisplacementMap::loadDisplacementMaps(const std::string &filename, const int timeStep)
//{
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
//}


 
// Returns the individual axis displacement value from their location in the grid.
double DisplacementMap::getDisplacementFromGridX(const int x, const int y, const int z)
{
	//return (*displacement_gridX)[(array_index)a][(array_index)b][(array_index)c];
    
    return Get_displacement_X_TRealMatrix(x + X_PML_OFFSET,
                                          y + Y_PML_OFFSET,
                                          z + Z_PML_OFFSET);
}

double DisplacementMap::getDisplacementFromGridY(const int x, const int y, const int z)
{
	//return (*displacement_gridY)[(array_index)a][(array_index)b][(array_index)c];
    
    return Get_displacement_Y_TRealMatrix(x + X_PML_OFFSET,
                                          y + Y_PML_OFFSET,
                                          z + Z_PML_OFFSET);
}

double DisplacementMap::getDisplacementFromGridZ(const int x, const int y, const int z)
{
	//return (*displacement_gridZ)[(array_index)a][(array_index)b][(array_index)c];
    
    return Get_displacement_Z_TRealMatrix(x + X_PML_OFFSET,
                                          y + Y_PML_OFFSET,
                                          z + Z_PML_OFFSET);
}



/// Invert the phase of the displacement data 180 degrees by multiplying through the matrix by -1.
void DisplacementMap::Invert_phase(void)
{
    float * raw_data_x = displacement_map_x->GetRawData();
    float * raw_data_y = displacement_map_y->GetRawData();
    float * raw_data_z = displacement_map_z->GetRawData();

    /// NOTE:
    /// - The total element count is the same for each axial component of the displacement.
    const size_t sensor_size  = displacement_map_x->GetTotalElementCount();

    for (size_t i = 0; i < sensor_size; i++)
    {
        /// Perform the inversion.
        raw_data_x[i] *= -1.0f;
        raw_data_y[i] *= -1.0f;
        raw_data_z[i] *= -1.0f;
    }
}





