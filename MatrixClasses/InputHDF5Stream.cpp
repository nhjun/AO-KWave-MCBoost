/**
 * @file        InputHDF5Stream.h
 * @author      Jacob Staley              \n
 *              MIRA Institute, University of Twente, BMPI    \n
 *              j.w.staley@utwente.nl   \n
 *
 * @brief       The header file of the class loading saved data for
 *              use with the AO simulation, or later processing
 *
 * @version     kspaceFirstOrder3D 2.13
 * @date        30, June 2013 (created) \n
 *
 *
 * @section License
 * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org).\n
 * Copyright (C) 2012 Jiri Jaros and Bradley Treeby
 *
 * This file is part of k-Wave. k-Wave is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * k-Wave is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with k-Wave. If not, see <http://www.gnu.org/licenses/>.
 */



#include <MatrixClasses/InputHDF5Stream.h>


using namespace std;

//----------------------------------------------------------------------------//
//                              Constants                                     //
//----------------------------------------------------------------------------//


//------------------------------------------------------------------------------

/**
 * Open the input stream and the dataset
 *
 */
void TInputHDF5Stream::SetHDF5File(THDF5_File & HDF5_File)
{
    this->HDF5_File = &HDF5_File;
    
    Position = TDimensionSizes(0,0,0);

}// end of SetHDF5File
//------------------------------------------------------------------------------


/**
 * Close the output stream and the dataset
 *
 */
void TInputHDF5Stream::CloseStream()
{
    
    HDF5_File->CloseDataset(HDF5_Dataset_id);
    HDF5_Dataset_id  = H5I_BADID;
    
}// end of CloseStream
//------------------------------------------------------------------------------



/**
 * Add data into the stream (usually one time step sensor data)
 *
 * @param [in] SourceMatrix - Matrix from where to pick the values
 * @param [in] Index        - Index used to pick the values
 * @param [in, out] TempBuffer - Temp buffer to make the data block contiguous
 *
 */
void TInputHDF5Stream::ReadData(const char * DatasetName, float * TempBuffer)
{
    
    
//#ifndef __NO_OMP__
//#pragma omp parallel for if (Index.GetTotalElementCount() > 1e6)
//#endif
//    for (size_t i = 0; i < Index.GetTotalElementCount(); i++) {
//        TempBuffer[i] = SourceMatrix[Index[i]];
//    }
    
    
    
    
///    TDimensionSizes ScalarSizes(1,1,1);
///    long int Nx, Ny, Nz;
///    HDF5_File.ReadCompleteDataset(Nt_Name, ScalarSizes, &Nt);
    
    HDF5_Dataset_id = HDF5_File->OpenDataset(DatasetName);
    TDimensionSizes Dims = HDF5_File->GetDatasetDimensionSizes(DatasetName);
    
    /// The size of the block of data for this timestep.
    TDimensionSizes BlockSize(Dims.X,1,1);
    
    /// Dims.X is the total number of sensor elements.
    Position.X = Dims.X;
    
    /// Read the appropriate hyperslab from the HDF5 file.
    /// Note: What ends up getting passed to ReadHyperSlab is,
    ///       Position  = {1, timestep, Nx*Ny*Nz}
    ///       BlockSize = {1, 1, Nx*Ny*Nz}
    //Position.Y = 2;
    HDF5_File->ReadHyperSlab(HDF5_Dataset_id, Position, BlockSize, TempBuffer);
    
    /// Position.Y tracks the timestep to load.
    Position.Y++;
    
    HDF5_File->CloseDataset(HDF5_Dataset_id);
    
}// end of ReadData
//------------------------------------------------------------------------------



/**
 * Destructor
 *
 */
TInputHDF5Stream::~TInputHDF5Stream()
{
    
    if (HDF5_Dataset_id != -1) HDF5_File->CloseDataset(HDF5_Dataset_id);
    
}// end of Destructor
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//                              Implementation                                //
//                             protected methods                              //
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//                              Implementation                                //
//                              private methods                               //
//----------------------------------------------------------------------------//

/// -------------------------------------------------