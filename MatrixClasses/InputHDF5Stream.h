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

#ifndef INPUTHDF5STREAM_H
#define	INPUTHDF5STREAM_H

#include <string>

#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/LongMatrix.h>

#include <HDF5/HDF5_File.h>

using namespace std;

/**
 * @class TInputHDF5Stream
 * @brief Input stream for loading sensor data.
 */
class TInputHDF5Stream {
public:
    /// Constructor
    TInputHDF5Stream() : HDF5_File(NULL), HDF5_Dataset_id(H5I_BADID), Position(0,0,0) {};
    
    /// Create stream
    //virtual void CreateStream(THDF5_File & HDF5_File, const char * DatasetName,
    //                          const TDimensionSizes & TotalSize, const TDimensionSizes & ChunkSize,
    //                          const int CompressionLevel);
    /// Close stream
    virtual void CloseStream();
    
    /// Add data into stream
    virtual void ReadData   (TRealMatrix& Source_matrix, TLongMatrix& Index, float * TempBuffer);
    
    /// Destructor
    virtual ~TInputHDF5Stream();
private:
    /**
     *  Copy constructor is not allowed for public
     * @param src
     */
    TInputHDF5Stream(const TInputHDF5Stream& src);
    
    /// operator = is not allowed for public
    TInputHDF5Stream & operator = (const TInputHDF5Stream& src);
    
    /// HDF5 file handle
    THDF5_File* HDF5_File;
    /// HDF5 dataset handle
    hid_t       HDF5_Dataset_id;
    
    /// Position in the dataset
    TDimensionSizes Position;
    
};// end of TInputHDF5Stream

#endif	/* INPUTREALSTREAM_H */

