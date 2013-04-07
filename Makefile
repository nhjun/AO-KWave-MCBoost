# /**
# * @File        Makefile
# * @Author      Jiri Jaros
# * @Affiliation CECS,ANU, Australia
# * @Email       jiri.jaros@anu.edu.au
# * @Comments    Linux makefile
# * 
# * @Tool        K-Space 2.13
# * @Created     24 September 2012, 10:57 AM
# * @LastModif   02 October 2012, 12:10 PM
#
# * @License: 
# * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org)
# * Copyright (C) 2012 Jiri Jaros and Bradley Treeby
# * 
# * This file is part of the k-Wave. k-Wave is free software: you can redistribute it 
# * and/or modify it under the terms of the GNU Lesser General Public License as 
# * published by the Free Software Foundation, either version 3 of the License, 
# * or (at your option) any later version.
# * 
# * k-Wave is distributed in the hope that it will be useful, but 
# * WITHOUT ANY WARRANTY; without even the implied warranty of 
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# * See the GNU Lesser General Public License for more details. 
# * 
# * You should have received a copy of the GNU Lesser General Public License 
# * along with k-Wave. If not, see <http://www.gnu.org/licenses/>.
# */

#################################################################################
#	The source codes can be compiled ONLY under Linux x64 			#
#	by GNU g++ 4.3 and newer OR Intel Compiler icpc 11			#
#	The code uses STATIC Linking by default					#
#										#
#	Necessary libraries:							#
#		- FFTW 3.0 and newer OR Intel MKL 11 and newer			#
#		- HDF5 version 1.8 and newer					#
#										#
#	How to compile libraries						#
#		- FFTW : download from "http://www.fftw.org/" 			#
#			run configure script with following parameters:		#
#			--enable-float --enable-sse --enable-openmp		#
#		- MKL  : Only if not using FFTW. Dowload from 			#
#			 "http://software.intel.com/en-us/intel-mkl" 	        #
#		- HDF5 : download from http://www.hdfgroup.org/HDF5/ 		#
#			run configure script with these parameters:	        #		
#			--enable-hl 						#
#									      	#
#	This makefile enables static and dynamic compilation using 		#
#		1) GNU gcc + FFTW  						#
#		2) Intel compiler + Intel MKL					#
#									      	#
#									      	#
#################################################################################


#################################################################################
#	  Set following flags based on your compiler and library paths 		#
#################################################################################

# Select compiler
COMPILER = GNU
#COMPILER = Intel

# static lining is deafult
#LINKING = STATIC
LINKING = DYNAMIC

#set up paths: FFT_DIR for FFTW or MKL for MKL
FFT_DIR=/usr/local
MKL_DIR=/opt/intel/composer_xe_2011_sp1/mkl
HDF5_DIR=/Users/betty/Desktop/Software/hdf5-1.8.10/hdf5

#################################################################################



############################## GNU g++ + FFTW ###################################
ifeq ($(COMPILER),GNU)
  CXX	   = /usr/local/bin/g++

  CPU_FLAGS = -mtune=native -msse4.2 -m64

  #Generic CPU (any intel and AMD 64b CPU)
  #CPU_FLAGS = -msse2 -m64

  #Intel Core i7 
  #CPU_FLAGS=-msse4.2 -m64

  # CFLAGS for running 
  #------------------------
  CXXFLAGS = -O3 -fopenmp $(CPU_FLAGS) -ffast-math -fassociative-math -Wall \
		     -I$(HDF5_DIR)/include -I$(FFT_DIR)/include -I .
  

  # CFLAGS for debugging
  #------------------------
  #CXXFLAGS = -O0 -fopenmp $(CPU_FLAGS) -Wall -g -I$(HDF5_DIR)/include -I$(FFT_DIR)/include -I .

  ifeq ($(LINKING),STATIC)

	LDFLAGS  = -fopenmp $(CPU_FLAGS) -L$(HDF5_DIR)/lib -L$(FFT_DIR)/lib -static

	LIBS     = $(FFT_DIR)/lib/libfftw3f.a  		\
		   $(FFT_DIR)/lib/libfftw3f_omp.a 	\
	  	   $(HDF5_DIR)/lib/libhdf5_hl.a 	\
		   $(HDF5_DIR)/lib/libhdf5.a 		\
	   	   -lz
    else	
	LDFLAGS  = -fopenmp $(CPU_FLAGS) -L$(HDF5_DIR)/lib -L$(FFT_DIR)/lib
	LIBS     = -lfftw3f  -lfftw3f_omp -lhdf5 -lhdf5_hl -lboost_thread -lboost_system
  endif
endif


########################### Intel Compiler + MKL ################################
ifeq ($(COMPILER),Intel)
	
  CXX	   = icpc

  #Generic CPU (any intel and AMD 64b CPU)
  CPU_FLAGS=-xsse2 -m64
  
  #Intel Core i7 
#  CPU_FLAGS=-xsse4.2 -m64
  
		
  CXXFLAGS = -fast -openmp  $(CPU_FLAGS) -Wall \
  	     -I$(HDF5_DIR)/include -I$(MKL_DIR)/include -I$(MKL_DIR)/include/fftw -I .

  ifeq ($(LINKING),STATIC)
	LDFLAGS  = -fast -openmp $(CPU_FLAGS) -L$(HDF5_DIR)/lib -mkl=parallel

	LIBS = 	$(HDF5_DIR)/lib/libhdf5_hl.a 			\
		$(HDF5_DIR)/lib/libhdf5.a 			\
		-lz
  else	
	LDFLAGS  = -openmp $(CPU_FLAGS) -L$(HDF5_DIR)/lib  -mkl=parallel
  	LIBS     = -lhdf5 -lhdf5_hl -lz
  endif
endif


################################# Compile #####################################

MC_SRCS=MC-Boost/$(wildcard *.cpp)
OBJS=$(MC_SRCS:.cpp=.o)

TARGET		= kspaceFirstOrder3D-OMP_MC-Boost

all:		$(TARGET)	


$(TARGET):	main.o 					\
		HDF5/HDF5_File.o			\
		KSpaceSolver/KSpaceFirstOrder3DSolver.o	\
		MatrixClasses/BaseFloatMatrix.o		\
		MatrixClasses/BaseLongMatrix.o		\
		MatrixClasses/ComplexMatrix.o		\
		MatrixClasses/FFTWComplexMatrix.o	\
		MatrixClasses/LongMatrix.o		\
		MatrixClasses/MatrixContainer.o		\
		MatrixClasses/OutputHDF5Stream.o	\
		MatrixClasses/RealMatrix.o		\
		MatrixClasses/UXYZ_SGXYZMatrix.o	\
		Parameters/CommandLineParameters.o	\
		Parameters/Parameters.o			\
        	AO-Sim/AO_sim.o				\
		MC-Boost/absorber.o                     \
		MC-Boost/boundary.o			\
		MC-Boost/circularDetector.o		\
		MC-Boost/cylinderAbsorber.o		\
		MC-Boost/detector.o			\
		MC-Boost/displacementMap.o		\
		MC-Boost/layer.o			\
		MC-Boost/logger.o			\
		MC-Boost/medium.o			\
		MC-Boost/photon.o			\
		MC-Boost/pressureMap.o			\
		MC-Boost/refractiveMap.o		\
		MC-Boost/sphereAbsorber.o		\
		MC-Boost/vector3D.o			\
		MC-Boost/MC_Boost.o			\
		MC-Boost/RNG.o				\
		

	$(CXX) $(LDFLAGS) main.o 			\
		HDF5/HDF5_File.o			\
		KSpaceSolver/KSpaceFirstOrder3DSolver.o	\
		MatrixClasses/BaseFloatMatrix.o		\
		MatrixClasses/BaseLongMatrix.o		\
		MatrixClasses/ComplexMatrix.o		\
		MatrixClasses/FFTWComplexMatrix.o	\
		MatrixClasses/LongMatrix.o          	\
		MatrixClasses/MatrixContainer.o		\
		MatrixClasses/OutputHDF5Stream.o	\
		MatrixClasses/RealMatrix.o          	\
		MatrixClasses/UXYZ_SGXYZMatrix.o	\
		Parameters/CommandLineParameters.o	\
		Parameters/Parameters.o             	\
		AO-sim/AO_sim.o				\
        	MC-Boost/absorber.o                     \
                MC-Boost/boundary.o                     \
                MC-Boost/circularDetector.o             \
                MC-Boost/cylinderAbsorber.o            	\
                MC-Boost/detector.o                     \
                MC-Boost/displacementMap.o              \
                MC-Boost/layer.o                        \
                MC-Boost/logger.o                       \
		MC-Boost/medium.o			\
                MC-Boost/photon.o                       \
                MC-Boost/pressureMap.o                  \
                MC-Boost/refractiveMap.o                \
                MC-Boost/sphereAbsorber.o               \
                MC-Boost/vector3D.o                     \
		MC-Boost/MC_Boost.o			\
		MC-Boost/RNG.o				\
		$(LIBS)                             	\
		-o $@


$(TARGET).o : $(TARGET).cpp
	$(CXX) $(CXXFLAGS) -c $(TARGET).cpp


clean:
	rm -f *.o HDF5/*.o KSpaceSolver/*.o MatrixClasses/*.o Parameters/*.o MC-Boost/*.o AO-Sim/*.o $(TARGET)


