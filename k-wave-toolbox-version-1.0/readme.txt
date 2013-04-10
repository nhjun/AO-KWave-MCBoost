____________________________________________________________________________

                                  k-Wave

                    A MATLAB toolbox for the time-domain 
                     simulation of acoustic wave fields
____________________________________________________________________________

PRODUCT OVERVIEW
____________________________________________________________________________

k-Wave is a freely available 3rd party toolbox for MATLAB developed for the 
time-domain simulation of acoustic wave propagation. It uses a grid-based 
iterative k-space pseudospectral solution method to coupled first-order 
acoustic equations. It has been designed to make tissue realistic 
simulations of ultrasound and photoacoustics fast and easy to use.

   Version 1.0, Released 13th November 2012
   Written by Bradley Treeby and Ben Cox

   Tested using:
      Mac OS X Lion: MATLAB R2011a
      Windows 7 64-bit: MATLAB R2010b, R2011a, R2011b, R2012a, R2012b
      Windows XP 32-bit: MATLAB R2007a (see Troubleshooting notes), R2009a
      Linux (Fedora 14, Ubuntu 10.04): R2010b, R2011a, R2011b

Please report bugs and suggestions on http://www.k-wave.org/forum
The toolbox may be downloaded from http://www.k-wave.org/download.php

NOTE: The photoacoustic reconstruction functions kspaceLineRecon and 
kspacePlaneRecon (all toolbox versions) do not currently work with R2012b.
____________________________________________________________________________

INSTALLATION INSTRUCTIONS
____________________________________________________________________________

WINDOWS: 

  1. Save and unpack the k-Wave zip file to a suitable folder, eg. 

     C:\Program Files\MATLAB\<version>\toolbox\k-Wave Toolbox

  2. Add this folder to the MATLAB path. This can be done by selecting 
     "Set Path" using the dropdown menus (File, Set Path) then clicking 
     "Add Folder", selecting the k-Wave Toolbox folder and clicking "save". 
     Alternatively, this can be done by adding the line

     addpath('<pathname>\k-Wave Toolbox');

     eg. addpath('C:\Program Files\MATLAB\R2009b\toolbox\k-Wave Toolbox');

     to the startup.m file. If no startup.m file exists, create one and save
     it in the MATLAB startup directory or any directory in the MATLAB path. 

  3. Restart Matlab, or run the command

     >> startup

     at the command line to add the path to the current Matlab session.

  4. Open the help browser by clicking on the blue question mark icon on 
     the menu bar, then select 'k-Wave Toolbox' under the contents tab. To 
     get started, try running some of the examples by first selecting them 
     and then clicking on 'run the file'.


LINUX: 

  1. Save and unpack the k-Wave zip file to a suitable folder, eg. 

     ~/Matlab/k-Wave

  2. Add this folder to the MATLAB path. This can be done by adding the 
     line

     addpath('<pathname>/k-Wave Toolbox');

     eg. addpath('~/Matlab/k-Wave');

     to the startup.m file. If no startup.m file exists, create one and save
     it in the MATLAB startup directory or any directory in the MATLAB path.

  3. Restart Matlab, or run the command

     >> startup

     at the command line to add the path to the current Matlab session.

  4. Open the help browser by clicking on the blue question mark icon on 
     the menu bar, then select 'k-Wave Toolbox' under the contents tab. To 
     get started, try running some of the examples by first selecting them 
     and then clicking on 'run the file'.
____________________________________________________________________________

RELEASE NOTES
____________________________________________________________________________

New Features and Changes: 
- 3D simulations can now be run using an optimised C++ code 
- the data recorded by the sensor mask is now set using 
  sensor.record = {'p', 'u', 'p_final', ...} 
- time-varying and averaged acoustic intensity can be returned by setting 
  sensor.record = {'I', 'I_avg', ...} 
- the inputs 'ReturnVelocity' and sensor.record_mode and the output 
  field_data have been deprecated 
- the time index at which the sensor starts recording can be set using 
  sensor.record_start_index 
- 'DataCast' can now be set to 'gpuArray-single' or 'gpuArray-double' to 
  automatically run simulations on a graphics processing unit (GPU) using 
  the Parallel Computing Toolbox (R2012a or later) 
- output variables are now created and returned using the data type 
  specified by 'DataCast' 
- output variables can be returned in double precision by setting 
  'DataRecast' to true 
- Cartesian sensor masks with linear interpolation are now supported for 
  all 'DataCast' settings and the interpolation speed during runtime has 
  been significantly improved 
- 'CartInterp' set to 'linear' is now the default setting for 1D, 2D, and 3D
  simulations 
- the calculation of kgrid.k within kWaveGrid has been modified to improve 
  memory efficiency 
- spectrum has been renamed to spect to avoid a naming conflict with the 
  MATLAB signal processing toolbox, and now accepts matrix inputs 
- makeTime now uses the maximum sound speed to compute dt, and the minimum 
  sound speed to compute t_end 
- unmaskSensorData function inputs have changed 

Bug Fixes: 
- bug fix in using 'PlotLayout' set to true when medium.sound_speed or 
  medium.density are scalar (generated plot error) 
- bug fix in list of data cast variables for 'DataCast' set to 'GPUsingle' 
- bug fix in setting the source scale parameter for time varying pressure 
  sources in heterogeneous media with p_mode set to 'dirichlet' 
- bug fix in setting the axis limits in kspaceFirstOrder2D with 'MeshPlot' 
  set to true 
- bug fix in makeGrid when defining the wavenumber variables kgrid.k, 
  kgrid.kx (etc) for grid sizes with an odd number of grid points 
- bug fix in makeSphere for grid sizes with an odd number of grid points 
- bug fix in getWin for type set to 'Gaussian' or 'Kaiser' with a user 
  input for 'Param' (param value was not used) 
- bug fix in the scan_line method of kWaveTransducer for scan lines formed 
  using a negative steering angle (generated error) 
- bug fix in grid2cart for 1D inputs (list of Cartesian points returned in 
  the wrong direction) 
- bug fix in computing axis units in kspaceSecondOrder 

New Functions: 
- beamPlot 
- focus 
- getFDMatrix 
- gradientSpect 
- gradientFD 
- kspaceFirstOrder3DC 
- makeLine 
- makeSphericalSection 
- reorderSensorData 
- scaleFig 
- scanConversion 
- writeMatrix 

New Examples: 
- Simulating B-mode Images Using A Phased Array
____________________________________________________________________________

LICENSE
____________________________________________________________________________

k-Wave (c) 2009-2012 Bradley Treeby and Ben Cox

The k-Wave toolbox is distributed by the copyright owners under the terms of
the GNU Lesser General Public License (LGPL) which is a set of additional 
permissions added to the GNU General Public License (GPL). The full text of 
both licenses is included with the toolbox in the folder 'license'.

The licence places copyleft restrictions on the k-Wave toolbox. Essentially,
anyone can use the software for any purpose (commercial or non-commercial), 
the source code for the toolbox is freely available, and anyone can 
redistribute the software (in its original form or modified) as long as the
distributed product comes with the full source code and is also licensed 
under the LGPL. You can make private modified versions of the toolbox 
without any obligation to divulge the modifications so long as the modified
software is not distributed to anyone else. The copyleft restrictions only 
apply directly to the toolbox, but not to other (non-derivative) software 
that simply links to or uses the toolbox. 

k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more 
details (http://www.gnu.org/licenses/lgpl.html). 

If you find the toolbox useful for your academic work, please consider 
citing:

B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation and 
reconstruction of photoacoustic wave-fields," J. Biomed. Opt., vol. 15, no. 
2, p. 021314, 2010.

and/or

B. E. Treeby, J. Jaros, A. P. Rendell, and B. T. Cox, "Modeling nonlinear 
ultrasound propagation in heterogeneous media with power law absorption 
using a k-space pseudospectral method," J. Acoust. Soc. Am., vol. 131, no. 6, 
pp. 4324-4336, 2012.

along with any other relevant publications. The first paper gives an overview 
of the toolbox with applications in photoacoustics, and the second describes 
the nonlinear ultrasound model and the C++ code. 
____________________________________________________________________________