#ifndef MEDIUM_H
#define MEDIUM_H

#include <vector>
#include <string>
#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
#include <fstream>

#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

#include <MC-Boost/kwave_struct.h>
#include <MC-Boost/layer.h>
#include <MC-Boost/voxel_struct.h>


//#include <MatrixClasses/MatrixContainer.h>
#include <MatrixClasses/RealMatrix.h>
//#include <MatrixClasses/ComplexMatrix.h>
#include <MatrixClasses/LongMatrix.h>
//#include <MatrixClasses/OutputHDF5Stream.h>
//#include <MatrixClasses/UXYZ_SGXYZMatrix.h>
//#include <MatrixClasses/FFTWComplexMatrix.h>




// Maximum number of bins that hold absorption values.
const int MAX_BINS = 101;


// Forward declaration of PressureMap and DisplacementMap objects.
class PressureMap;
class RefractiveMap;
class DisplacementMap;
class Detector;
class Vector3d;
class Photon;




//typedef struct {
//    // Create a pointer to a PressureMap object.  For use with
//	// modeling acousto-optics.
//	PressureMap * pmap;
//
//	// Create a pointer to a RefractiveMap object.  For use with
//	// modeling acousto-optics.
//	RefractiveMap * nmap;
//    
//    // Create a pointer to a Displacement object.  For use with
//    // modeling acousto-optics.
//    DisplacementMap * dmap;
//    
//    // Frequency of the transducer used.
//    double  transducerFreq;
//    
//    // The wavenumber of the ultrasound.
//    double  waveNumber;
//    
//    /// The density of the medium.
//    double  density;
//    
//    /// The speed-of-sound of the medium
//    double  speed_of_sound;
//    
//    /// The total number of elements in the sensor mask.
//    long    sensor_mask_index_size;
//
//    // Number of time steps in the simulation.
//    int totalTimeSteps;
//    
//    
//} kWaveSim;



// Medium is a container object that holds one or many layer objects that the
// photon is propagated through.  This allows easy simulation of heterogeneous 
// media with Monte Carlo simulations.
class Medium
{
	
public:

	friend class Photon;

	// A structure that holds K-Wave simulation data and attributes
	// This is public since the members of the struct are objects with
	// their own private methods.  This is simply a convenient container.
	kWaveSim kwave;

	Medium(void);
    Medium(const double x, const double y, const double z);
	~Medium();
    
    // Common initializations for the Medium object.  Called from constructors.
    void    initCommon(void);

    // Set acoustic properties of medium.
    // NOTE: 'eta' is the pezio-optical coefficient
    void	setDensitySOSPezio(const double density, const double SOS, const double eta);


    // Set the background refractive index.
    void	setBackgroundRefractiveIndex(const double n_background) {this->background_refractive_index = n_background;}
	
	// Add some portion of the photon's energy that was lost at this interaction
	// point (i.e. due to absorption) to the medium's grid.
	void	absorbEnergy(const double z, const double energy);
	
	// Same as above, only the argument is an array of absorbed energy values
	// that is copied entirely to the Medium.
	void	absorbEnergy(const double *energy_array);

	// Print the grid for this medium.
	void	printGrid(const int num_photons);
    
    // Set the number of time steps that occurred in the K-Wave simulation.
    void    setKWaveTimeSteps(const int timeSteps) {kwave.totalTimeSteps = timeSteps;}
	
	// Add a layer to the medium.
	void	addLayer(Layer *layer);
    void    addLayer(Layer_Properties props);
    
    // Add a detector to the medium.
    void    addDetector(Detector *detector);
    
    // See if photon has crossed the detector plane.
    int    photonHitDetectorPlane(const boost::shared_ptr<Vector3d> p0);
	
	// Add a pressure map object that holds pressure generated from K-Wave
    void 	addPressureMap(PressureMap *p_map);
    
    
    /// Assigns the current pressure from k-Wave to the medium during run-time.
    /// NOTE: There is no offline processing, pressure matrices are passed in as they are obtained from 'AO_sim' or loaded
    ///       in from a previous run of kWave.
    void    Create_refractive_map(TRealMatrix * refractive_total);
    /// Used with bending of photon trajectories.
    void    Create_refractive_map(TRealMatrix * refractive_x,
                                  TRealMatrix * refractive_y,
                                  TRealMatrix * refractive_z);
    
    /// Assigns the current velocity from k-Wave to the medium during run-time.
    /// NOTE: There is no offline processing, velocity matrices are passed in as they are obtained from 'AO_sim', or loaded
    ///       in from a previous run of kWave.
    void    Create_displacement_map(TRealMatrix * disp_x,
                                    TRealMatrix * disp_y,
                                    TRealMatrix * disp_z);

//    /// Takes the refractive map (total) computed from KSpaceSolver directly, or loaded in from a previous run of kWave.
//    void    Assign_refractive_map(TRealMatrix * refractive_total);
//    
//    /// Takes the refractive maps (gradient) computed from KSpaceSolver directly, or loaded in from a previous run of kWave.
//    void    Assign_refractive_map(TRealMatrix * refractive_x,
//                                  TRealMatrix * refractive_y,
//                                  TRealMatrix * refractive_z);
//    
//    /// Takes the displacement maps computed from KSpaceSolver directly, or loaded in from a previous run of kWave.
//    void    Assign_displacement_map(TRealMatrix * disp_x,
//                                    TRealMatrix * disp_y,
//                                    TRealMatrix * disp_z);
    
    // Add a refractive map object that holds refractive index values generated from k-Wave pressures.
    void	addRefractiveMap(RefractiveMap *n_map);

    // Add a displacement map object that holds pressure generated from k-Wave
    void    addDisplacementMap(DisplacementMap *d_map);
	
	// Returns the absorption coefficient (mu_a) for a given depth (i.e. a layer).
	double	getLayerAbsorptionCoeff(double depth);
	
	// Returns the scattering coefficient (mu_s) for a given depth (i.e. a layer).
	double	getLayerScatterCoeff(double depth);
	
	// Return the anisotropy ('g') value for a given depth (i.e. a layer).
	double	getAnisotropyFromDepth(double depth);
	
	// Return layer from depth passed in.
	Layer * getLayerFromDepth(double depth);

	// Return the layer above the current layer.
	Layer * getLayerAboveCurrent(Layer *currentLayer);

	// Return the layer below the current layer.
	Layer * getLayerBelowCurrent(double depth);

    // Return the max depth of the medium.
    double 	getDepth() {return depth;}
    
    // Return the refractive index of the medium.
    //double  getRefractiveIndex(void) {return refractive_index;}
    
    // Write photon coordinates to file.
    void 	writePhotonCoords(std::vector<double> &coords);

    // Write photon exit locations and phases to file.
    void	writeExitCoordsAndLength(std::vector<double> &coords_phase);

    // Write the photon exit locations, phase and weight to file.
    void	writeExitCoordsLengthWeight(std::vector<double> &coords_phase_weight);

    // Return the bounds of the medium.
    double getXbound(void) {return x_bound;}
    double getYbound(void) {return y_bound;}
    double getZbound(void) {return z_bound;}
    
    
    // Set Number of voxels in the medium.
    void Set_Nx(size_t Nx) {this->Nx = Nx;}
    void Set_Ny(size_t Ny) {this->Ny = Ny;}
    void Set_Nz(size_t Nz) {this->Nz = Nz;}
    
    // Set the size of voxels in the medium.
    void Set_dx(double dx) {this->dx = dx;}
    void Set_dy(double dy) {this->dy = dy;}
    void Set_dz(double dz) {this->dz = dz;}
    
    // Set the overall dimensions of the grid (e.g. 2x2x2 cm).
    void Set_dim_x(double x) {this->x_bound = x;}
    void Set_dim_y(double y) {this->y_bound = y;}
    void Set_dim_z(double z) {this->z_bound = z;}
    
    
    /// Set the size of the PML used in the k-Wave simulation.
    /// This is used to correctly offset into the displacement
    /// and refractive index grids.
    void Set_X_PML(size_t x_pml_offset) {X_PML_OFFSET = x_pml_offset;};
    void Set_Y_PML(size_t y_pml_offset) {Y_PML_OFFSET = y_pml_offset;};
    void Set_Z_PML(size_t z_pml_offset) {Z_PML_OFFSET = z_pml_offset;};
    
    
    // Get dimensions of grid.
    size_t Get_Nx() {return Nx;}
    size_t Get_Ny() {return Ny;}
    size_t Get_Nz() {return Nz;}
    
    double Get_dx() {return dx;}
    double Get_dy() {return dy;}
    double Get_dz() {return dz;}
	
private:
	
	// The arrays that hold the weights dropped during interaction points.
	//double	Cplanar[MAX_BINS];		// Planar photon concentration.
	double *Cplanar;

	
    // The total depth of the medium (meters).
    double depth;
    double x_bound,
           y_bound,
           z_bound;
	
	// Create a STL vector to hold the layers of the medium.
    std::vector<Layer *> p_layers;
    
    // Create a STL vector to hold the detectors in the medium.
    std::vector<Detector *> p_detectors;
    
	// Mutex to serialize access to the sensor array.
	boost::mutex m_sensor_mutex;
    
    /// Mutex for layers.
    boost::mutex m_layer_above_mutex;
    boost::mutex m_layer_below_mutex;
    
	// Mutex to serialize access to the data file that is written
	// by photons.
	boost::mutex m_data_file_mutex;


	// File for dumping photon paths to.  Used in the Photon class.
    std::ofstream coords_file;

	// File for dumping data regarding exit location, path length, weight, etc.
	// to file for post processing in matlab.
    std::ofstream photon_data_file;
    
    // The refrective index outside of the medium.  We assume air.
    double refractive_index;
    
    // The unmodulated (i.e. background) refractive index of the medium.
    double background_refractive_index;


    // The voxel sizes of the medium.  This will match the size of the k-Wave simulation voxel size
    // and is only set here for convenience and later use in the 'Photon' class.
    double dx, dy, dz;
    size_t Nx, Ny, Nz;
    
    /// Create a small container for voxel size and numbers.
    VoxelAttributes voxel_dims;
    
    
    /// The size (voxels) of the PML used in k-Wave simulation.
    /// Used to correctly index into displacement/refractive grids.
    size_t X_PML_OFFSET;
    size_t Y_PML_OFFSET;
    size_t Z_PML_OFFSET;

    // The density of the medium.
    //double density;

    // The acoustic speed-of-sound of the medium.
    //double speed_of_sound;

    // The adiabatic piezo-optical coefficient of the medium.
    //float pezio_optical_coeff;
    

};

#endif	// MEDIUM_H

