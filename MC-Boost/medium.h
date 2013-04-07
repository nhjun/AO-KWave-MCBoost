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

    // Set the density of the medium.
    //void	setDensity(const double density) {kwave.density = density;}

    // Set the Speed-of-Sound of the medium.
    //void	setSOS(const double sos) {kwave.speed_of_sound = sos;}

    // Set the pezio-optical coefficient of the medium.
    //void	setPezioOpticCoeff(const double eta) {this->pezio_optical_coeff = eta;}

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
    /// NOTE: There is no offline processing, pressure matrices are passed in as they are obtained from 'AO_sim'.
    void    Create_refractive_map(TRealMatrix * pressure,
                                  TRealMatrix * rhox,
                                  TRealMatrix * rhoy,
                                  TRealMatrix * rhoz,
                                  TRealMatrix * rho0,
                                  TRealMatrix * c2,
                                  float pezio_optical_coeff);
    
    /// Assigns the current velocity from k-Wave to the medium during run-time.
    /// NOTE: There is no offline processing, velocity matrices are passed in as they are obtained from 'AO_sim'.
    void    Create_displacement_map(TRealMatrix * ux,
                                    TRealMatrix * uy,
                                    TRealMatrix * uz,
                                    float US_freq,
                                    float dt);
    
    
    /// Return the velocity from the voxel coordinates.
    /// Voxel indices are specified from Monte-carlo simulation coordinates.
    ///float   Get_velocity_from_voxel(voxel_x, voxel_y, voxel_z);
    
    /// Return the displacement from the voxel coordinates.
    /// Voxel indices are specified from Monte-carlo simulation coordinates.
    ///float   Get_displacement_from_voxel(voxel_x, voxel_y, voxel_z);
    
    
    // Add a refractive map object that holds refractive index values generated from k-Wave pressures.
    void	addRefractiveMap(RefractiveMap *n_map);

    // Add a displacement map object that holds pressure generated from k-Wave
    void    addDisplacementMap(DisplacementMap *d_map);
    
	// Load the pressure data generated from K-Wave (at simulation time step 'dt') into the pressure map object.
	void 	loadPressure(std::string &filename, const int dt);

	// Load the pressure data generated from K-Wave if only the file name is given.
	void	loadPressure(std::string &filename);
    
	// Load the displacement data generated from K-Wave (at simulation time step 'dt') into the displacement map object.
    void    loadDisplacements(std::string &filename, const int dt);
    
    // Loads a pressure file from k-Wave generated data and calculates the displacement, essentially converting the data to displacements
    // representative of simulated pressures.
    // TODO: Implement this method so that simulations can account for varying attributes over time due to heating from ultrasound.
//    void	loadDisplacementsFromPressure(std::string &filename, const int dt,
//    									  const double density,
//    									  const double speed_of_sound,
//    									  const double pezio_optical_coeff,
//    									  const double background_refractive_index);

    void	loadDisplacementsFromPressure(std::string &filename, const int dt);

    // Load the pressure data generated from k-Wave and convert it to refractive index values.
    // NOTE: eta is the pezio-optical coefficient of the medium.
    void 	loadRefractiveMap(std::string &filename, const double density, const double sos, const double n_background, const double eta);
    void    loadRefractiveMap(std::string &filename, const double density, const double sos, const double eta, const double n_background, const int dt);

	// Return the pressure from the pressure grid based on cartesian coordinates
	// of the current location of the photon in the medium.
	double	getPressureFromCartCoords(double x, double z, double y);
    
    // Return the pressure from the pressure grid based on the location of the photon.
    double  getPressureFromPhotonLocation(const boost::shared_ptr<Vector3d> photonCoords);
    

    // Return the displacement vector coordinates from the location of the photon in the medium.
    boost::shared_ptr<Vector3d> getDisplacementFromPhotonLocation(const boost::shared_ptr<Vector3d> photonCoords);

	// Return the pressure from the pressure grid based on array index into
	// the matrix.
	double	getPressureFromGridCoords(int x, int z, int y);

	// Return the grid where absorption was accumulated.
	double * getPlanarGrid() {return Cplanar;}

	// Assign the array which will hold the planar absorbance values.
	void	setPlanarArray(double *planar);
	
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

    // The density of the medium.
    //double density;

    // The acoustic speed-of-sound of the medium.
    //double speed_of_sound;

    // The adiabatic piezo-optical coefficient of the medium.
    //float pezio_optical_coeff;
    

};

#endif	// MEDIUM_H

