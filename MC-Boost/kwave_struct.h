#ifndef __K_Wave_C__kwave_struct__
#define __K_Wave_C__kwave_struct__

#include <MC-Boost/pressureMap.h>
#include <MC-Boost/refractiveMap.h>
#include <MC-Boost/displacementMap.h>


typedef struct {
    // Create a pointer to a PressureMap object.  For use with
	// modeling acousto-optics.
	PressureMap * pmap;
    
	// Create a pointer to a RefractiveMap object.  For use with
	// modeling acousto-optics.
	RefractiveMap * nmap;
    
    // Create a pointer to a Displacement object.  For use with
    // modeling acousto-optics.
    DisplacementMap * dmap;
    
    // Frequency of the ultrasound used.
    float  US_freq;
    
    // The wavenumber of the ultrasound.
    double  waveNumber;
    
    /// The density of the medium.
    ///double  density;
    
    /// The speed-of-sound of the medium
    //double  speed_of_sound;
    
    /// The total number of elements in the sensor mask.
    long    sensor_mask_index_size;
    
    // Number of time steps in the simulation.
    int totalTimeSteps;
    
    /// Time step (i.e. dt used in k-Wave for updating calculations).
    float dt;
    
    
} kWaveSim;








#endif