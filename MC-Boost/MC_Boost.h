//
//  MC_Boost.h
//  K_Wave_C
//
//  Created by jacob on 11/30/12.
//  Copyright (c) 2012 BMPI. All rights reserved.
//

#ifndef __K_Wave_C__MC_Boost__
#define __K_Wave_C__MC_Boost__

#include <MC-Boost/RNG.h>
#include <MC-Boost/coordinates.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <boost/thread/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout;
using std::endl;
using std::string;

#include "photon.h"

// Forward class declarations.
class   Medium;




class MC_Boost
{
public:
    // Default constructor is.
    MC_Boost();
    ~MC_Boost();
    
    
    // Get code name
    string GetCodeName() {return "MC-Boost v1.0"; };

	/// Get the number of CPU threads used in the simulation.
	size_t	Get_CPU_threads() {return NUM_THREADS;};

	    
    // Set number of threads.
    void    Set_num_threads(size_t num_threads);
    
    // Set the total number of photon packets to simulate.
    void    Set_num_photons(size_t num_photons_to_simulate);

    
    // Generates RNG seeds for photon paths that were detected.
    void    Generate_MC_RNG_seeds(Medium *m, coords LaserInjectionCoords);
    
    
    // Load all the seeds that produced paths that left the medium through the exit aperture.
    void    Load_exit_RNG_seeds();

    /// Set the filename to load RNG seeds from.
    void    Set_RNG_seed_file(const std::string rng_file)
    {
        rng_seed_file = rng_file;
    }

    
    /// Run the monte-carlo simulation using the produced seeds.
    void    Run_MC_sim_timestep(Medium *m, coords LaserInjectionCoords, int timestep);
    
    
    /// Set which mechanisms (if any at all) of acousto-optics to simulate.
    void    Simulate_displacement(bool flag)        {Params.DISPLACE            = flag;};
    void    Simulate_refractive_gradient(bool flag) {Params.REFRACTIVE_GRADIENT = flag;};
    void    Simulate_refractive_total(bool flag)    {Params.REFRACTIVE_TOTAL    = flag;};
    void    Simulate_modulation_depth(bool flag)    {Params.MODULATION_DEPTH    = flag;};
    
    /// Set whether or not to save seeds for this run.
    void    Save_RNG_seeds(bool flag)               {Params.SAVE_SEEDS = flag;};
    
    /// Set whether or not to use seeds that were previously generated.
    void    Use_RNG_seeds(bool flag)                {Params.USE_SEEDS = flag;};


    
private:
    // Number of threads that are running the simulation (max = boost::thread::hardware_concurrency()).
    // NOTE: NUM_THREADS == NUM_PHOTON_OBJS
    size_t  NUM_THREADS;
    size_t  NUM_PHOTON_OBJS;
    
    // Total amount of photons that will be simulated.
    size_t  MAX_NUM_PHOTONS;
    
    // Vector containing all the seeds that produced photon paths that left through the exit aperture.
    RNG_seed_vector exit_seeds;

    // File to which the seeds are written.
    std::string     rng_seed_file;
    
    /// File that holds the exit data, which depends on what is set to be collected.
    std::string     exit_data_file;
    
    // Booleans that dictate toggles on various mechanisms during the simulation.
    // DISPLACE             => Displace the scattering events from ultrasound pressure.
    // REFRACTIVE_GRADIENT  => Create curved trajectories based on pressure induced refractive index gradients.
    // REFRACTIVE_TOTAL     => Accumlate phase based on spatially varying refractive index values on straight light paths.
    // SAVE_SEEDS           => Save seeds that produced paths that made through the exit aperture (i.e. the detector).
    bool    DISPLACE;
    bool    REFRACTIVE_GRADIENT;
    bool    REFRACTIVE_TOTAL;
    bool    SAVE_SEEDS;
    
    MC_Parameters Params;
    
    
   

    
};

#endif /* defined(__K_Wave_C__MC_Boost__) */
