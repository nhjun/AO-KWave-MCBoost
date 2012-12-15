//
//  MC_Boost.h
//  K_Wave_C
//
//  Created by jacob on 11/30/12.
//  Copyright (c) 2012 BMPI. All rights reserved.
//

#ifndef __K_Wave_C__MC_Boost__
#define __K_Wave_C__MC_Boost__

#include "RNG.h"
#include "coordinates.h"
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

// Forward class declarations.
class Photon;
class Medium;



typedef std::vector<RNGSeeds> RNG_seed_vector;



class MC_Boost
{
public:
    // Default constructor is.
    MC_Boost();
    ~MC_Boost();
    
    
    // Get code name
    string GetCodeName() {return "MC-Boost v1.0"; };
    
    // Set number of threads.
    void    Set_num_threads(size_t num_threads);
    
    // Set the total number of photon packets to simulate.
    void    Set_num_photons(size_t num_photons_to_simulate);

    
    // Generates RNG seeds for photon paths that were detected.
    void    Generate_RNG_seeds(Medium *m, coords LaserInjectionCoords);
    
    
    // Load all the seeds that produced paths that left the medium through the exit aperture.
    void    Load_exit_RNG_seeds();
    
    
    /// Run the monte-carlo simulation using the produced seeds.
    void    Run_seeded_MC_sim(Medium *m, coords LaserInjectionCoords);
    
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
    
    // Booleans that dictate toggles on various mechanisms during the simulation.
    // DISPLACE             => Displace the scattering events from ultrasound pressure.
    // REFRACTIVE_GRADIENT  => Create curved trajectories based on pressure induced refractive index gradients.
    // SAVE_SEEDS           => Save seeds that produced paths that made through the exit aperture (i.e. the detector).
    bool    DISPLACE;
    bool    REFRACTIVE_GRADIENT;
    bool    SAVE_SEEDS;
    
    
   

    
};

#endif /* defined(__K_Wave_C__MC_Boost__) */
