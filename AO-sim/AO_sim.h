//
//  AO_sim.h
//  K_Wave_C
//
//  Created by jacob on 12/3/12.
//  Copyright (c) 2012 BMPI. All rights reserved.
//

#ifndef __K_Wave_C__AO_sim__
#define __K_Wave_C__AO_sim__

/**
 * Includes for k-Wave simulation
 */
#include <KSpaceSolver/KSpaceFirstOrder3DSolver.h>

/**
 * Includes for threaded Monte-Carlo simulation
 */
#include <MC-Boost/MC_Boost.h>
#include <MC-Boost/medium.h>
#include <MC-Boost/layer.h>
#include <MC-Boost/circularDetector.h>
#include <MC-Boost/coordinates.h>



#include <exception>
#include <iostream>
using namespace std;




class AO_Sim
{
public:
    AO_Sim();
    ~AO_Sim();
    
    
    
    /// Run the actual acousto-optic simulation.
    void        Run_acousto_optics_sim(TParameters * Parameters);
    
    
    /********************** k-Wave ****************************************************/
    
    /// Print k-Wave code and license.
    void        Print_kWave_code_and_license()
                {
                    assert (KSpaceSolver != NULL);
                    KSpaceSolver->PrintFullNameCodeAndLicense(stdout);
                }
    
    /// Print k-Wave header (i.e. version number).
    string      Print_kWave_header()
                {
                    assert (KSpaceSolver != NULL);
                    return KSpaceSolver->GetCodeName();
                }
    
    
    /// Print the k-Wave simulation parameters.
    void        Print_kWave_sim_params()
                {
                    assert (KSpaceSolver != NULL);
                    KSpaceSolver->PrintParametersOfSimulation(stdout);
                }
    
    void        kWave_allocate_memory();
    
    
    

    
    
    /********************** MC-Boost *************************************************/
    
    /// Print MC-Boost header (i.e. version number).
    string  Print_MCBoost_header() {return da_boost->GetCodeName();}
    
    
    /// Creates the grid for the monte-carlo simulation based upon the dimensions
    /// used in the k-Wave simulation.
    void    Create_MC_grid(TParameters * parameters);
    
    

    
    
    /// Add a layer to the medium, which defines the optical properties of that layer
    /// within the medium.
    void    Add_layer_MC_medium(Layer * layer)
            {
                assert(m_medium != NULL);
                m_medium->addLayer(layer);
            }
    void    Add_layer_MC_medium(Layer_Properties &props)
            {
                assert(m_medium != NULL);
                m_medium->addLayer(new Layer(props));
            }
    
    
    /// Add a detector to the medium, which serves as a way to track photons that leave
    /// the medium.  It acts as an exit aperture.
    void    Add_circular_detector_MC_medium(Detector_Properties &props)
            {
                assert(m_medium != NULL);
                m_medium->addDetector(new CircularDetector(props));
            }
    
    
    /// Set the location at which light will be injected into the medium.
    void    Set_laser_injection_coords(coords laser_coords)
            {
                m_Laser_injection_coords.x = laser_coords.x;
                m_Laser_injection_coords.y = laser_coords.y;
                m_Laser_injection_coords.z = laser_coords.z;
            }
    
    
    /// Set the number of threads that will be run for the monte-carlo simulation.
    void    Set_num_MC_threads(size_t num_threads)
            {
                da_boost->Set_num_threads(num_threads);
            }
    
    
    /// Set the number of photons to simulate.
    void    Set_num_photons(size_t num_photons)
            {
                da_boost->Set_num_photons(num_photons);
            }
    
    
    /// Run the monte-carlo simulation once to save seeds that produced photons that left
    /// the medium through the exit aperture.
    void    Generate_exit_seeds()
            {
                da_boost->Generate_RNG_seeds(m_medium, m_Laser_injection_coords);
            }
    
    
    /// Load the seeds that were saved to file.  These seeds produced paths of photons
    /// that eventually found their way out of the exit aperture.
    void    Load_generated_seeds()
            {
                /// Load the seeds for the monte-carlo simulation.
                da_boost->Load_exit_RNG_seeds();
            }
    
    /// Return the z-axis depth.
    double  Get_MC_Zaxis_depth() {return m_medium->getZbound();}
    double  Get_MC_Yaxis_depth() {return m_medium->getYbound();}
    double  Get_MC_Xaxis_depth() {return m_medium->getXbound();}
    
    
private:
    /// The medium shared between the monte-carlo and k-wave simulations.
    Medium      * m_medium;
    
    /// Instance of the monte-carlo simulation.
    MC_Boost    * da_boost;
    
    /// Holds the coordinates where light is injected into the medium.
    coords      m_Laser_injection_coords;
    
    // Instance of the k-Wave simulation.
    TKSpaceFirstOrder3DSolver * KSpaceSolver;

    
};

#endif /* defined(__K_Wave_C__AO_sim__) */
