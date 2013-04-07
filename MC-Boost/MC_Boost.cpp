//
//  MC_Boost.cpp
//  K_Wave_C
//
//  Created by jacob on 11/30/12.
//  Copyright (c) 2012 BMPI. All rights reserved.
//


// This serves to encapsulate the threaded Monte-Carlo simulation
// into a single object containing all relevant parts to perform the simulation.


#include "MC_Boost.h"
#include "Medium.h"
#include "Photon.h"
#include "logger.h"
#include "vector3d.h"
#include <boost/lexical_cast.hpp>



MC_Boost::MC_Boost(void)
{
    // Default initializations.
    NUM_THREADS = -1;
    NUM_PHOTON_OBJS = NUM_THREADS;
    MAX_NUM_PHOTONS = -1;    
    
    // Set the various mechanisms to false as default.
    DISPLACE = REFRACTIVE_GRADIENT = SAVE_SEEDS = false;
    
    // The file that the seeds are written to after calling 'Generate_RNG_seeds'.
    rng_seed_file = "./Data/seeds_for_exit.dat";
    
    /// Holds data collected when a photon exits through the aperture.
    exit_data_file = "./Data/photon-exit-data";
}


MC_Boost::~MC_Boost(void)
{
    
}


void
MC_Boost::Set_num_photons(size_t num_photons_to_simulate)
{
    assert(num_photons_to_simulate > 0);
    
    MAX_NUM_PHOTONS = num_photons_to_simulate;
    
}

// Sets the number of threads to run for the monte-carlo simulation.
// This is also the number of 'Photon' objects that will be created, which each thread is responsible for.
void
MC_Boost::Set_num_threads(size_t num_threads_to_run)
{
    // Assert that we are not trying to use more cores than possible.
    assert(num_threads_to_run <= boost::thread::hardware_concurrency());
    
    NUM_PHOTON_OBJS = NUM_THREADS = num_threads_to_run;
    
}


// Runs monte-carlo simulation to save seeds that generated paths of photons that
// were detected through the detection aperture.
void
MC_Boost::Generate_RNG_seeds(Medium * medium, coords LaserInjectionCoords)
{
    // Ensure we have something to simulate.
    assert (MAX_NUM_PHOTONS > 0);
    assert (NUM_THREADS > 0);

    // and a medium in which to run the simulation.
    assert (medium != NULL);
 
    cout << "........... Generating Monte-Carlo Exit Seeds .........";
    cout.flush();
    
    
    // The logger is a singleton.  To bypass any problems with using singletons in a multi-threaded application
	// initialization occurs in main before any threads are spawned.
	//
    Logger::getInstance()->createRNGSeedFile(rng_seed_file);
    
    
    // Initialize the C++ RNG.
    srand(13);
    
    
    Photon *photons[NUM_PHOTON_OBJS];
    RNG *rng[NUM_PHOTON_OBJS];
	boost::thread *threads[NUM_THREADS];
    
    
    for (size_t i = 0; i < NUM_THREADS; i++)
    {
      
        // Create an instance for the pseudo-random number generator for each photon object.
        //RNG rng = RNG();
        //RNG *rng = new RNG();
        rng[i] = new RNG();
        
        // Initialize the RNG with the independent seeds.
        (rng[i])->initRNG(rand() + 128*(i+1),
                          rand() + 128*(i+1),
                          rand() + 128*(i+1),
                          rand() + 128*(i+1));
        //(rng[i])->initRNG(817501462, 3451853807, 1786904099, 2960543463);
        //(rng[i])->initRNG(384053651, 3206279130, 2743324459, 1296909233);
        
        photons[i] = new Photon();
        threads[i] = new boost::thread(&Photon::injectPhoton, photons[i], medium,
                                       MAX_NUM_PHOTONS/NUM_THREADS, rng[i],
                                       LaserInjectionCoords, DISPLACE, REFRACTIVE_GRADIENT, SAVE_SEEDS);
//                threads[i] = new boost::thread(&Photon::injectPhoton, photons[i], medium,
//                                               1, rng,
//                                               LaserInjectionCoords, DISPLACE, REFRACTIVE_GRADIENT, SAVE_SEEDS);
//
        
    }
    
    // Join all created threads once they have done their work.
	//
	for (size_t i = 0; i < NUM_THREADS; i++)
	{
        if ((threads[i])->joinable())
            (threads[i])->join();
	}
    
    // Clean up memory.
    //
    for (size_t j = 0; j < NUM_THREADS; j++)
    {
        delete photons[j];
        delete threads[j];
        delete rng[j];
        
        photons[j] = NULL;
        threads[j] = NULL;
        rng[j]     = NULL;
        
    }
    
    cout << "... done\n";
    cout << "Detected: " << Logger::getInstance()->getNumDetectedPhotons() << " photons\n";
    
}



/// Run the monte-carlo simulation using the produced seeds.
void
MC_Boost::Run_seeded_MC_sim_timestep(Medium *medium, coords LaserInjectionCoords, int time)
{
    
    // Ensure we have something to simulate.
    assert (MAX_NUM_PHOTONS > 0);
    assert (NUM_THREADS > 0);
    assert (exit_seeds.size() > 0);
    // and a medium in which to run the simulation.
    assert (medium != NULL);
    
    
    /// Create the photon and thread objects.
    Photon *photons[NUM_PHOTON_OBJS];
    RNG *rng[NUM_PHOTON_OBJS];
	boost::thread *threads[NUM_THREADS];
    
    
    /// Since each thread only runs a set of seeds that produced an exit photon,
    /// we set the iteration to only once.  That way the thread only produces that
    /// exit photon without resetting and running more, which might not exit through
    /// the aperture.  We can only guarantee the set of seeds will reproduce the already
    /// detected photon, not others.
    const int ITERATIONS = 1;
    
    
    /// Define what will happen in the AO simulation based on whether objects have been created that define
    /// the possible mechanisms.
    
    /// The logger is a singleton.  To bypass any problems with using singletons in a multi-threaded application
	/// initialization occurs in main before any threads are spawned.
	///
    std::string exit_data_per_timestep = exit_data_file +
                                         "_" +
                                         boost::lexical_cast<std::string>(time) +
                                         ".dat";
    Logger::getInstance()->openExitFile(exit_data_per_timestep);
    
    
    /// Run the simulation for each set of seeds that produced an exit photon.
    /// Seeds are handed off to threads until they are exhausted.  We don't want
    /// to generate an large amount of threads at one time though, otherwise throughput
    /// suffers.  So, we launch a thread group based on the number of processors on the
    /// machine, and consume the seeds until all photons have been launched.  
    size_t t_spawned = 0;
    size_t i = 0;
    RNG_seed_vector::iterator iter = exit_seeds.begin();
    boost::thread_group tgroup;
    
    
    /// Continue until all seeds have been simulated by threads.
    while (t_spawned < exit_seeds.size())
    {
        /// Only launch a group of threads at a time, based upon the number of cores on the machine.
        for (i = 0; i < NUM_THREADS; )
        {
            
            //cout << iter->s1 << " " << iter->s2 << " " << iter->s3 << " " << iter->s4 << "\n";
            // Create an instance for the pseudo-random number generator for each photon object.
            rng[i] = new RNG();
            // Initialize the RNG with the independent seeds.
            (rng[i])->initRNG(iter->s1,
                              iter->s2,
                              iter->s3,
                              iter->s4);
            
            photons[i] = new Photon();
            tgroup.create_thread(boost::bind(&Photon::injectPhoton, photons[i], medium,
                                             ITERATIONS, rng[i],
                                             LaserInjectionCoords, DISPLACE, REFRACTIVE_GRADIENT, SAVE_SEEDS));

            
            ++i;
            ++iter;
            ++t_spawned;
            if (t_spawned >= exit_seeds.size())
                break;
        }
        
        tgroup.join_all();
        
        
        
        // Clean up memory.
        //
        for (size_t k = 0; k < i; k++)
        {
            delete photons[k];
            delete rng[k];
            
            photons[k] = NULL;
            rng[k]     = NULL;
            
        }
    } // end WHILE_LOOP(exit_seeds.size())
    
    cout << "Detected: " << Logger::getInstance()->getNumDetectedPhotons() << " photons\n";
}



void
MC_Boost::Load_exit_RNG_seeds()
{
    // Open the file that has the seeds that produced photons that exited through the aperture.
	std::ifstream rng_seed_stream;
	rng_seed_stream.open(rng_seed_file.c_str());
	
    
    if (!rng_seed_stream)
	{
		cout << "!!! ERROR: Could not open (" << rng_seed_file << ") !!!\n";
		exit(1);
	}
    
    
    cout << "........... Reading Monte-Carlo Exit Seeds ...........";

    
    // Read in data from the stream.
    while (rng_seed_stream.good())
	{
		RNGSeeds temp;
		rng_seed_stream >> temp.s1;
		rng_seed_stream >> temp.s2;
		rng_seed_stream >> temp.s3;
		rng_seed_stream >> temp.s4;
        exit_seeds.push_back(temp);
	}
    
    
	// Remove the last row of seeds as they are all zeros due to
	// the way the while loop above terminates.
	exit_seeds.erase(exit_seeds.end()-1);
    
    
    rng_seed_stream.close();
    
    cout << ".... done\n";
    cout << "Number seeds read: " << exit_seeds.size() << "\n";

    
}






