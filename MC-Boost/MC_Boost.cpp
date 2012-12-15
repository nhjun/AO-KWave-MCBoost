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


MC_Boost::MC_Boost(void)
{
    // Default initializations.
    NUM_THREADS = -1;
    NUM_PHOTON_OBJS = NUM_THREADS;
    MAX_NUM_PHOTONS = -1;    
    
    // Set the various mechanisms to false as default.
    DISPLACE = REFRACTIVE_GRADIENT = SAVE_SEEDS = false;
    
    // The file that the seeds are written to after calling 'Generate_RNG_seeds'.
    rng_seed_file = "seeds_for_exit.txt";
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
    
    
    // The logger is a singleton.  To bypass any problems with using singletons in a multi-threaded application
	// initialization occurs in main before any threads are spawned.
	//
	std::string exit_data_file;
	Logger::getInstance()->createRNGSeedFile(rng_seed_file);
    
    // Set the number of threads (and number of photons to simultaneously simulate).
    //Set_num_threads(boost::thread::hardware_concurrency());
    
    // Initialize the C++ RNG.
    srand(13);
    
    // Define what will happen in the simulation.
    DISPLACE = false;
    REFRACTIVE_GRADIENT = false;
    SAVE_SEEDS = true;
    
    
    
    Photon *photons[NUM_PHOTON_OBJS];
	boost::thread *threads[NUM_THREADS];
    
    
    
    for (size_t i = 0; i < NUM_THREADS; i++)
    {
      
        // Create an instance for the pseudo-random number generator for each photon object.
        RNG rng = RNG();
        // Initialize the RNG with the independent seeds.
        rng.initRNG(rand() + 128*(i+1),
                     rand() + 128*(i+1),
                     rand() + 128*(i+1),
                     rand() + 128*(i+1));
        
        
        photons[i] = new Photon();
        threads[i] = new boost::thread(&Photon::injectPhoton, photons[i], medium,
                                       MAX_NUM_PHOTONS/NUM_THREADS, rng,
                                       LaserInjectionCoords, DISPLACE, REFRACTIVE_GRADIENT, SAVE_SEEDS);
        
    }
    
    // Join all created threads once they have done their work.
	//
	for (size_t i = 0; i < NUM_THREADS; i++)
	{
		threads[i]->join();
	}
    
    // Clean up memory.
    //
    for (size_t j = 0; j < NUM_THREADS; j++)
    {
        delete photons[j];
        delete threads[j];
        
        photons[j] = NULL;
        threads[j] = NULL;
        
    }
    
    cout << "... done\n";
    cout << "Detected: " << Logger::getInstance()->getNumDetectedPhotons() << "\n";
    
}



/// Run the monte-carlo simulation using the produced seeds.
void
MC_Boost::Run_seeded_MC_sim(Medium *medium, coords LaserInjectionCoords)
{
    // Ensure we have something to simulate.
    assert (MAX_NUM_PHOTONS > 0);
    assert (NUM_THREADS > 0);
    assert (exit_seeds.size() > 0);
    // and a medium in which to run the simulation.
    assert (medium != NULL);

    
    
    
    // Set the number of threads (and number of photons to simultaneously simulate).
    //Set_num_threads(boost::thread::hardware_concurrency());
    
    
    
    // Define what will happen in the simulation.
    DISPLACE = false;
    REFRACTIVE_GRADIENT = true;
    SAVE_SEEDS = false;
    
    
    
    Photon *photons[NUM_PHOTON_OBJS];
	boost::thread *threads[NUM_THREADS];
    
    
    /// Since each thread only runs a set of seeds that produced an exit photon,
    /// we set the iteration to only once.  That way the thread only produces that
    /// exit photon without resetting and running more, which might not exit.
    const int iterations = 1;
    
    
    /// Run the simulation for each set of seeds that produced an exit photon.
    for (RNG_seed_vector::iterator it = exit_seeds.begin(); it != exit_seeds.end(); it++)
    {
        for (size_t i = 0; i < NUM_THREADS; i++)
        {
            
            // Create an instance for the pseudo-random number generator for each photon object.
            RNG rng = RNG();
            // Initialize the RNG with the independent seeds.
            rng.initRNG(it->s1,
                        it->s2,
                        it->s3,
                        it->s4);
            
            photons[i] = new Photon();
            threads[i] = new boost::thread(&Photon::injectPhoton, photons[i], medium,
                                           iterations, rng,
                                           LaserInjectionCoords, DISPLACE, REFRACTIVE_GRADIENT, SAVE_SEEDS);
            
        }
        
        // Join all created threads once they have done their work.
        //
        for (size_t i = 0; i < NUM_THREADS; i++)
        {
            threads[i]->join();
        }
        
        // Clean up memory.
        //
        for (size_t j = 0; j < NUM_THREADS; j++)
        {
            delete photons[j];
            delete threads[j];
            
            photons[j] = NULL;
            threads[j] = NULL;
            
        }
    } // end FOR_LOOP(exit_seeds.size())
    
    cout << "... done\n";
    cout << "Detected: " << Logger::getInstance()->getNumDetectedPhotons() << "\n";
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






