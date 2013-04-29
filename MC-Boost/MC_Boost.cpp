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
#include <MC-Boost/medium.h>
#include <MC-Boost/photon.h>
#include <MC-Boost/logger.h>
#include <MC-Boost/vector3D.h>
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

    
 
    cout << "\n........... Generating Monte-Carlo Exit Seeds .........";
    cout.flush();
    
    
    // The logger is a singleton.  To bypass any problems with using singletons in a multi-threaded application
	// initialization occurs in main before any threads are spawned.
	//
    Logger::getInstance()->createRNGSeedFile(rng_seed_file);
    
    
    // Initialize the C++ RNG.
    srand(13);
    
    
    /// Generate seeds from only 1 thread, ensuring the seeds used in the threaded propagation are not
	/// correlated.
    size_t ONLY_ONE_THREAD = 1;
	size_t THREAD_CNT = ONLY_ONE_THREAD;
    
    
    Photon *photons[NUM_PHOTON_OBJS];
	boost::thread *threads[THREAD_CNT];
    
	

	/// Create a std::vector that contains objects of type 'RNGSeed', which will be the seeds fed
	/// to the RNG during execution.
	RNG_seed_vector *seeds = new RNG_seed_vector;
	RNGSeeds temp;

    for (size_t i = 0; i < THREAD_CNT; ++i)
    {

		temp.s1 = rand() + 128*(i+1);
		temp.s2 = rand() + 128*(i+1);
		temp.s3 = rand() + 128*(i+1);
		temp.s4 = rand() + 128*(i+1);
		seeds->push_back(temp);
        
		cout << temp.s1 << " " << temp.s2 << " " << temp.s3 << " " << temp.s4 << '\n';
        photons[i] = new Photon();
//        threads[i] = new boost::thread(&Photon::injectPhoton, photons[i], medium,
//                                       MAX_NUM_PHOTONS/THREAD_CNT, rng[i],
//                                       LaserInjectionCoords, DISPLACE, REFRACTIVE_GRADIENT, SAVE_SEEDS);

		threads[i] = new boost::thread(&Photon::TESTING, photons[i], medium,
                                       MAX_NUM_PHOTONS/THREAD_CNT, seeds,
                                       LaserInjectionCoords, DISPLACE, REFRACTIVE_GRADIENT, SAVE_SEEDS);

    }
    
    // Join all created threads once they have done their work.
	//
	for (size_t i = 0; i < THREAD_CNT; i++)
	{
        if ((threads[i])->joinable())
            (threads[i])->join();
	}
    
    // Clean up memory.
    //
    for (size_t j = 0; j < THREAD_CNT; j++)
    {
        delete photons[j];
        delete threads[j];
        
        photons[j] = NULL;
        threads[j] = NULL;
    }
    
    delete seeds;
    
    cout << "... done\n"
		 << "Simulated: " << MAX_NUM_PHOTONS << '\n'
    	 << "Detected: " << Logger::getInstance()->Get_num_detected_seed_photons() << " photons\n";
    
}



/// Run the monte-carlo simulation using the produced seeds.
void
MC_Boost::Run_seeded_MC_sim_timestep(Medium *medium, coords LaserInjectionCoords, int time)
{
    
    /// Ensure we have something to simulate.
    ///
    assert (MAX_NUM_PHOTONS > 0);
    assert (NUM_THREADS > 0);
    ///assert (exit_seeds.size() > 0);
    /// and a medium in which to run the simulation.
    ///
    assert (medium != NULL);
    
    
    /// The logger is a singleton.  To bypass any problems with using singletons in a multi-threaded application
	/// initialization occurs in main before any threads are spawned.
	///
    std::string exit_data_per_timestep = exit_data_file +
                                         "_" +
                                         boost::lexical_cast<std::string>(time) +
                                         ".dat";
    Logger::getInstance()->openExitFile(exit_data_per_timestep);
    
    
//    /// Run the simulation for each set of seeds that produced an exit photon.
//    /// Seeds are handed off to threads until they are exhausted.  We don't want
//    /// to generate an large amount of threads at one time though, otherwise throughput
//    /// suffers.  So, we launch a thread group based on the number of processors on the
//    /// machine, and consume the seeds until all photons have been launched.
//    size_t i = 0;
//    size_t k = 0;
//    RNG_seed_vector::iterator iter = exit_seeds.begin();
//	
//	
//    /// Only launch a group of threads at a time, based upon the number of cores on the machine.
//    NUM_THREADS = 1;
//    RNG_seed_vector *seeds[NUM_THREADS];
//    RNGSeeds temp;
//    boost::thread *threads[NUM_THREADS];
//    Photon *photons[NUM_THREADS];
//     
//    for (i = 0; i < NUM_THREADS; ++i)
//    {
////        seeds[i] = new RNG_seed_vector;
////        
////        for (k = 0; k < (exit_seeds.size()/NUM_THREADS); k++, iter++)
////        {
////            temp.s1 = iter->s1;
////            temp.s2 = iter->s2;
////            temp.s3 = iter->s3;
////            temp.s4 = iter->s4;
////            
////            (seeds[i])->push_back(temp);
////        }
////        photons[i] = new Photon();
////        threads[i] = new boost::thread(&Photon::TESTING, photons[i], medium,
////                                       exit_seeds.size()/NUM_THREADS, seeds[i],
////                                       LaserInjectionCoords, DISPLACE, REFRACTIVE_GRADIENT, SAVE_SEEDS);
//        
//        temp.s1 = rand() + 128*(i+1);
//		temp.s2 = rand() + 128*(i+1);
//		temp.s3 = rand() + 128*(i+1);
//		temp.s4 = rand() + 128*(i+1);
//		(seeds[i])->push_back(temp);
//        
//		cout << temp.s1 << " " << temp.s2 << " " << temp.s3 << " " << temp.s4 << '\n';
//        photons[i] = new Photon();
//        threads[i] = new boost::thread(&Photon::TESTING, photons[i], medium,
//                                       exit_seeds.size()/NUM_THREADS, seeds[i],
//                                       LaserInjectionCoords, DISPLACE, REFRACTIVE_GRADIENT, SAVE_SEEDS);
//        
//
//        
//    }
//    
//
//	// Join all created threads once they have done their work.
//	//
//	for (i = 0; i < NUM_THREADS; i++)
//	{
//        if ((threads[i])->joinable())
//            (threads[i])->join();
//	}
//    
//    // Clean up memory.
//    //
//    for (size_t j = 0; j < NUM_THREADS; j++)
//    {
//        delete photons[j];
//        delete threads[j];
//        //delete seeds[j];
//        
//        photons[j] = NULL;
//        threads[j] = NULL;
//        //seeds[j]   = NULL;
//    }

    // Initialize the C++ RNG.
    srand(13);
    
	size_t THREAD_CNT = NUM_THREADS;
    
    
    Photon *photons[NUM_PHOTON_OBJS];
	boost::thread *threads[THREAD_CNT];
    
	
    
	/// Create a std::vector that contains objects of type 'RNGSeed', which will be the seeds fed
	/// to the RNG during execution.
	RNG_seed_vector *seeds[THREAD_CNT];
	RNGSeeds temp;
    
    for (size_t i = 0; i < THREAD_CNT; ++i)
    {
		seeds[i] = new RNG_seed_vector;        

		temp.s1 = rand() + 128*(i+1);
		temp.s2 = rand()*rand() + 128*(i+1);
		temp.s3 = rand()*rand()*rand() + 128*(i+1);
		temp.s4 = rand()*rand()*rand()*rand() + 128*(i+1);
		
		(seeds[i])->push_back(temp);
        
		///cout << temp.s1 << " " << temp.s2 << " " << temp.s3 << " " << temp.s4 << '\n';
        photons[i] = new Photon();
        //        threads[i] = new boost::thread(&Photon::injectPhoton, photons[i], medium,
        //                                       MAX_NUM_PHOTONS/THREAD_CNT, rng[i],
        //                                       LaserInjectionCoords, DISPLACE, REFRACTIVE_GRADIENT, SAVE_SEEDS);
        
		threads[i] = new boost::thread(&Photon::TESTING, photons[i], medium,
                                       MAX_NUM_PHOTONS/THREAD_CNT, seeds[i],
                                       LaserInjectionCoords, DISPLACE, REFRACTIVE_GRADIENT, SAVE_SEEDS);
        
    }
    
    // Join all created threads once they have done their work.
	//
	for (size_t i = 0; i < THREAD_CNT; i++)
	{
        if ((threads[i])->joinable())
            (threads[i])->join();
	}
    
    // Clean up memory.
    //
    for (size_t j = 0; j < THREAD_CNT; j++)
    {
        delete photons[j];
        delete threads[j];
		delete seeds[j];
        
        photons[j] = NULL;
        threads[j] = NULL;
		seeds[j]   = NULL;
    }
    
    
    
    cout << "... done\n";
    cout << "Detected: " << Logger::getInstance()->Get_num_exited_photons() << " photons\n";
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
    
    
    cout << "\n........... Reading Monte-Carlo Exit Seeds ...........";

    
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






