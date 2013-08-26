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
    DISPLACE = REFRACTIVE_GRADIENT = REFRACTIVE_TOTAL = SAVE_SEEDS = false;

    /// Set mechanisms.
    Params.DISPLACE             = false;
    Params.REFRACTIVE_TOTAL     = false;
    Params.REFRACTIVE_GRADIENT  = false;
    Params.MODULATION_DEPTH     = false;
    Params.SAVE_SEEDS           = false;
    Params.USE_SEEDS            = false;

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
MC_Boost::Generate_MC_RNG_seeds(Medium * medium, coords LaserInjectionCoords)
{
    

    cout << "\n........... Generating Monte-Carlo Exit Seeds .........\n";
    cout.flush();


    // The logger is a singleton.  To bypass any problems with using singletons in a multi-threaded application
	// initialization occurs in main before any threads are spawned.
	//
    Logger::getInstance()->createRNGSeedFile(rng_seed_file);


    // Initialize the C++ RNG.
    srand(13);


    /// Generate seeds from only 1 thread, ensuring the seeds used in the threaded propagation are not
    /// correlated.
    /// NOTE:
    /// - This allows for multiple threads to generate the seeds if 'THREAD_CNT' is changed, but it
    ///   allows for the possibility of producing correlated "random" streams if care isn't taken
    ///   to properly account for this fact.  Some options are to implement "leap frogging" or other
    ///   techniques to thread safe the RNG.
    size_t THREAD_CNT = 1;



    Photon *photons[NUM_PHOTON_OBJS];
	boost::thread *threads[THREAD_CNT];



	/// Create a std::vector that contains objects of type 'RNGSeed', which will be the seeds fed
    /// to the RNG during execution.
    RNG_seed_vector *seeds[NUM_THREADS];
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
        threads[i] = new boost::thread(&Photon::Inject_photon, photons[i], medium,
                                       MAX_NUM_PHOTONS, seeds[i],
                                       LaserInjectionCoords, Params);

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



    cout << "... done\n"
		 << "Simulated: " << MAX_NUM_PHOTONS << '\n'
    	 << "Detected: " << Logger::getInstance()->Get_num_detected_seed_photons() << " photons\n";

}



/// Run the monte-carlo simulation using the produced seeds.
void
MC_Boost::Run_MC_sim_timestep(Medium *medium, coords LaserInjectionCoords, int time)
{
    
    /// Update the user on what is going to simulated and set/create instances of needed objects.
    if (Params.DISPLACE)
    {
        cout << "Displacement enabled\n";
    }
    if (Params.REFRACTIVE_TOTAL)
    {
        cout << "Refractive total enabled\n";
    }
    if (Params.MODULATION_DEPTH)
    {
        cout << "Modulation depth enabled\n";
    }
    if (Params.SAVE_SEEDS)
    {
        cout << "\nSaving seeds enabled\n";
        
        /// Generate the seeds by running monte-carlo simulations with one thread.  This removes any possible
        /// errors of the RNG producing duplicate numbers, thus removing any potential correlation.
        Generate_MC_RNG_seeds(medium, LaserInjectionCoords);
        
        /// After seeds are generated, we return.  If the seeds are to be used during this run, let the caller
        /// take care of that in main().
        return;
    }
    cout.flush();

    /// The logger is a singleton.  To bypass any problems with using singletons in a multi-threaded application
	/// initialization occurs in main before any threads are spawned.
    if (Params.DISPLACE || Params.REFRACTIVE_TOTAL || Params.REFRACTIVE_GRADIENT || Params.MODULATION_DEPTH)
    {
        std::string exit_data_per_timestep = exit_data_file +
                                            "_" +
                                            boost::lexical_cast<std::string>(time) +
                                            ".dat";
        Logger::getInstance()->openExitFile(exit_data_per_timestep);
    }

    
    
    /// Create a std::vector that contains objects of type 'RNGSeed', which will be the seeds fed
    /// to the RNG during execution.
    RNG_seed_vector *seeds[NUM_THREADS];
    RNGSeeds temp;
    
    /// Create an array of photon and thread pointers, which are used to launch Boost threads.
    Photon *photons[NUM_PHOTON_OBJS];
    boost::thread *threads[NUM_THREADS];
    
    
    /// A decision is to be made at this point, and that is whether or not this simulation wants to use
    /// pre-generated seeds for the RNG, or just run normally (generate seeds and then use them, which
    /// has the potential for producing correlated random streams if running more than one thread).
    /// Depending on the parameter set, choose the correct option.
    ///
    if (Params.USE_SEEDS)
    {
        cout << "\nUsing seeds enabled\n";
        /// This is the choice of loading in precomputed seeds, which can then be used in a highly
        /// threaded scenario.
        
        /// Load in the precomputed seeds.
        Load_exit_RNG_seeds();
        
        /// To overcome the problem of continously spawning and killing threads (large overhead), we simply divide
        /// all of the seeds evenly between processor cores.  As an example, we saved 10K seeds (i.e. 10K photon bundles)
        /// and this machine has 10 CPU cores, then each core gets 1K seeds to run.
        int i;
        int k;
        RNG_seed_vector::iterator iter = exit_seeds.begin();
        
        
        /// For each thread (i.e. CPU core we are using) we load in a portion of the seeds for the RNG
        /// that runs in each 'Photon' object.
        for (i = 0; i < NUM_THREADS; ++i)
        {
            seeds[i] = new RNG_seed_vector;
            
            for (k = 0; k < (exit_seeds.size()/NUM_THREADS); k++, iter++)
            {
                temp.s1 = iter->s1;
                temp.s2 = iter->s2;
                temp.s3 = iter->s3;
                temp.s4 = iter->s4;
                
                (seeds[i])->push_back(temp);
            }
            
            photons[i] = new Photon();
            threads[i] = new boost::thread(&Photon::Inject_photon, photons[i], medium,
                                           (seeds[i])->size(), seeds[i],
                                           LaserInjectionCoords, Params);
            
        } /// End FOR LOOP (NUM_THREADS)
        
        
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
            delete seeds[j];
            
            photons[j] = NULL;
            threads[j] = NULL;
            seeds[j]   = NULL;
        }
        
    }
    else
    {
        /// Do not load any sees, simply run and produce them as normal.
        /// NOTE:
        /// - Not guaranteed to produce uncorrelated streams of "random" numbers from the RNG.
        
        /// Initialize the C++ RNG.
        /// NOTE:
        /// - This is only used to randomly seed my RNG.
        srand(13);
        
        
        for (size_t i = 0; i < NUM_THREADS; ++i)
        {
            seeds[i] = new RNG_seed_vector;
            
            temp.s1 = rand() + 128*(i+1);
            temp.s2 = rand()*rand() + 128*(i+1);
            temp.s3 = rand()*rand()*rand() + 128*(i+1);
            temp.s4 = rand()*rand()*rand()*rand() + 128*(i+1);
            
            (seeds[i])->push_back(temp);
            
            ///cout << temp.s1 << " " << temp.s2 << " " << temp.s3 << " " << temp.s4 << '\n';
            photons[i] = new Photon();
            threads[i] = new boost::thread(&Photon::Inject_photon, photons[i], medium,
                                           MAX_NUM_PHOTONS/NUM_THREADS, seeds[i],
                                           LaserInjectionCoords, Params);
            
            //        threads[i] = new boost::thread(&Photon::TESTING, photons[i], medium,
            //                                        exit_seeds.size()/NUM_THREADS, seeds[i],
            //                                        LaserInjectionCoords, Params);
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
            delete seeds[j];
            
            photons[j] = NULL;
            threads[j] = NULL;
            seeds[j]   = NULL;
        }
        
    } /// End of if (Params.USE_SEEDS)


    
    
    /// Notify the user of what was detected if any of the acousto-optics mechanisms were turned on.
    cout << "... done\n";
    if (Params.DISPLACE || Params.REFRACTIVE_TOTAL || Params.REFRACTIVE_GRADIENT || Params.MODULATION_DEPTH)
    {
        cout << "Detected: " << Logger::getInstance()->Get_num_exited_photons() << " photons\n";
    }
    
    
    
    //    /// Run the simulation for each set of seeds that produced an exit photon.
    //    /// Seeds are handed off to threads until they are exhausted.  We don't want
    //    /// to generate a large amount of threads at one time though, otherwise throughput
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






