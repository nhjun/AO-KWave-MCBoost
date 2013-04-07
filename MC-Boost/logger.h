//
//  logger.h
//  Xcode
//
//  Created by jacob on 7/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//


// Logger singleton.
// NOTE:  Construction is NOT thread-safe, must be initialized in main before any threads are spawned.
#ifndef LOGGER_H
#define LOGGER_H

#include <ctime>
#include <fstream>
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <boost/thread/mutex.hpp>
#include <boost/lexical_cast.hpp>



// Forward decleration of objects.
class Photon;
class Vector3d;



class Logger 
{
public:
    static Logger * getInstance(void);
    
    void openExitFile(const std::string &filename);
    void createRNGSeedFile(const std::string &filename);
    void openAbsorberFile(const std::string &filename);
    void openTOFFile(const std::string &filename);
    void Open_vel_disp_file(const std::string &filename);
    
    void write(double val);
    void writeExitData(const boost::shared_ptr<Vector3d> vectorCoords);
    void writeExitData(const boost::shared_ptr<Vector3d> vectorCoords,
                       const double weight,
                       bool tagged);
    void writeExitData(const boost::shared_ptr<Vector3d> photonVector,
                       const double weight);
    void writeExitData(const boost::shared_ptr<Vector3d> photonVector,
                       const double weight,
                       const double transmissionAngle);
    
    
    void writeAbsorberData(const double absorbedWeight);
    
    // XXX: Finish me
    void writeAbsorberData(const double absorbedWeight,
                           const double theta,
                           const double phi);

    // Writes the photon's weight, transmission angle, modulated path length through the medium,
    // and its exit location on the exit detector window.
    //
    void writeWeightAngleLengthCoords(const double exitWeight,
                                      const double transmissionAngle,
                                      const double modulatedPathLength,
                                      const boost::shared_ptr<Vector3d> photonVector);
    void writeWeightAngleLengthCoords(const double exitWeight,
    								  const double dirx,
    								  const double diry,
    								  const double dirz,
    								  const double modulatedPathLength,
    								  const boost::shared_ptr<Vector3d> photonVector);
    void writeWeightAngleLengthCoords(Photon &p);
    
    
    /// Writes the velocity and the corresponding displacement (at a specific location) obtained from the AO_sim calculations.
    /// Allows the comparison between the two for identifying correct behavior.
    void    Write_velocity_displacement(float ux, float uy, float uz,
                                        float disp_x, float disp_y, float disp_z);  /// For 3-D axis.
    void    Write_velocity_displacement(float u,
                                        float disp);                             /// Along a single axis.


    // XXX:
    // - Does this introduce race conditions by pointing to a threaded object that could
    //   potentially have data changing in obscure ways?  Unsure, but each object is given
    //   it's own CPU "core" to run on, which means any object's state between switches (which
    //   there should be none (i.e. context switching) in a perfect world since threads == cores)
    //   should be coherent.
    //
    void    writePhoton(Photon *p);
    
    
    // Writes the seed that generated the random events that lead this photon to escape through
    // the aperture.
    //
    void    writeRNGSeeds(const unsigned int s1, const unsigned int s2,
    					const unsigned int s3, const unsigned int s4);
    
    
    // Returns the number of photons that were detected through the exit-aperture.
    //
    int     getNumDetectedPhotons(void) {return photon_count;}
    

    // Writes the time-of-flight value for the photon bundle when it exits the medium.
    void    writeTOFData(const double tof);
    
    
    // Starts the clock for timing.
    void    startClock() {start = clock();}
    
    // Ends the clock for timing.
    void    endClock() {end = ((double)clock() - start) / CLOCKS_PER_SEC;}
    
    // Returns the current time.
    std::string     getCurrTime(void);

private:
    Logger();                            // default constructor is private
    Logger(const Logger&){};             // copy constructor is private
    ~Logger();
    
    Logger& operator=(const Logger& rhs) {};  // assignment operator is private
    
    static Logger * pInstance;
    
    // The output streams associated with data for the photon and data for
    // the absorbers.
    ofstream exit_data_stream;              // Photon exit from medium data stream.
    ofstream absorber_data_stream;          // Absorber stream.
    ofstream rng_seed_stream;               // Random-number-generator stream.
    ofstream tof_stream;                    // Time-of-flight stream.
    ofstream velocity_displacement_stream;  // Velocity and displacement stream.

    
    // Tracks how many photons were detected through the aperture.
    int photon_count;
    
    boost::mutex m_mutex;
    boost::mutex m_tof_mutex;
    
    
    // Used to append time/date to saved data files and update execution time.
    time_t epoch;
    struct tm *ptr_ts;
    clock_t start, end;
};

#endif
