//
//  logger.cpp
//  Xcode
//
//  Created by jacob on 7/18/11.
//  Copyright 2011 BMPI. All rights reserved.
//

#include "vector3D.h"
#include "photon.h"
#include "logger.h"
#include <cmath>
using std::cos;



Logger * Logger::pInstance = 0;



Logger::Logger()
{
    photon_count = 0;
}


Logger::~Logger()
{
    exit_data_stream.close();
    absorber_data_stream.close();
    rng_seed_stream.close();
    tof_stream.close();
}

Logger * Logger::getInstance(void)
{
    if (!pInstance)
    {
        pInstance = new Logger();
    }
    
    return pInstance;
}


void Logger::openExitFile(const std::string &filename)
{
    photon_count = 0;
    // Ensure file stream is not already open.
    if (exit_data_stream.is_open())
        exit_data_stream.close();
    
    exit_data_stream.open(filename.c_str());
    if (!exit_data_stream)
    {
    	cout << "!!! ERROR: Could not open '" << filename << "' for writing.  Check directory structure.\n";
    	exit(1);
    }
}


void Logger::openTOFFile(const std::string &filename)
{
	// Ensure file stream is not already open.
	if (tof_stream.is_open())
		tof_stream.close();

	tof_stream.open(filename.c_str());
	if (!tof_stream)
	{
		cout << "!!! ERROR: Could not open '" << filename << "' for writing.  Check directory structure.\n";
		exit(1);
	}
}




void Logger::createRNGSeedFile(const std::string &filename)
{
    photon_count = 0;
    rng_seed_stream.open(filename.c_str());
    if (!rng_seed_stream)
    {
        cout << "!!! ERROR: Could not open '" << filename << "' for writing !!!\n";
        exit(1);
    }
}

void Logger::openAbsorberFile(const std::string &filename)
{
    // Ensure file stream is not already open.
    if (absorber_data_stream.is_open())
        absorber_data_stream.close();
    
    absorber_data_stream.open(filename.c_str());
    if (!absorber_data_stream)
    {
    	cout << "!!! ERROR: Could not open '" << filename << "' for writing.  Check directory structure.\n";
    	exit(1);
    }
}



/// File for writing velocity and calculated displacements obtained from a location in the medium.
/// Used for debugging.
void Logger::Open_vel_disp_file(const std::string &filename)
{
    if (velocity_displacement_stream.is_open())
        velocity_displacement_stream.close();
    
    velocity_displacement_stream.open(filename.c_str());
    if (!velocity_displacement_stream)
    {
    	cout << "!!! ERROR: Could not open '" << filename << "' for writing.  Check directory structure.\n";
    	exit(1);
    }
}


void Logger::write(double val)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);

    exit_data_stream << "val = " << val << endl;
}

void Logger::writeExitData(const boost::shared_ptr<Vector3d> photonVector)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);   
    photon_count++;
    exit_data_stream << photonVector;
    exit_data_stream << "\n";
    exit_data_stream.flush();
}


void Logger::writeExitData(const boost::shared_ptr<Vector3d> photonVector,
                   const double weight,
                   bool tagged)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);
    photon_count++;    
    // Write out the location (x,y,z), exit angle (theta), weight of photon, and whether it was
    // tagged.
    exit_data_stream << tagged << "," 
                     << weight << ","
                     << photonVector->getDirZ() << ","
                     << photonVector << "\n";
    
    exit_data_stream.flush();
}



void Logger::writeExitData(const boost::shared_ptr<Vector3d> photonVector,
                           const double weight)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);
    photon_count++;
    // Write out the location (x,y,z), exit angle (theta), weight of photon, and whether it was
    // tagged.
    exit_data_stream << weight << "," 
                     << photonVector->getDirZ() << "," 
                     << photonVector << "\n";
    
}



void Logger::writeExitData(const boost::shared_ptr<Vector3d> photonVector,
                           const double exitWeight,
                           const double transmissionAngle
                           )
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);
    photon_count++;
    // Write out the location (x,y,z), transmission angle (theta), weight of photon
    exit_data_stream << exitWeight << "," 
                     << transmissionAngle << "," 
                     << photonVector << "\n";
    
    exit_data_stream.flush();
}



void Logger::writeWeightAngleLengthCoords(const double exitWeight,
                                          const double transmissionAngle,
                                          const double modulatedPathLength,
                                          const boost::shared_ptr<Vector3d> photonLocation)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);
    photon_count++;
    // Write out the location (x,y,z), transmission angle (theta), weight of photon
    exit_data_stream << exitWeight << "," 
                     << transmissionAngle << ","
                     << modulatedPathLength << ","
                     << photonLocation << "\n";
    
    exit_data_stream.flush();
    
}

void Logger::writeWeightAngleLengthCoords(Photon &p)
{
	boost::mutex::scoped_lock lock(m_mutex);
    photon_count++;
	exit_data_stream << p.weight << " "
					 << p.currLocation->getDirX() << " "
					 << p.currLocation->getDirY() << " "
					 << p.transmission_angle << " "
					 << p.displaced_optical_path_length << " "
					 //<< p.unmodulated_optical_path_length << " "
					 << p.currLocation->location.x << " "
					 << p.currLocation->location.y << " "
					 << p.currLocation->location.z << "\n";
	exit_data_stream.flush();
}
                                          
                                  
void Logger::writePhoton(Photon *p)
{
    // Grab the lock to ensure that the logger doesn't get interrupted by a thread
    // in the middle of a write, causing the output to be corrupted.
    boost::mutex::scoped_lock lock(m_mutex);
    photon_count++;
    cout << "Logger::writePhoton() stub\n";
}


void Logger::writeAbsorberData(const double absorbedWeight)
{
    boost::mutex::scoped_lock lock(m_mutex);

    absorber_data_stream << absorbedWeight << "\n";

      
    absorber_data_stream.flush();
}


void Logger::writeRNGSeeds(const unsigned int s1, const unsigned int s2,
							const unsigned int s3, const unsigned int s4)
{
    boost::mutex::scoped_lock lock(m_mutex);
    photon_count++;
    rng_seed_stream << s1 << " " <<
                       s2 << " " <<
                       s3 << " " <<
                       s4 << " " << "\n";
    rng_seed_stream.flush();
}


void Logger::writeTOFData(const double tof)
{

	boost::mutex::scoped_lock lock(m_tof_mutex);
    photon_count++;
	tof_stream << tof << " \n";
	tof_stream.flush();
}



/// For 3-D axis.
void Logger::Write_velocity_displacement(float ux, float uy, float uz,
                                         float disp_x, float disp_y, float disp_z)
{
    boost::mutex::scoped_lock lock(m_mutex);
    velocity_displacement_stream << ux << " "
                                 << uy << " "
                                 << uz << " "
                                 << disp_x << " "
                                 << disp_y << " "
                                 << disp_z << "\n";
    velocity_displacement_stream.flush();
    
}


/// Along a single axis.
void Logger::Write_velocity_displacement(float u,
                                         float disp)
{
    boost::mutex::scoped_lock lock(m_mutex);
    velocity_displacement_stream << u << " " << disp << "\n";
    velocity_displacement_stream.flush();
    
}



// Returns the current time.
std::string
Logger::getCurrTime(void)
{
    
	// Set current time variable to be used with naming data files that are saved from the simulations.
	epoch = time(NULL);
	ptr_ts = localtime(&epoch);
    
	return (boost::lexical_cast<std::string>(ptr_ts->tm_hour) + "_" +
			boost::lexical_cast<std::string>(ptr_ts->tm_min) + "_" +
			boost::lexical_cast<std::string>(ptr_ts->tm_sec));
}

 


