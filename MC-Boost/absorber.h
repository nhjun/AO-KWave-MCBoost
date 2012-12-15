//
//  absorber.h
//  Xcode
//
//  Created by jacob on 6/20/11.
//
#ifndef ABSORBER_H
#define ABSORBER_H

#include "vectorMath.h"
using namespace VectorMath;
#include <boost/thread/mutex.hpp>


// Forward declerations of objects.
class Vector3d;
class Logger;


class Absorber 
{
public:
    Absorber(const double x, const double y, const double z);
    Absorber(const Vector3d &location);
    Absorber(const boost::shared_ptr<Vector3d> location);
    void InitCommon(void);
    virtual ~Absorber();
    

    virtual bool hitAbsorberBoundary(const boost::shared_ptr<Vector3d> photonVector) = 0;
    virtual bool inAbsorber(const boost::shared_ptr<Vector3d> photonVector) = 0;
    virtual bool crossedAbsorber(const boost::shared_ptr<Vector3d> currPoint, 
                                 const boost::shared_ptr<Vector3d> prevPoint) = 0;
    

    double getAbsorberAbsorptionCoeff(void) {return this->mu_a;}
    double getAbsorberScatteringCoeff(void) {return this->mu_s;}
    

    void setAbsorberAbsorptionCoeff(const double mu_a) {this->mu_a = mu_a;}
    void setAbsorberScatterCoeff(const double mu_s) {this->mu_s = mu_s;}
    
    // Update absorber weight.
    void updateAbsorbedWeight(const double absorbed);
    
    // Write the absorber data out to file to be used in post-processing.
    void writeData(void);
    
    
protected:
    // The optical properties of the absorber.
    double mu_a;
    double mu_s;
    double refractive_index;
    double anisotropy;
    double absorbedWeight;
    
    // The coordinates of the center point of the absorber in the medium.
    boost::shared_ptr<Vector3d> center;
    
    // Create a mutex to serialize access when updating the absorbed weight.
    // in this absorber.
    boost::mutex m_mutex;
};

#endif
