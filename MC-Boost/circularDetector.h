//
//  circularDetector.h
//  Xcode
//
//  Created by jacob on 7/20/11.
//  Copyright 2011 BMPI. All rights reserved.
//


#ifndef CIRCULARDETECTOR_H
#define CIRCULARDETECTOR_H

#include "detector.h"

class CircularDetector : public Detector
{
public:
    CircularDetector(const Detector_Properties &props);
    CircularDetector(const double radius, const Vector3d &centerPoint);
    CircularDetector(const double radius, const boost::shared_ptr<Vector3d> centerPoint);
    ~CircularDetector();

    virtual bool photonPassedThroughDetector(const boost::shared_ptr<Vector3d> p0,
                                             const boost::shared_ptr<Vector3d> p1);
    virtual bool photonHitDetector(const boost::shared_ptr<Vector3d> p0);
    virtual void savePhotonExitCoordinates(const boost::shared_ptr<Vector3d> exitCoords);
    void savePhotonExitData(const boost::shared_ptr<Vector3d> exitCoords,
                       const double weight,
                       const bool tagged);
    virtual void savePhotonExitWeight(void);
    
    
private:
    // Radius of the detector. (meters)
    double radius;
};

#endif
