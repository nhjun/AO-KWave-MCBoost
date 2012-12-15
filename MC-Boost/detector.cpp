//
//  detector.cpp
//  Xcode
//
//  Created by jacob on 7/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "detector.h"


Detector::Detector(const double x, const double y, const double z)
{
    center.location.x = x;
    center.location.y = y;
    center.location.z = z;
    
    // initialize the vector normal to the plane to have
    // direction since it is only used to know the direction
    // and not location that the plane faces.
    normalVector.withDirection();
}



Detector::Detector(const Vector3d &centerPoint)
{
    center.location.x = centerPoint.location.x;
    center.location.y = centerPoint.location.y;
    center.location.z = centerPoint.location.z;
    
    // initialize the vector normal to the plane to have
    // direction since it is only used to know the direction
    // and not location that the plane faces.
    normalVector.withDirection();
}

Detector::Detector(const boost::shared_ptr<Vector3d> centerPoint)
{
    center.location.x = centerPoint->location.x;
    center.location.y = centerPoint->location.y;
    center.location.z = centerPoint->location.z;
    
    
    // initialize the vector normal to the plane to have
    // direction since it is only used to know the direction
    // and not location that the plane faces.
    normalVector.withDirection();
}

Detector::~Detector()
{
    
}