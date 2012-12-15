//
//  CylinderAbsorber.cpp
//  Xcode
//
//  Created by jacob on 6/20/11.
//
#include "absorber.h"
#include "cylinderAbsorber.h"


// Cylinder constructor takes the radius and the top and bottom coordinates
// of the cylinder.
CylinderAbsorber::CylinderAbsorber(const double radius, const boost::shared_ptr<Vector3d> top,
                                                        const boost::shared_ptr<Vector3d> bottom)
:Absorber(center)
{
    // Calculate the length of the cylinder.
    this->length = sqrt(pow(top->location.x - bottom->location.x, 2) +
                        pow(top->location.y - bottom->location.y, 2) +
                        pow(top->location.z - bottom->location.z, 2));
    
    this->radius = radius;
    this->topCenter = (*top);
    this->bottomCenter = (*bottom);
    
}


CylinderAbsorber::~CylinderAbsorber()
{
    // STUB
}


bool CylinderAbsorber::crossedAbsorber(const boost::shared_ptr<Vector3d> photonVector)
{
    // STUB
}


bool CylinderAbsorber::hitAbsorberBoundary(const boost::shared_ptr<Vector3d> photonVector)
{
    // STUB
}

// Need to draw a line from top and bottom disc of the cylinder.
// If the length of the line from the photon to that line is 
// larger than the cylinder's radius, then return false.  Else
// return true.
// Also need to check if the photon had passed through the cylinder.
bool CylinderAbsorber::inAbsorber(const boost::shared_ptr<Vector3d> photonVector)
{
    //
}


void CylinderAbsorber::cartesianToCylindrical(void)
{
    // STUB
}