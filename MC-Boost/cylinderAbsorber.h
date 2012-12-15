//
//  CylinderAbsorber.h
//  Xcode
//
//  Created by jacob on 6/20/11.
//

#ifndef CYLINDERABSORBER_H
#define CYLINDERABSORBER_H

#include <cmath>

class Absorber;

class CylinderAbsorber : public Absorber 
{
public:
    CylinderAbsorber(const double radius, const double x, const double y, const double z);
    CylinderAbsorber(const double radius, const boost::shared_ptr<Vector3d> top,
                                        const boost::shared_ptr<Vector3d>bottom); 
    ~CylinderAbsorber();
    
    virtual bool hitAbsorberBoundary(const boost::shared_ptr<Vector3d> photonVector);
    virtual bool inAbsorber(const boost::shared_ptr<Vector3d> photonVector);
    virtual bool crossedAbsorber(const boost::shared_ptr<Vector3d> photonVector);
    
    
    void    cartesianToCylindrical(void);
    
private:
    // The height and radius of the cylindrical absorber.
    double length;
    double radius;
    
    // Cartesian coordinates of the center location of one end
    // of the cyclinder.
    Vector3d topCenter;       
    
    // Cartesian coordinates of center location of the other end
    // of the cyclinder.
    Vector3d bottomCenter;
};

#endif