//
//  sphere.h
//  Xcode
//
//  Created by jacob on 6/20/11.
//

#ifndef SPHEREABSORBER_H
#define SPHEREABSORBER_H

class Absorber;

class SphereAbsorber : public Absorber
{
public:
    SphereAbsorber(const double radius, const double x, const double y, const double z);
    SphereAbsorber(const double radius, const Vector3d &c);
    SphereAbsorber(const double radius, const boost::shared_ptr<Vector3d> c);
    ~SphereAbsorber();
    
    virtual bool hitAbsorberBoundary(const boost::shared_ptr<Vector3d> photonVector);
    virtual bool inAbsorber(const boost::shared_ptr<Vector3d> photonVector);
    virtual bool crossedAbsorber(const boost::shared_ptr<Vector3d> A,
                                 const boost::shared_ptr<Vector3d> B);

    
    // Check if photon is within the radius of the absorber.
    bool inSphereVolume(const boost::shared_ptr<Vector3d> photonVector);
    sphereCoords    cartesianToSpherical(const coords &c);
    
    
private:
    double radius;
    
};

#endif // SPHERE_H