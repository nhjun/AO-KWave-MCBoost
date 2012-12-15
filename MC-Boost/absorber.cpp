//
//  absorber.cpp
//  Xcode
//
//  Created by jacob on 6/20/11.
//

#include "absorber.h"
#include "vector3D.h"
#include "logger.h"



Absorber::Absorber(const double x, const double y, const double z)
{
    center = boost::shared_ptr<Vector3d> (new Vector3d);
    center->location.x = x;
    center->location.y = y;
    center->location.z = z;

    InitCommon();
}

Absorber::Absorber(const Vector3d &c)
{
    center = boost::shared_ptr<Vector3d> (new Vector3d);
    center->location.x = c.location.x;
    center->location.y = c.location.y;
    center->location.z = c.location.z;
    
    InitCommon();
}

Absorber::Absorber(const boost::shared_ptr<Vector3d> c)
{
    center = boost::shared_ptr<Vector3d> (new Vector3d);
    center->location.x = c->location.x;
    center->location.y = c->location.y;
    center->location.z = c->location.z;
    
    InitCommon();
}

void Absorber::InitCommon(void)
{
    absorbedWeight = 0.0f;
}

void Absorber::updateAbsorbedWeight(const double absorbed)
{
    boost::mutex::scoped_lock lock(m_mutex);
    this->absorbedWeight += absorbed;
}


void Absorber::writeData(void)
{
    Logger::getInstance()->writeAbsorberData(absorbedWeight);
}

Absorber::~Absorber()
{

    // STUB
}

