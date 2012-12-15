//
//  vectorMath.h
//  Xcode
//
//  Created by jacob on 7/12/11.
//  Copyright 2011 BMPI. All rights reserved.
//


#ifndef VECTORMATH_H
#define VECTORMATH_H

#include "vector3D.h"
#include "coordinates.h"
#include "boost/smart_ptr.hpp"
#include <cmath>

class Vector3d;

namespace VectorMath 
{

    inline boost::shared_ptr<Vector3d> crossProduct(const boost::shared_ptr<Vector3d> A, const boost::shared_ptr<Vector3d> B)
    {
        boost::shared_ptr<Vector3d> result(new Vector3d);
        
        result->location.x = (A->location.y*B->location.z) - (B->location.y*A->location.z);
        result->location.y = (B->location.x*A->location.z) - (A->location.x*B->location.z);
        result->location.z = (A->location.x*B->location.y) - (A->location.y*B->location.x);
        
        return result;
    };
    
    inline Vector3d crossProduct(const Vector3d &A, const Vector3d &B)
    {
        Vector3d result;
        
        result.location.x = (A.location.y*B.location.z) - (B.location.y*A.location.z);
        result.location.y = (B.location.x*B.location.z) - (B.location.x*A.location.z);
        result.location.z = (A.location.x*B.location.y) - (A.location.y*B.location.x);
        
        return result;
    };
    
    
    
    inline double dotProduct(const boost::shared_ptr<Vector3d> A, const boost::shared_ptr<Vector3d> B)
    {
		return (A->location.x*B->location.x +
                A->location.y*B->location.y +
                A->location.z*B->location.z);
    };
    
    inline double dotProduct(const Vector3d &A, const Vector3d &B)
    {
        return (A.location.x*B.location.x +
                A.location.y*B.location.y +
                A.location.z*B.location.z);
    }
    
    
    
    inline double Distance(const boost::shared_ptr<Vector3d> A, const boost::shared_ptr<Vector3d> B)
    {
        double disX = A->location.x - B->location.x;
		double disY = A->location.y - B->location.y;
		double disZ = A->location.z - B->location.z;
        
		return sqrt(disX*disX + disY*disY + disZ*disZ);
    };
    
    inline double Distance(const Vector3d &A, const Vector3d &B)
    {
        double disX = A.location.x - B.location.x;
		double disY = A.location.y - B.location.y;
		double disZ = A.location.z - B.location.z;
        
		return sqrt(disX*disX + disY*disY + disZ*disZ);
    }
    
    inline Vector3d subtractCoords(const Vector3d &A, const Vector3d &B)
    {
        Vector3d result;
        
        result.location.x = A.location.x - B.location.x;
        result.location.y = A.location.y - B.location.y;
        result.location.z = A.location.z - B.location.z;
        
        return result;
    }
    
    
    // XXX: Is this correct?
    inline coords multiplyCoords(const coords &A, const coords &B)
    {
        coords result;
        
        result.x = A.x * B.x;
        result.y = A.y * B.y;
        result.z = A.z * B.z;
        
        return result;
    }
    
    inline double Length(const coords &A)
    {
        return sqrt(A.x*A.x + A.y*A.y + A.z*A.z);
    }
    
    
    inline double Length(const boost::shared_ptr<Vector3d> A)
    {
        double x = A->location.x;
        double y = A->location.y;
        double z = A->location.z;
        
        return sqrt(x*x + y*y + z*z);
    };
    
    
    inline void Normalize(boost::shared_ptr<Vector3d> A)
    {
        double size = Length(A);
        
        A->location.x /= size;
        A->location.y /= size;
        A->location.z /= size;
    };

}




#endif // VECTORMATH_H
