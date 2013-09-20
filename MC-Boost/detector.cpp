//
//  detector.cpp
//  Xcode
//
//  Created by jacob on 7/18/11.
//
#include <iostream>
using std::cout;
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


Detector::Detector(const Detector_Properties &props)
{
    center.location.x = props.x_coord;
    center.location.y = props.y_coord;
    center.location.z = props.z_coord;

	/// Notify where the detector is added to the medium.
    cout << "-----------------------------------------------------\n"
         << "Adding a detector to the medium /\n"
         << "--------------------------------\n"
         << " Location: [x=" << center.location.x << ", y="
							<< center.location.y << ", z="
						    << center.location.z << "] (meters)\n";

    
    // initialize the vector normal to the plane to have
    // direction since it is only used to know the direction
    // and not location that the plane faces.
    normalVector.withDirection();

	/// Set the plane on which the detector lays.
	if (props.xy_plane)
	{
		setDetectorPlaneXY();
        cout << " plane: x-y\n";
	}
	else if (props.xz_plane)
	{	
		setDetectorPlaneXZ();
        cout << " plane: x-z\n";
	}
	else if (props.yz_plane)
	{
		setDetectorPlaneYZ();
        cout << " plane: y-z\n";
	}
	else
	{
		cout << "!!!ERROR: Detector plane has not been defined.\n";
		assert((xy_plane == true) || (xz_plane == true) || (yz_plane == true));  // One plane must be set.
	}
		
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
