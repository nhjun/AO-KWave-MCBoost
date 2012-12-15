#ifndef COORDINATES_H
#define COORDINATES_H

// Structure for containing cartesian coordinates.
typedef struct {
	double x;
	double y;
	double z;
}  coords;


// Structure for containing spherical coordinates.
typedef struct {
    double r;
    double theta;
    double phi;
} sphereCoords;



// Structure for containing the direction cosines of the photon.
typedef struct {
    double x;
    double y;
    double z;
} directionCos;



#endif  // COORDINATES_H
