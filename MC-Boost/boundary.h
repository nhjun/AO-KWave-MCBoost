/*
 *  boundary.h
 */

#ifndef BOUNDARY_H
#define BOUNDARY_H

class Boundary 
{
private:
	// Total depth of our boundary for all layers combined.
	double depth;
	
	// Width of the boundary for all layers (assuming value given from -x to +x, -y to +y)
	double width;
	
	// Cumulative absorption value of photons that pass through the boundary
	double boundary_energy;

public:
	Boundary(void);
	Boundary(double depth, double width);
	

	// Getters
	double getWidth(void) {return width;}	
	double getDepth(void) {return depth;}
	
	// Set 'boundary_energy' to value when photon crosses the boundary
	void setBoundaryEnergy(double photon_weight) {boundary_energy += photon_weight;}
};

#endif //BOUNDARY_H
