#include "boundary.h"

Boundary::Boundary(void)
{
	depth = 0;
	width = 0;
	boundary_energy = 0;
	
}

Boundary::Boundary(double depth, double width)
{
	this->depth = depth;
	this->width = width;
	this->boundary_energy = 0;
}


