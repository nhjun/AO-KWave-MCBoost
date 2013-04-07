//
//  RNG.cpp
//  K_Wave_C
//
//  Created by jacob on 11/30/12.
//  Copyright (c) 2012 BMPI. All rights reserved.
//

#include "RNG.h"
#include <iostream>
using std::cout;


RNG::RNG()
{
    
}


RNG::~RNG()
{
    
    
}

void
RNG::initRNG(unsigned int state1, unsigned int state2,
                        unsigned int state3, unsigned int state4)
{
    z1 = state1;
    z2 = state2;
    z3 = state3;
    z4 = state4;
    
//    cout << "\nz1=" << z1 << "\n"
//         << "z2=" << z2 << "\n"
//         << "z3=" << z3 << "\n"
//         << "z4=" << z4 << "\n\n";
}



unsigned int
RNG::TausStep(unsigned int &z, int s1, int s2, int s3, unsigned long M)
{
	unsigned int b=(((z << s1) ^ z) >> s2);
	z = (((z & M) << s3) ^ b);
	return z;
}

unsigned int
RNG::LCGStep(unsigned int &z, unsigned int A, unsigned long C)
{
	return z=(A*z+C);
}


double
RNG::HybridTaus(void)
{
	// Combined period is lcm(p1,p2,p3,p4)~ 2^121
	return 2.3283064365387e-10 * (              // Periods for the RNG.
                                  TausStep(z1, 13, 19, 12, 4294967294UL) 	^  // p1=2^31-1
                                  TausStep(z2, 2, 25, 4, 4294967288UL) 	^    // p2=2^30-1
                                  TausStep(z3, 3, 11, 17, 4294967280UL) 	^   // p3=2^28-1
                                  LCGStep(z4, 1664525, 1013904223UL)        // p4=2^32
                                  );
}

double
RNG::getRandNum()
{
    return HybridTaus();
}



RNGSeeds
RNG::getState()
{
    RNGSeeds curr_state;
    
    curr_state.s1 = z1;
    curr_state.s2 = z2;
    curr_state.s3 = z3;
    curr_state.s4 = z4;
    
    return curr_state;
}


//// Overloaded operators for working with 'RNGSeeds' structures.
RNG& RNG::operator=(const RNG &rhs)
{
    
    z1 = rhs.z1;
    z2 = rhs.z2;
    z3 = rhs.z3;
    z4 = rhs.z4;
    
    return *this;
}


