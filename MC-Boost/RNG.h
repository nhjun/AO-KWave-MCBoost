//
//  RNG.h
//  K_Wave_C
//
//  Created by jacob on 11/30/12.
//  Copyright (c) 2012 BMPI. All rights reserved.
//

#ifndef __K_Wave_C__RNG__
#define __K_Wave_C__RNG__

#include <vector>



typedef struct {
	unsigned int s1;
	unsigned int s2;
	unsigned int s3;
	unsigned int s4;
} RNGSeeds;



typedef std::vector<RNGSeeds> RNG_seed_vector;



class RNG
{
public:
    RNG();
    ~RNG();
    
    void            initRNG(unsigned int state1, unsigned int state2,
                           unsigned int state3, unsigned int state4);
    unsigned int    TausStep(unsigned int &z, int s1, int s2, int s3, unsigned long M);
    unsigned int    LCGStep(unsigned int &z, unsigned int A, unsigned long C);
    double          getRandNum();
    double          HybridTaus();
    
    RNGSeeds        getState();
    
    // Overloaded operators for working with 'RNGSeeds' structures.
    RNG& operator=(const RNG &rhs);
    
    
private:
    // Used with the thread safe RNG to track state.
	unsigned int z1, z2, z3, z4;
};

#endif //__K_Wave_C__RNG__
