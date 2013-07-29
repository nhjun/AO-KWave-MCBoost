//
//  multikey.h
//  K_Wave_C
//
//  Created by jacob on 07/28/13.
//  Copyright (c) 2013 BMPI. All rights reserved.
//


#ifndef _MULTIKEY_H_
#define _MULTIKEY_H_

#include <stdint.h>


class MultiKey
{
public:
    uint32_t key1;
    uint32_t key2;
    uint32_t key3;
    uint32_t key4;
    
    
    MultiKey(uint32_t k1, uint32_t k2, uint32_t k3, uint32_t k4);
    ~MultiKey();
        
    
    bool operator<(const MultiKey &right) const
    {
        if ( key1 == right.key1 )
        {
            if ( key2 == right.key2 )
            {
                if ( key3 == right.key3 )
                {
                    return key4 < right.key4;
                }
                else
                {
                    return key3 < right.key3;
                    
                }
             }
             else
             {
                 return key2 < right.key2;
             }
        }
        else
        {
            return key1 < right.key1;
            
        }
    }
                    

                    
};

#endif // _MULTIKEY_H_


