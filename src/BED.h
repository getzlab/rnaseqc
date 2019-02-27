//
//  BED.hpp
//  IntervalTree
//
//  Created by Aaron Graubert on 7/11/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#ifndef BED_h
#define BED_h

#include "GTF.h"

namespace rnaseqc {
    struct bedException : public std::exception {
        std::string error;
        bedException(std::string msg) : error(msg) {};
    };
    
    std::ifstream& extractBED(std::ifstream&, Feature&);
}
#endif /* BED_h */
