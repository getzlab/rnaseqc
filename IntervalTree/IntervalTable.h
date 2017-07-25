//
//  IntervalTable.hpp
//  IntervalTree
//
//  Created by Aaron Graubert on 7/12/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#ifndef IntervalTable_h
#define IntervalTable_h

#include "Interval.h"
#include <map>
#include <list>
#include <set>
#include <vector>

typedef std::vector<std::list<Feature*>> featureMap;

class IntervalTable
{
    std::list<Feature> features;
    
public:
    IntervalTable();
    void insert(Feature&);
    featureMap* loadRange(coord, unsigned int);
    featureMap* loadRange();
    resultSet* intersect(Feature&);
};

resultSet* intersectFeature(Feature&, featureMap*);
#endif /* IntervalTable_h */
