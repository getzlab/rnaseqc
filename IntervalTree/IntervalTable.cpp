//
//  IntervalTable.cpp
//  IntervalTree
//
//  Created by Aaron Graubert on 7/12/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#include "IntervalTable.h"

IntervalTable::IntervalTable() : features()
{
}

void IntervalTable::insert(Feature &next)
{
    this->features.push_back(next);
}

featureMap* IntervalTable::loadRange(coord start, unsigned int window)
{
    featureMap *output = new featureMap();
    output->reserve(window+1);
    std::list<Feature*> featuresInWindow;
    auto current = this->features.begin();
    for (unsigned int i = 0u; i<= window; ++i)
    {
        output->push_back(std::list<Feature*>());
        while (current->start >= start + i && current != this->features.end()) featuresInWindow.push_back(&(*current++));
        auto f = featuresInWindow.begin();
        while (f != featuresInWindow.end())
        {
            if ((**f).end > start + i) f = featuresInWindow.erase(f);
            else
            {
                output->rbegin()->push_back(*f);
            }
        }
    }
    return output;
}

featureMap* IntervalTable::loadRange()
{
    unsigned int endpoint = this->features.rbegin()->end;
    featureMap *output = new featureMap();
    output->reserve(endpoint+1);
    std::list<Feature*> featuresInWindow;
    auto current = this->features.begin();
    for (unsigned int i = 0u; i<= endpoint; ++i)
    {
        output->push_back(std::list<Feature*>());
        while (current->start >= i && current != this->features.end()) featuresInWindow.push_back(&(*current++));
        auto f = featuresInWindow.begin();
        while (f != featuresInWindow.end())
        {
            if ((**f).end > i) f = featuresInWindow.erase(f);
            else
            {
                output->rbegin()->push_back(*f);
            }
        }
    }
    return output;
}

resultSet* IntervalTable::intersect(Feature &target)
{
    auto range = this->loadRange(target.start, 20000u);
    resultSet* output =  intersectFeature(target, range);
    delete range;
    return output;
}

resultSet* intersectFeature(Feature &target, featureMap *range)
{
    resultSet *output = new resultSet(compIntervalStart);
    for (coord i = target.start; i <= target.end; ++i)
    {
        std::list<Feature*> index = range->at(i);
        if (index.size())
        {
//            while (current.first != current.second) output->insert((current.first++)->second);
            for (auto current = index.begin(); current != index.end(); ++current) output->insert(**current);
        }
    }
    return output;
}
