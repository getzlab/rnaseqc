//
//  Expression.hpp
//  IntervalTree
//
//  Created by Aaron Graubert on 8/2/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#ifndef Expression_h
#define Expression_h

#include "GTF.h"
#include "Metrics.h"
#include <vector>
#include <list>
#include <set>
#include <iostream>
#include <bamtools/api/BamReader.h>
#include <bamtools/api/BamAlignment.h>

//Utility functions
template <typename T>
bool isIn(std::set<T> &s, T member)
{
    return s.size() ? s.find(member) != s.end() : false;
}

unsigned int extractBlocks(BamTools::BamAlignment&, std::vector<Feature>&, unsigned short);
//unsigned int legacyExtractBlocks(BamTools::BamAlignment&, std::vector<Feature>&, unsigned short);
std::list<Feature>* intersectBlock(Feature&, std::list<Feature>&);
void trimFeatures(BamTools::BamAlignment&, std::list<Feature> &);

//Metrics functions
unsigned int fragmentSizeMetrics(unsigned int, std::map<unsigned short, std::list<Feature>>*, std::map<std::string, std::string>&, std::list<long long>&, BamTools::SamSequenceDictionary&, std::vector<Feature>&, BamTools::BamAlignment&);

void exonAlignmentMetrics(unsigned int, std::map<unsigned short, std::list<Feature>>&, Metrics&, BamTools::SamSequenceDictionary&, std::map<std::string, double>&, std::map<std::string, double>&, std::vector<Feature>&, BamTools::BamAlignment&, unsigned int);

void legacyExonAlignmentMetrics(unsigned int, std::map<unsigned short, std::list<Feature>>&, Metrics&, BamTools::SamSequenceDictionary&, std::map<std::string, double>&, std::map<std::string, double>&, std::vector<Feature>&, BamTools::BamAlignment&, unsigned int);

#endif /* Expression_h */
