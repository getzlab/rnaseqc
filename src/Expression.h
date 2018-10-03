//
//  Expression.hpp
//  IntervalTree
//
//  Created by Aaron Graubert on 8/2/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#ifndef Expression_h
#define Expression_h

#include "Metrics.h"
#include "BamReader.h"
#include <vector>
#include <set>
#include <iostream>

//Utility functions
unsigned int extractBlocks(Alignment&, std::vector<Feature>&, chrom, bool);
//unsigned int legacyExtractBlocks(BamTools::BamAlignment&, std::vector<Feature>&, chrom);
std::list<Feature>* intersectBlock(Feature&, std::list<Feature>&);
void trimFeatures(Alignment&, std::list<Feature>&);
void trimFeatures(Alignment&, std::list<Feature>&, BaseCoverage&);
void dropFeatures(std::list<Feature>&, BaseCoverage&);

//Metrics functions
unsigned int fragmentSizeMetrics(unsigned int, std::map<chrom, std::list<Feature>>*, std::map<std::string, std::string>&, std::map<long long, unsigned long>&,std::vector<Feature>&, Alignment&, SeqLib::BamHeader&);

void exonAlignmentMetrics(unsigned int, std::map<chrom, std::list<Feature>>&, Metrics&, std::vector<Feature>&, Alignment&, SeqLib::BamHeader&, unsigned int, unsigned short, BaseCoverage&);

void legacyExonAlignmentMetrics(unsigned int, std::map<chrom, std::list<Feature>>&, Metrics&, std::vector<Feature>&, Alignment&, SeqLib::BamHeader&, unsigned int, unsigned short, BaseCoverage&);

#endif /* Expression_h */
