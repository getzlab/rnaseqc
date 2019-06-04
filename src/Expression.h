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

namespace rnaseqc {
    //Utility functions
    unsigned int extractBlocks(Alignment&, std::vector<Feature>&, chrom, bool);
    //unsigned int legacyExtractBlocks(BamTools::BamAlignment&, std::vector<Feature>&, chrom);
    std::list<Feature>* intersectBlock(Feature&, std::list<Feature>&);
    void trimFeatures(Alignment&, std::list<Feature>&);
    void trimFeatures(Alignment&, std::list<Feature>&, BaseCoverage&);
    void dropFeatures(std::list<Feature>&, BaseCoverage&);
    
    // Definitions for fragment tracking
    typedef std::tuple<std::string, coord> FragmentMateEntry; // Used to record mate end point
    const std::size_t EXON = 0, ENDPOS = 1;
    
    //Metrics functions
    unsigned int fragmentSizeMetrics(unsigned int, std::map<chrom, std::list<Feature>>*, std::map<std::string, FragmentMateEntry>&, std::map<long long, unsigned long>&,std::vector<Feature>&, Alignment&, SeqLib::HeaderSequenceVector&);
    
    void exonAlignmentMetrics(std::map<chrom, std::list<Feature>>&, Metrics&, std::vector<Feature>&, Alignment&, SeqLib::HeaderSequenceVector&, unsigned int, Strand, BaseCoverage&, const bool, const bool);
    
    void legacyExonAlignmentMetrics(unsigned int, std::map<chrom, std::list<Feature>>&, Metrics&, std::vector<Feature>&, Alignment&, SeqLib::HeaderSequenceVector&, unsigned int, Strand, BaseCoverage&, const bool, const bool);
    
    Strand feature_strand(Alignment&, Strand);
}

#endif /* Expression_h */
