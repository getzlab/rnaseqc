//
//  GTF.hpp
//  IntervalTree
//
//  Created by Aaron Graubert on 6/28/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#ifndef GTF_h
#define GTF_h

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <vector>
#include <sstream>
#include "Fasta.h"

namespace rnaseqc {
    struct gtfException : public std::exception {
        std::string error;
        gtfException(std::string msg) : error(msg) {};
    };
    
    enum FeatureType {Gene, Transcript, Exon, Other};
    
    struct Feature {
        //Represents arbitrary genome features
        coord start, end;
        chrom chromosome;
        Strand strand;
        FeatureType type;
        std::string feature_id, gene_id, transcript_type;
        bool ribosomal;
    };
    
    //For comparing features
    bool operator==(const Feature &a, const Feature &b);
    bool compIntervalStart(const Feature&, const Feature&);
    bool compIntervalEnd(const Feature&, const Feature&);
    bool intersectPoint(const Feature&, const coord);
    bool intersectInterval(const Feature&, const Feature&);
    int partialIntersect(const Feature&, const Feature&);

    struct FeatureSpan {
        chrom chromosome;
        coord start, length;
    };
    
    
    extern std::map<std::string, std::string> geneNames, geneSeqs;
extern std::map<std::string, coord> geneLengths, geneCodingLengths;
    extern std::map<std::string, FeatureSpan> exonLengths;
    extern std::vector<std::string> geneList, exonList;
    extern std::map<std::string, std::vector<std::string>> exonsForGene;
    
    std::ifstream& operator>>(std::ifstream&, Feature&);
    std::map<std::string,std::string>& parseAttributes(std::string&, std::map<std::string,std::string>&);
}

#endif /* GTF_h */
