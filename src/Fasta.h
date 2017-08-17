//
//  RefSeq.hpp
//  IntervalTree
//
//  Created by Aaron Graubert on 5/23/18.
//  Copyright © 2018 Aaron Graubert. All rights reserved.
//

#ifndef Fasta_h
#define Fasta_h

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <list>
#include <bioio.hpp>
#include <exception>

struct fileException : public std::exception {
    std::string error;
    fileException(std::string msg) : error(msg) {};
};

struct invalidContigException : public std::exception {
    std::string error;
    invalidContigException(std::string msg) : error(msg) {};
};

typedef long long coord;
typedef unsigned long indexType;
typedef unsigned short chrom;

static const double PAGE_SIZE = 1e6;
static const unsigned short CACHE_SIZE = 10u;

extern std::map<std::string, chrom> chromosomes;
chrom chromosomeMap(std::string);

class Fasta {
    bool isOpen;
    std::ifstream reader;
    std::unordered_map<indexType, std::string> pageCache;
    std::list<indexType> lru;
    std::unordered_map<chrom, bioio::FastaContigIndex> contigIndex;
    void updateLRU(indexType);
    indexType pageForContig(chrom);
    std::string readSeq(chrom, coord);
    unsigned long calls, misses;
public:
    Fasta() : isOpen(), reader(), pageCache(), lru(), contigIndex(), calls(), misses() {};
    ~Fasta();
    void open(std::string&);
    std::string getSeq(chrom, coord, coord);
    std::string getSeq(chrom, coord, coord, bool);
    indexType pageForCoord(chrom, coord);
    coord pageOffset(indexType);
    
};

double gc(std::string&);

#endif /* Fasta_h */
