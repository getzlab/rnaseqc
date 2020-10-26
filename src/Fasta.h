//
//  Fasta.h
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

namespace rnaseqc {
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
    
    static const double PAGE_SIZE = 1e6; // Size of each cache page (in bases)
    static const unsigned short CACHE_SIZE = 10u; // How many pages are stored in the cache
    
    extern std::map<std::string, chrom> chromosomes;
    
    enum Strand {Forward, Reverse, Unknown};
    chrom chromosomeMap(std::string);
    
    class Fasta {
        // Represents an entire fasta file
        // Uses the bioio library for quickly retrieving sequences
        // Uses an internal LRU cache to minimize required reading
        bool _open;
        std::ifstream reader;
        std::unordered_map<indexType, std::string> pageCache;
        std::list<indexType> lru;
        std::unordered_map<chrom, bioio::FastaContigIndex> contigIndex;
        void updateLRU(indexType);
        indexType pageForContig(chrom);
        std::string readSeq(chrom, coord);
        unsigned long calls, misses;
    public:
        Fasta() : _open(), reader(), pageCache(), lru(), contigIndex(), calls(), misses() {};
        ~Fasta();
        void open(std::string&);
        std::string getSeq(chrom, coord, coord);
        std::string getSeq(chrom, coord, coord, Strand);
        indexType pageForCoord(chrom, coord);
        coord pageOffset(indexType);
        bool isOpen() const;
        
    };
    
    double gc(std::string&);
}

#endif /* Fasta_h */
