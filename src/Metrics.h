//
//  Metrics.hpp
//  IntervalTree
//
//  Created by Aaron Graubert on 7/5/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#ifndef Metrics_h
#define Metrics_h

#include "GTF.h"
#include <map>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <exception>

class Metrics {
    std::map<std::string, unsigned long> counter;
public:
    Metrics() : counter(){};
    void increment(std::string);
    void increment(std::string, int);
    unsigned long get(std::string);
    double frac(std::string, std::string);
    friend std::ofstream& operator<<(std::ofstream&, Metrics&);
};

class Collector {
    std::map<std::string, std::vector<std::pair<std::string, double> > > data;
    std::map<std::string, double> *target;
    bool dirty;
    double total;
public:
    Collector(std::map<std::string, double> *dataTarget) : data(), target(dataTarget), dirty(false), total(0.0)
    {
        
    }
    void add(const std::string&, const std::string&, const double);
    void collect(const std::string&);
    void collectSingle(const std::string&); //for legacy exon detection
    bool queryGene(const std::string&);
    bool isDirty();
    double sum();
};

struct CoverageEntry {
    coord offset;
    unsigned int length;
    std::string transcript_id, feature_id;
};

class BaseCoverage {
    std::map<std::string, std::vector<CoverageEntry> > coverage, cache;
    //Coverage is EID -> Entry<GID>
    //Cache is GID -> Entry<EID>
    std::ofstream writer;
public:
    BaseCoverage(const std::string &filename) : coverage(), cache(), writer(filename)
    {
        if (!this->writer) throw std::runtime_error("Unable to open BaseCoverage output file");
    }
    
    void add(const Feature&, const coord, const coord); //Adds to the cache
    void commit(const std::string&); //moves one gene out of the cache to permanent storage
    void reset(); //Empties the cache
    void dump(const Feature&); //Dumps one exon to the tmp file
    void close(); //Flush and close the ofstream
};

class BiasCounter {
    const int offset;
    const int windowSize;
    const unsigned long geneLength;
    std::map<std::string, unsigned long> fiveEnd;
    std::map<std::string, unsigned long> threeEnd;
public:
    BiasCounter(int offset, int windowSize, unsigned long geneLength) : offset(offset), windowSize(windowSize), geneLength(geneLength), fiveEnd(), threeEnd()
    {
        
    }
    
    void checkBias(Feature&, Feature&);
    double getBias(const std::string&);
};


std::ofstream& operator<<(std::ofstream&, Metrics&);
#endif /* Metrics_h */
