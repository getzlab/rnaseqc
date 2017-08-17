//
//  Metrics.hpp
//  IntervalTree
//
//  Created by Aaron Graubert on 7/5/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#ifndef Metrics_h
#define Metrics_h

#include <map>
#include <fstream>
#include <string>
#include <vector>
#include <utility>

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

std::ofstream& operator<<(std::ofstream&, Metrics&);
#endif /* Metrics_h */
