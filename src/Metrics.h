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
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <tuple>
#include <list>
#include <unordered_set>
#include <iterator>

namespace rnaseqc {
    class Metrics;
}

std::ofstream& operator<<(std::ofstream&, rnaseqc::Metrics&);

namespace rnaseqc {
    class Metrics {
        // For storing arbitrary counters
        std::map<std::string, unsigned long> counter;
    public:
        Metrics() : counter(){};
        void increment(std::string);
        void increment(std::string, int);
        unsigned long get(std::string);
        double frac(std::string, std::string);
        friend std::ofstream& ::operator<<(std::ofstream&, Metrics&);
    };
    
    class Collector {
        // For temporarily holding coverage on a read before we're ready to commit that coverage to a gene
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
        // Represents a single segment of aligned read bases for base-coverage computation
        coord offset;
        unsigned int length;
        std::string feature_id;
    };
    
    class BiasCounter {
        // For counting 3'/5' bias coverage
        const int offset;
        const int windowSize;
        const unsigned long geneLength;
        const unsigned int detectionThreshold;
        unsigned int countedGenes;
        std::map<std::string, unsigned long> fiveEnd;
        std::map<std::string, unsigned long> threeEnd;
    public:
        BiasCounter(int offset, int windowSize, unsigned long geneLength, unsigned int detectionThreshold) : offset(offset), windowSize(windowSize), geneLength(geneLength), detectionThreshold(detectionThreshold), countedGenes(0), fiveEnd(), threeEnd()
        {
            
        }
        
        void computeBias(const Feature&, std::vector<unsigned long>&);
        unsigned int countGenes() const;
        double getBias(const std::string&);
        const unsigned int getThreshold() const {
            return this->detectionThreshold;
        }
    };
    
    class BaseCoverage {
        // For computing per-base coverage of genes
        std::map<std::string, std::vector<CoverageEntry> > cache; //GID -> Entry<EID> tmp cache as exon hits are recorded
        std::map<std::string, std::vector<unsigned long> > coverage; //EID -> Coverage vector for exons still in window
        std::ofstream writer;
        const unsigned int mask_size;
        std::list<double> exonCVs, geneMeans, geneStds, geneCVs;
        BiasCounter &bias;
        std::unordered_set<std::string> seen;
        BaseCoverage(const BaseCoverage&) = delete; //No!
    public:
        BaseCoverage(const std::string &filename, const unsigned int mask, bool openFile, BiasCounter &biasCounter) : coverage(), cache(), writer(openFile ? filename : "/dev/null"), mask_size(mask), exonCVs(), geneMeans(), geneStds(), geneCVs(), bias(biasCounter), seen()
        {
            if ((!this->writer.is_open()) && openFile) throw std::runtime_error("Unable to open BaseCoverage output file");
            this->writer << "gene_id\tcoverage_mean\tcoverage_std\tcoverage_CV" << std::endl;
        }
        
        void add(const Feature&, const coord, const coord); //Adds to the cache
        void commit(const std::string&); //moves one gene out of the cache and adds hits to exon coverage vector
        void reset(); //Empties the cache
        //    void clearCoverage(); //empties out data that won't be used
        void compute(const Feature&); //Computes the per-base coverage for all transcripts in the gene
        void close(); //Flush and close the ofstream
        BiasCounter& getBiasCounter() const {
            return this->bias;
        }
        std::list<double>& getExonCVs() {
            return this->exonCVs;
        }
        std::list<double>& getGeneMeans() {
            return this->geneMeans;
        }
        std::list<double>& getGeneStds() {
            return this->geneStds;
        }
        std::list<double>& getGeneCVs() {
            return this->geneCVs;
        }
    };

    template <typename T> void sortContainer(T &data) {
        std::sort(data.begin(), data.end());
    }

    template <typename T> void sortContainer(std::list<T> &data) {
        data.sort();
    }
    
    template <typename T> double computeMedian(unsigned long size, T &&iterator)
    {
        if (size <= 0) // Couldn't decide if it would make sense to just report a median of 0. This seemed safer
            throw std::range_error("Cannot compute median of an empty list");
        for (unsigned long midpoint = (size - 1) / 2; midpoint > 0; --midpoint) ++iterator;
        if (size % 2)
        {
            double value = static_cast<double>(*(iterator++));
            return (value + static_cast<double>(*iterator)) / 2.0;
        }
        return static_cast<double>(*iterator);
    }

    typedef std::tuple<double, double, double, double> statsTuple;
    
    enum StatIdx {avg = 0, med = 1, std = 2, mad = 3};
    
    template <typename T>
    statsTuple getStatistics(T &data) {
        if (data.size()) {
            double avg = 0.0, std = 0.0;
            std::vector<double> deviations;
            sortContainer(data);
            const double size = data.size();
            double median = computeMedian(size, data.begin());
            for (auto element = data.begin(); element != data.end(); ++element) {
                avg += static_cast<double>(*element) / size;
                deviations.push_back(fabs(static_cast<double>(*element) - median));
            }
            sortContainer(deviations);
            double medDev = computeMedian(deviations.size(), deviations.begin()) * 1.4826;
            for (auto element = data.begin(); element != data.end(); ++element)
                std += pow(static_cast<double>(*element) - avg, 2.0) / size;
            std = pow(std, 0.5);
            return statsTuple(avg, median, std, medDev);
        }
        return statsTuple(NAN, NAN, NAN, NAN);
    }
    
    extern std::map<std::string, double> uniqueGeneCounts, geneCounts, exonCounts, geneFragmentCounts; //counters for read coverage of genes and exons
    extern std::map<std::string, std::unordered_set<std::string> > fragmentTracker; // tracks fragments encountered by each gene
}

#endif /* Metrics_h */
