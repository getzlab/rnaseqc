//
//  Metrics.cpp
//  IntervalTree
//
//  Created by Aaron Graubert on 7/5/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#include "Metrics.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include <unordered_set>

std::map<std::string, double> geneCounts, exonCounts; //counters for read coverage of genes and exons

std::tuple<double, double, double> computeCoverage(std::ofstream&, const Feature&, const unsigned int, const std::map<std::string, std::vector<unsigned long> >&, std::list<double>&, BiasCounter&);

void add_range(std::vector<unsigned long>&, coord, unsigned int);

void Metrics::increment(std::string key)
{
    this->counter[key]++;
}

void Metrics::increment(std::string key, int n)
{
    this->counter[key] += n;
}

unsigned long Metrics::get(std::string key)
{
    return this->counter[key];
}

double Metrics::frac(std::string a, std::string b)
{
    return static_cast<double>(this->get(a)) / this->get(b);
}

std::ofstream& operator<<(std::ofstream &stream, Metrics &counter)
{
    std::vector<std::string> keys =  {
        //"Alternative Alignments",
        //"Chimeric Reads",
        "Duplicate Reads",
        "End 1 Antisense",
        "End 2 Antisense",
        "End 1 Bases",
        "End 2 Bases",
        "End 1 Mapped Reads",
        "End 2 Mapped Reads",
        "End 1 Mismatches",
        "End 2 Mismatches",
        "End 1 Sense",
        "End 2 Sense",
        "Exonic Reads", //reads
        "Failed Vendor QC",
        // "Filtered by tag: ", //FIXME
        "Intergenic Reads",
        "Intragenic Reads",
        "Intron/Exon Disqualified Reads", //capitalize
        "Intronic Reads",
        "Low quality reads",
        "Mapped Duplicate Reads",
        "Mapped Reads",
        "Mapped Unique Reads",
        "Mismatched Bases", //bases
        "Reads excluded from exon counts",
        "Reads used for Intron/Exon counts",
        "rRNA Reads",
        "Split Reads",
        "Total Bases",
        "Total Mapped Pairs",
        "Total Reads",
        "Unique Mapping, Vendor QC Passed Reads",
        "Unpaired Reads"
    };
    stream << "Alternative Alignments\t" << counter.get("Alternative Alignments") << std::endl;
    stream << "Chimeric Reads\t";
    if (counter.get("Chimeric Reads_tag")) stream << counter.get("Chimeric Reads_tag") << std::endl;
    else stream << counter.get("Chimeric Reads_contig") << std::endl;
    for (int i = 0; i < keys.size(); ++i) stream << keys[i] << "\t" << counter.get(keys[i]) << std::endl;
    auto beg = counter.counter.begin();
    auto end = counter.counter.end();
    while (beg != end)
    {
        if( beg->first.length() > 17 && beg->first.substr(0,17) == "Filtered by tag: ")
        {
            stream << beg->first << "\t" << beg->second << std::endl;
        }
        ++beg;
    }
    return stream;
}

void Collector::add(const std::string &gene_id, const std::string &exon_id, const double coverage)
{
    if (coverage > 0)
    {
        this->data[gene_id].push_back(std::pair<std::string, double>(exon_id, coverage));
        this->dirty = true;
    }
}

void Collector::collect(const std::string &gene_id)
{
    for (auto entry = this->data[gene_id].begin(); entry != this->data[gene_id].end(); ++entry)
    {
        (*this->target)[entry->first] += entry->second;
        this->total += entry->second;
    }
}

void Collector::collectSingle(const std::string &gene_id)
{
    for (auto entry = this->data[gene_id].begin(); entry != this->data[gene_id].end(); ++entry)
    {
        (*this->target)[entry->first] += 1.0;
    }
}

bool Collector::queryGene(const std::string &gene_id)
{
    return static_cast<bool>(this->data[gene_id].size());
}

bool Collector::isDirty()
{
    return this->dirty;
}

double Collector::sum()
{
    return this->total;
}

void BaseCoverage::add(const Feature &exon, const coord start, const coord end) //Adds to the cache
{
    CoverageEntry tmp;
    tmp.offset = start - exon.start;
    tmp.length = end - start;
    tmp.feature_id = exon.feature_id;
    this->cache[exon.gene_id].push_back(tmp);
}

void BaseCoverage::commit(const std::string &gene_id) //moves one gene out of the cache to permanent storage
{
    if (this->seen.count(gene_id))
    {
        std::cerr << "Gene encountered after computing coverage " << gene_id << std::endl;
        return;
    }
    auto beg = this->cache[gene_id].begin();
    auto end = this->cache[gene_id].end();
    while (beg != end)
    {
        if (this->coverage.find(beg->feature_id) == this->coverage.end()) this->coverage[beg->feature_id] = std::vector<unsigned long>(exonLengths[beg->feature_id], 0ul);
        add_range(this->coverage[beg->feature_id], beg->offset, beg->length);
        ++beg;
    }
}

void BaseCoverage::reset() //Empties the cache
{
    this->cache.clear();
}

//void BaseCoverage::clearCoverage()
//{
//    std::cerr << "Force dump of " << this->coverage.size() << " genes (" << this->coverage.begin()->first << ")" <<std::endl;
//    this->coverage.clear();
//}

void BaseCoverage::compute(const Feature &gene) //computes per-base coverage of the gene
{
    //Coverage is stored in EID -> coverage vector
    //First iterate over all exons of the gene and ensure they're filled
    //That way, stiching the exons will result in a complete transcript even for exons which haven't been seen
    for (auto exon_id = exonsForGene[gene.feature_id].begin(); exon_id != exonsForGene[gene.feature_id].end(); ++exon_id)
        if (this->coverage.find(*exon_id) == this->coverage.end()) this->coverage[*exon_id] = std::vector<unsigned long>(exonLengths[*exon_id], 0ul);
    //then compute coverage for the gene
    std::tuple<double, double, double> results = computeCoverage(this->writer, gene, this->mask_size, this->coverage, this->exonCVs, this->bias);
    this->geneMeans.push_back(std::get<0>(results));
    this->geneStds.push_back(std::get<1>(results));
    this->geneCVs.push_back(std::get<2>(results));
    //Now clean out the coverage map to save memory
    for (auto exon_id = exonsForGene[gene.feature_id].begin(); exon_id != exonsForGene[gene.feature_id].end(); ++exon_id)
         this->coverage.erase(*exon_id);
    this->seen.insert(gene.feature_id);
}

void BaseCoverage::close()
{
    this->writer.flush();
    this->writer.close();
}

void BiasCounter::computeBias(const Feature &gene, std::vector<unsigned long> &coverage)
{
    if (coverage.size() < this->geneLength) return; //Must meet minimum length req
    double lcov = 0.0, rcov = 0.0, windowSize = static_cast<double>(this->windowSize);
    for (unsigned int i = this->offset; i < this->offset + this->windowSize && i < coverage.size(); ++i)
        lcov += static_cast<double>(coverage[i]) / windowSize;
    for (int i = coverage.size() - (1 + this->offset); i <= 0 && i > coverage.size() - (this->offset + this->windowSize); --i)
        rcov += static_cast<double>(coverage[i]) / windowSize;
    this->threeEnd[gene.feature_id] += gene.strand == 1 ? rcov : lcov;
    this->fiveEnd[gene.feature_id] += gene.strand == 1 ? lcov : rcov;
}

double BiasCounter::getBias(const std::string &geneID)
{
    double cov5 = this->fiveEnd[geneID];
    double cov3 = this->threeEnd[geneID];
    if (cov5 + cov3 > 0.0) return cov3 / (cov5 + cov3);
    return -1.0;
}


void add_range(std::vector<unsigned long> &coverage, coord offset, unsigned int length)
{
    //    if (offset + length >= coverage.size()) coverage.resize(offset + length, 0ul);
    //    for (coord i = offset; offset < offset + length; ++i) coverage[i] = coverage[i] + 1;
    //    for (coord i = 0; i < offset + length; ++i)
    //    {
    //        unsigned long x = i >= offset ? 1ul : 0ul;
    //        if (i >= coverage.size()) coverage.push_back(x);
    //        else coverage[i] = coverage[i] + x;
    //    }
    for (coord i = offset; i < offset + length; ++i) coverage[i] += 1ul;
}

std::tuple<double, double, double> computeCoverage(std::ofstream &writer, const Feature &gene, const unsigned int mask_size, const std::map<std::string, std::vector<unsigned long> > &coverage, std::list<double> &totalExonCV, BiasCounter &bias)
{
    std::vector<std::vector<bool> > coverageMask;
    std::vector<unsigned long> geneCoverage;
    unsigned int maskRemainder = mask_size;
    for (unsigned int i = 0; i < exonsForGene[gene.feature_id].size(); ++i)
    {
        coverageMask.push_back(std::vector<bool>(exonLengths[exonsForGene[gene.feature_id][i]], true)); //First store a pre-filled mask for the exon
        for (unsigned int j = 0; j < coverageMask.back().size() && maskRemainder; ++j, --maskRemainder) //now, remove coverage from the front of the exon until either it, or the mask size is depleted
            coverageMask.back()[j] = false;
    }
    maskRemainder = mask_size; //reset the exon mask to mask out the end
    for (int i = exonsForGene[gene.feature_id].size() - 1; i >= 0 && maskRemainder; --i) //repeat the process, masking out regions from the back until the mask size is depleted
        for (int j = coverageMask[i].size() - 1; j >= 0 && maskRemainder; --j, --maskRemainder)
            coverageMask[i][j] = false;
    for (unsigned int i = 0; i < exonsForGene[gene.feature_id].size(); ++i)
    {
        const std::vector<unsigned long> &exon_coverage = coverage.at(exonsForGene[gene.feature_id][i]); //get the coverage vector for the current exon
        double exonMean = 0.0, exonStd = 0.0, exonSize = 0.0;
        std::vector<bool> mask = coverageMask[i];
        
//        assert(mask.size() == exon_coverage.size());
//        assert(exon_coverage.size() == exonLengths[exonsForGene[gene.feature_id][i]]);

        for (unsigned int j = 0; j < mask.size(); ++j) if (mask[j]) exonSize += 1.0; //count the remaining unmasked length of the exon
        if (exonSize > 0)
        {
//            std::list<unsigned int> effective_coverage;
            auto maskIter = mask.begin();
            for (auto start = exon_coverage.begin(); start != exon_coverage.end(); ++start)
                if (*(maskIter++))
                {
                    exonMean += static_cast<double>(*start) / exonSize;
//                    effective_coverage.push_back(*start);
                }
            maskIter = mask.begin();
            for (auto start = exon_coverage.begin(); start != exon_coverage.end(); ++start)
                if (*(maskIter++)) exonStd += pow(static_cast<double>(*start) - exonMean, 2.0) / exonSize;
            exonStd = pow(exonStd, 0.5);// / exonMean; //technically it's a CV now
//            effective_coverage.sort();
            //            assert(effective_coverage.size() == exonSize);
            exonStd /= exonMean; //now it's a CV
            if (!(std::isnan(exonStd) || std::isinf(exonStd)))
            {
                totalExonCV.push_back(exonStd);
                //            exonCV.push_back(log2(exonStd));
            }
//            biasCoverage.reserve(biasCoverage.size() + exon_coverage.size());
//            biasCoverage.insert(biasCoverage.end(), exon_coverage.begin(), exon_coverage.end());
        }
        geneCoverage.reserve(geneCoverage.size() + exon_coverage.size());
        geneCoverage.insert(geneCoverage.end(), exon_coverage.begin(), exon_coverage.end());
    }
//    assert(geneCoverage.size() == geneCodingLengths[gene.feature_id]);
    if (geneCoverage.size() != geneCodingLengths[gene.feature_id])
        std::cerr << "Coding span mismatch " << gene.feature_id << " " << geneCoverage.size() << " " << geneCodingLengths[gene.feature_id] << std::endl;
    
    if (geneCounts[gene.feature_id] > bias.getThreshold()) bias.computeBias(gene, geneCoverage);
    
    double avg = 0.0, std = 0.0;
    //    auto median = coverage.begin();
    if (mask_size)
    {
//        auto tmp_size = coverage.size();
        //to account for the mask, erase bases from the vector
        //If the mask is larger than the gene, just erase all of it
        //Otherwise, erase from (end-mask) -> end
        geneCoverage.erase((mask_size > geneCoverage.size() ? geneCoverage.begin() : geneCoverage.end() - mask_size), geneCoverage.end());
//        tmp_size = coverage.size();
        // If there is still coverage area, erase from the front to either the end (if the mask is larger than remaining coverage) or to (front + mask)
        if (geneCoverage.size()) geneCoverage.erase(geneCoverage.begin(), (mask_size > geneCoverage.size() ? geneCoverage.end() : geneCoverage.begin() + mask_size));
    }
    double size = static_cast<double>(geneCoverage.size());
    if (size > 0) //If there's still any coverage after applying the mask
    {
        for (auto beg = geneCoverage.begin(); beg != geneCoverage.end(); ++beg)
            avg += static_cast<double>(*beg) / size;
        for (auto base = geneCoverage.begin(); base != geneCoverage.end(); ++base)
            std += std::pow(static_cast<double>(*base) - avg, 2.0) / size;
        std = std::pow(std, 0.5);
        writer << gene.feature_id << "\t";
        writer << avg << "\t" << std << "\t" << (std / avg) << std::endl;
    }
    return std::make_tuple(avg, std, (std / avg));
}
