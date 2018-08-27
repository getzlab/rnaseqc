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

std::tuple<double, double, double> computeCoverage(std::ofstream&, const Feature&, const std::string&, const unsigned int, std::map<std::string, std::vector<CoverageEntry*> >&, std::list<double>&, BiasCounter&);

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
    tmp.transcript_id = exon.transcript_id;
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
        CoverageEntry tmp;
        tmp.offset = beg->offset;
        tmp.length = beg->length;
        tmp.transcript_id = beg->transcript_id;
        tmp.feature_id = beg->feature_id;
        this->coverage[gene_id].push_back(tmp);
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
    //in case of multiple transcripts per gene, support multiple TIDs
    //buffer all exons, and store an unordered set transcript list
    //then, after the gene has been exhausted into the buffer, compute coverage on each transcript
    std::unordered_set<std::string> transcripts;
    std::map<std::string, std::vector<CoverageEntry*>> coverage;
    for (auto entry = this->coverage[gene.feature_id].begin(); entry != this->coverage[gene.feature_id].end(); ++entry)
    {
        //memory allocation here
        //maybe switch to pointers?
        //prolly much smaller that way
        coverage[entry->feature_id].push_back(&(*entry));
        transcripts.insert(entry->transcript_id);
    }
    for (auto transcript = transcripts.begin(); transcript != transcripts.end(); ++transcript)
    {
        std::tuple<double, double, double> results = computeCoverage(this->writer, gene, *transcript, this->mask_size, coverage, this->exonCVs, this->bias);
        this->transcriptMeans.push_back(std::get<0>(results));
        this->transcriptStds.push_back(std::get<1>(results));
        this->transcriptCVs.push_back(std::get<2>(results));
    }
//    this->coverage[gene.feature_id].clear();
    this->coverage.erase(gene.feature_id);
    this->seen.insert(gene.feature_id);
}

void BaseCoverage::close()
{
    this->writer.flush();
    this->writer.close();
}

std::list<double>& BaseCoverage::getExonCVs()
{
    return this->exonCVs;
}

std::list<double>& BaseCoverage::getTranscriptMeans()
{
    return this->transcriptMeans;
}

std::list<double>& BaseCoverage::getTranscriptStds()
{
    return this->transcriptStds;
}

std::list<double>& BaseCoverage::getTranscriptCVs()
{
    return this->transcriptCVs;
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

std::tuple<double, double, double> computeCoverage(std::ofstream &writer, const Feature &gene, const std::string &transcript_id, const unsigned int mask_size, std::map<std::string, std::vector<CoverageEntry*> > &entries, std::list<double> &totalExonCV, BiasCounter &bias)
{
    std::vector<unsigned long> coverage;
    //    list<double> exonCV;
    std::vector<std::vector<bool> > coverageMask;
    unsigned int maskRemainder = mask_size;
    for (unsigned int i = 0; i < transcriptExons[transcript_id].size(); ++i)
    {
        coverageMask.push_back(std::vector<bool>(exonLengths[transcriptExons[transcript_id][i]], true));
        for (unsigned int j = 0; j < coverageMask.back().size() && maskRemainder; ++j, --maskRemainder)
            coverageMask.back()[j] = false;
    }
    maskRemainder = mask_size;
    for (int i = transcriptExons[transcript_id].size() - 1; i >= 0 && maskRemainder; --i)
        for (int j = coverageMask[i].size() - 1; j >= 0 && maskRemainder; --j, --maskRemainder)
            coverageMask[i][j] = false;
    for (unsigned int i = 0; i < transcriptExons[transcript_id].size(); ++i)
    {
        auto beg = entries[transcriptExons[transcript_id][i]].begin();
        auto end = entries[transcriptExons[transcript_id][i]].end();
        std::vector<unsigned long> exon_coverage(exonLengths[transcriptExons[transcript_id][i]], 0ul);
        while (beg != end)
        {
            add_range(exon_coverage, (**beg).offset, (**beg).length);
            ++beg;
        }
        double exonMean = 0.0, exonStd = 0.0, exonSize = 0.0;
        std::vector<bool> mask = coverageMask[i];
        
        //        assert(mask.size() == exon_coverage.size());
        for (unsigned int j = 0; j < mask.size(); ++j) if (mask[j]) exonSize += 1;
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
        coverage.reserve(coverage.size() + exon_coverage.size());
        coverage.insert(coverage.end(), exon_coverage.begin(), exon_coverage.end());
    }
    //    if (coverage.size() != transcriptCodingLengths[transcript_id])
    //        cerr << "Coding span mismatch " << transcript_id << " " << coverage.size() << " " << transcriptCodingLengths[transcript_id] << endl;
    //    auto current_cv = exonCV.begin();
    //    if (current_cv != exonCV.end())
    //    {
    //        double prev_cv = *(current_cv++);
    //        while (current_cv != exonCV.end())
    //        {
    //            deltaCV.push_back(fabs((*current_cv++) - prev_cv));
    //        }
    //    }
    
    bias.computeBias(gene, coverage);
    
    double avg = 0.0, std = 0.0;
    //    auto median = coverage.begin();
    if (mask_size)
    {
//        auto tmp_size = coverage.size();
        coverage.erase((mask_size > coverage.size() ? coverage.begin() : coverage.end() - mask_size), coverage.end());
//        tmp_size = coverage.size();
        if (coverage.size()) coverage.erase(coverage.begin(), (mask_size > coverage.size() ? coverage.end() : coverage.begin() + mask_size));
    }
    double size = static_cast<double>(coverage.size());
    if (size > 0)
    {
        for (auto beg = coverage.begin(); beg != coverage.end(); ++beg)
            avg += static_cast<double>(*beg) / size;
        for (auto base = coverage.begin(); base != coverage.end(); ++base)
            std += std::pow(static_cast<double>(*base) - avg, 2.0) / size;
        std = std::pow(std, 0.5);
        writer << gene.feature_id << "\t" << transcript_id << "\t";
        writer << avg << "\t" << std << "\t" << (std / avg) << std::endl;
    }
    return std::make_tuple(avg, std, (std / avg));
}
