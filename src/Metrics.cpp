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
#include <algorithm>
#include <iterator>

namespace rnaseqc {


    std::map<std::string, double> uniqueGeneCounts, geneCounts, exonCounts, geneFragmentCounts; //counters for read coverage of genes and exons

    std::map<std::string, std::unordered_set<std::string> > fragmentTracker; // tracks fragments encountered by each gene
    
    std::tuple<double, double, double> computeCoverage(Fasta&, std::ofstream&, const Feature&, const unsigned int, const std::map<std::string, std::vector<unsigned long> >&, std::map<std::string, ExonCoverage>&, BiasCounter&);

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

    // Add coverage to an exon
    void Collector::add(const std::string &gene_id, const std::string &exon_id, const double coverage)
    {
        if (coverage > 0)
        {
            this->data[gene_id].push_back(std::pair<std::string, double>(exon_id, coverage));
            this->dirty = true;
        }
    }

    //Commit all the exon coverage from this gene to the global exon coverage counter
    void Collector::collect(const std::string &gene_id)
    {
        for (auto entry = this->data[gene_id].begin(); entry != this->data[gene_id].end(); ++entry)
        {
            (*this->target)[entry->first] += entry->second;
            this->total += entry->second;
        }
    }

    //Legacy version of the above function. Ignores the actual coverage and reports a full read count
    void Collector::collectSingle(const std::string &gene_id)
    {
        for (auto entry = this->data[gene_id].begin(); entry != this->data[gene_id].end(); ++entry)
        {
            (*this->target)[entry->first] += 1.0;
        }
    }

    //Check if there is any coverage on any exon of this gene
    bool Collector::queryGene(const std::string &gene_id)
    {
        return static_cast<bool>(this->data[gene_id].size());
    }

    // Check if any coverage has been reported whatsoever
    bool Collector::isDirty()
    {
        return this->dirty;
    }

    //Get the sum of all coverage that was committed for this read (should always be <= 1)
    double Collector::sum()
    {
        return this->total;
    }

    //Adds coverage from one aligned segment of a read to this exon. Coverage feeds into cache until gene leaves search window
    void BaseCoverage::add(const Feature &exon, const coord start, const coord end)
    {
        CoverageEntry tmp;
        tmp.offset = start - exon.start;
        tmp.length = end - start;
        tmp.feature_id = exon.feature_id;
        this->cache[exon.gene_id].push_back(tmp);
    }

    //Commit the cached coverage to this gene after deciding to count the read towards the gene
    void BaseCoverage::commit(const std::string &gene_id)
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
            if (this->coverage.find(beg->feature_id) == this->coverage.end()) this->coverage[beg->feature_id] = std::vector<unsigned long>(exonLengths[beg->feature_id].length, 0ul);
            //Add each coverage entry to the per-base coverage vector for the exon
            //At this stage exons each have their own vectors.
            //During the compute() step, exons get stiched together
            add_range(this->coverage[beg->feature_id], beg->offset, beg->length);
            ++beg;
        }
    }

    void BaseCoverage::reset() //Empties the cache
    {
        this->cache.clear();
    }

    //computes per-base coverage of the gene
    void BaseCoverage::compute(const Feature &gene)
    {
        //Coverage is stored in EID -> coverage vector
        //First iterate over all exons of the gene and ensure they're filled
        //That way, stiching the exons will result in a complete transcript even for exons which haven't been seen
        for (auto exon_id = exonsForGene[gene.feature_id].begin(); exon_id != exonsForGene[gene.feature_id].end(); ++exon_id)
            if (this->coverage.find(*exon_id) == this->coverage.end()) this->coverage[*exon_id] = std::vector<unsigned long>(exonLengths[*exon_id].length, 0ul);
        //then compute coverage for the gene
        std::tuple<double, double, double> results = computeCoverage(this->fastaReader, this->writer, gene, this->mask_size, this->coverage, this->exonCoverage, this->bias);
        if (std::get<0>(results) != -1)
        {
            this->geneMeans.push_back(std::get<0>(results));
            this->geneStds.push_back(std::get<1>(results));
            this->geneCVs.push_back(std::get<2>(results));
        }
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

    //Compute 3'/5' bias based on genes' per-base coverage
    void BiasCounter::computeBias(const Feature &gene, std::vector<unsigned long> &coverage)
    {

        if (coverage.size() < this->geneLength) return; //Must meet minimum length req
        unsigned long peak = 0ul;
        unsigned peak_pos = 0;
        for (unsigned i = 0; i < coverage.size(); ++i) if (coverage[i] > peak)
        {
          peak_pos = i;
          peak = coverage[i];
        }
        auto coverageMedianPos = coverage.begin() + peak_pos;
        std::list<unsigned long> coveragePeakEntries;
        //First scroll half a window to the right of the peak (stop if we reach the end)
        for (int i = 0; i < this->windowSize/2 && coverageMedianPos != coverage.end(); ++i) ++coverageMedianPos;
        //Then scroll back 1 full window, adding entries to the list
        for (int i = 0; i < this->windowSize && coverageMedianPos != coverage.begin(); ++i) coveragePeakEntries.push_back(*(coverageMedianPos--));
        coveragePeakEntries.sort();
        double coveragePeakMedian = computeMedian(coveragePeakEntries.size(), coverageMedianPos);
        

        if (coveragePeakMedian >= 100) {
            std::vector<unsigned long> percentileContainer(coverage);
            std::sort(percentileContainer.begin(), percentileContainer.end());
            {
                auto xcursor = percentileContainer.begin();
                while (xcursor != percentileContainer.end() && (*xcursor) == 0ul) ++xcursor;
                percentileContainer.erase(percentileContainer.begin(), xcursor);
            }
            unsigned long lowerLimit = percentileContainer[percentileContainer.size()*0.05];
            unsigned long trimmed_length = 0ul;
            
            {
                auto cursor = coverage.begin();
                while (cursor != coverage.end() && (*cursor) <= lowerLimit)
                {
                    ++trimmed_length;
                    ++cursor;
                }
                coverage.erase(coverage.begin(), cursor);
            }
            {
                while (coverage.size() > 0 && coverage.back() <= lowerLimit) {
                    coverage.pop_back();
                    ++trimmed_length;
                }
            }

            if (coverage.size() >= this->geneLength)
            {
                double windowSize = static_cast<double>(this->windowSize);
                std::vector<double> lcov, rcov;
                lcov.reserve(this->windowSize);
                rcov.reserve(this->windowSize);
                for (unsigned int i = this->offset; i < this->offset + this->windowSize && i < coverage.size(); ++i)
                    lcov.push_back(static_cast<double>(coverage[i]));
                for (int i = coverage.size() - (this->windowSize + this->offset); i >= 0 && i < coverage.size() - this->offset; ++i)
                    rcov.push_back(static_cast<double>(coverage[i]));
                std::sort(lcov.begin(), lcov.end());
                std::sort(rcov.begin(), rcov.end());
                if (gene.strand == Strand::Forward)
                {
                    this->threeEnd[gene.feature_id] += computeMedian(rcov.size(), rcov.begin());
                    this->fiveEnd[gene.feature_id] += computeMedian(lcov.size(), lcov.begin());
                } else
                {
                    this->threeEnd[gene.feature_id] += computeMedian(lcov.size(), lcov.begin());
                    this->fiveEnd[gene.feature_id] += computeMedian(rcov.size(), rcov.begin());
                }
                
            }
            
        }
        

    }
    

    //Extract the bias for a gene
    double BiasCounter::getBias(const std::string &geneID)
    {
        double cov5 = this->fiveEnd[geneID];
        double cov3 = this->threeEnd[geneID];
        if (cov5 + cov3 > 0.0)
        {
            this->countedGenes++;
            return cov3 / (cov5 + cov3);
        }
        return -1.0;
    }
    
    unsigned int BiasCounter::countGenes() const
    {
        return this->countedGenes;
    }


    void add_range(std::vector<unsigned long> &coverage, coord offset, unsigned int length)
    {
        const size_t size = coverage.size();
        for (coord i = offset; i < offset + length && i < size; ++i) coverage[i] += 1ul;
        if (offset + length > size) std::cerr << "Error: Attempted to write more coverage than present on exon. Coverage-based metrics may be inaccurate. This may be a sign of an invalid bam or gtf entry" << std::endl;
    }

    //Compute exon coverage metrics, then stich exons together and compute gene coverage metrics
    std::tuple<double, double, double> computeCoverage(Fasta& fastaReader, std::ofstream &writer, const Feature &gene, const unsigned int mask_size, const std::map<std::string, std::vector<unsigned long> > &coverage, std::map<std::string, ExonCoverage>& totalExonCV, BiasCounter &bias)
    {
        std::vector<std::vector<bool> > coverageMask;
        std::vector<unsigned long> geneCoverage;
        unsigned int maskRemainder = mask_size;
        for (unsigned int i = 0; i < exonsForGene[gene.feature_id].size(); ++i)
        {
            coverageMask.push_back(std::vector<bool>(exonLengths[exonsForGene[gene.feature_id][i]].length, true)); //First store a pre-filled mask for the exon
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

            for (unsigned int j = 0; j < mask.size(); ++j) if (mask[j]) exonSize += 1.0; //count the remaining unmasked length of the exon
            if (exonSize > 0)
            {
                auto maskIter = mask.begin();
                for (auto start = exon_coverage.begin(); start != exon_coverage.end(); ++start)
                    if (*(maskIter++)) exonMean += static_cast<double>(*start) / exonSize;
                maskIter = mask.begin();
                for (auto start = exon_coverage.begin(); start != exon_coverage.end(); ++start)
                    if (*(maskIter++)) exonStd += pow(static_cast<double>(*start) - exonMean, 2.0) / exonSize;
                exonStd = pow(exonStd, 0.5);
                exonStd /= exonMean; //now it's a CV
                
                if (!(std::isnan(exonStd) || std::isinf(exonStd))) {
                    FeatureSpan exonPos = exonLengths[exonsForGene[gene.feature_id][i]];
                    if (fastaReader.hasContig(exonPos.chromosome)) {
                        std::string exonSeq = fastaReader.getSeq(exonPos.chromosome, exonPos.start, exonPos.start + exonPos.length);
                        totalExonCV[exonsForGene[gene.feature_id][i]] = {exonStd, gc(exonSeq)};
                    } else totalExonCV[exonsForGene[gene.feature_id][i]] = {exonStd, -1.0};
                }
            }
            // Reserve and append the exon vector to the growing gene vector
            geneCoverage.reserve(geneCoverage.size() + exon_coverage.size());
            geneCoverage.insert(geneCoverage.end(), exon_coverage.begin(), exon_coverage.end());
        }
        //at this point the gene coverage vector represents an UNMASKED, but complete transcript
        bias.computeBias(gene, geneCoverage); //no masking in bias
        double avg = 0.0, std = 0.0;
        // apply the mask to the full gene vector
        if (mask_size)
        {
            //to account for the mask, erase bases from the vector
            //If the mask is larger than the gene, just erase all of it
            //Otherwise, erase from (end-mask) -> end
            geneCoverage.erase((mask_size > geneCoverage.size() ? geneCoverage.begin() : geneCoverage.end() - mask_size), geneCoverage.end());
            // If there is still coverage area, erase from the front to either the end (if the mask is larger than remaining coverage) or to (front + mask)
            if (geneCoverage.size()) geneCoverage.erase(geneCoverage.begin(), (mask_size > geneCoverage.size() ? geneCoverage.end() : geneCoverage.begin() + mask_size));
        }
        double size = static_cast<double>(geneCoverage.size());
        writer << gene.feature_id << "\t";
        if (size > 0) //If there's still any coverage after applying the mask
        {
            for (auto beg = geneCoverage.begin(); beg != geneCoverage.end(); ++beg)
                avg += static_cast<double>(*beg) / size;
            for (auto base = geneCoverage.begin(); base != geneCoverage.end(); ++base)
                std += std::pow(static_cast<double>(*base) - avg, 2.0) / size;
            std = std::pow(std, 0.5);
            writer << avg << "\t" << std << "\t" << (std / avg) << std::endl;
            return std::make_tuple(avg, std, (std / avg));
        }
        writer << "0\t0\tnan" << std::endl;
        return std::make_tuple(-1, -1, -1);
    }


}

std::ofstream& operator<<(std::ofstream &stream, rnaseqc::Metrics &counter)
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
        "Exonic Reads",
        "Failed Vendor QC",
        "High Quality Reads",
        "Intergenic Reads",
        "Intragenic Reads",
        "Ambiguous Reads",
        "Intronic Reads",
        "Low Mapping Quality",
        "Low Quality Reads",
        "Mapped Duplicate Reads",
        "Mapped Reads",
        "Mapped Unique Reads",
        "Mismatched Bases",
        "Non-Globin Reads",
        "Non-Globin Duplicate Reads",
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
    stream << "Chimeric Fragments\t";
    if (counter.get("Chimeric Fragments_tag"))
    {
        stream << counter.get("Chimeric Fragments_tag") << std::endl;
        stream << "Chimeric Alignment Rate\t" << counter.frac("Chimeric Fragments_tag", "Total Mapped Pairs") << std::endl;
    }
    else
    {
        stream << counter.get("Chimeric Fragments_auto") << std::endl;
        stream << "Chimeric Alignment Rate\t" << counter.frac("Chimeric Fragments_auto", "Total Mapped Pairs") << std::endl;

    }
    for (int i = 0; i < keys.size(); ++i)
        if (keys[i] != "Split Reads" || counter.get("Split Reads"))
            stream << keys[i] << "\t" << counter.get(keys[i]) << std::endl;
    auto beg = counter.counter.begin();
    auto end = counter.counter.end();
    while (beg != end)
    {
        // Manually dump the counters for reads filtered by user supplied tags
        if( beg->first.length() > 17 && beg->first.substr(0,17) == "Filtered by tag: ")
        {
            stream << beg->first << "\t" << beg->second << std::endl;
        }
        ++beg;
    }
    return stream;
}
