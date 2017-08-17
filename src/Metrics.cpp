//
//  Metrics.cpp
//  IntervalTree
//
//  Created by Aaron Graubert on 7/5/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#include "Metrics.h"
#include <iostream>


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
    return (double) this->get(a) / this->get(b);
}

std::ofstream& operator<<(std::ofstream &stream, Metrics &counter)
{
    unsigned int NUM_KEYS = 33;
    std::string keys[] = {
        "Alternative Alignments",
        "Chimeric Pairs",
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
    for (int i = 0; i < NUM_KEYS; ++i) stream << keys[i] << "\t" << counter.get(keys[i]) << std::endl;
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
    return (bool) this->data[gene_id].size();
}

bool Collector::isDirty()
{
    return this->dirty;
}

double Collector::sum()
{
    return this->total;
}
