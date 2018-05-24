//
//  RefSeq.cpp
//  IntervalTree
//
//  Created by Aaron Graubert on 5/23/18.
//  Copyright © 2018 Aaron Graubert. All rights reserved.
//

#include "Fasta.h"
#include <boost/filesystem.hpp>
#include <cmath>
#include <utility>

std::map<std::string, chrom> chromosomes;

chrom chromosomeMap(std::string chr)
{
    auto entry = chromosomes.find(chr);
    if (entry == chromosomes.end())
    {
        chromosomes[chr] = chromosomes.size() + 1u;
    }
    return chromosomes[chr];
}

std::string getChromosomeName(chrom idx)
{
    for (auto entry = chromosomes.begin(); entry != chromosomes.end(); ++entry) if (entry->second == idx) return entry->first;
    throw invalidContigException("Invalid chromosome index");
}

void complement(std::string &sequence)
{
    std::string tmp = sequence;
    auto src = sequence.rbegin();
    unsigned int i = 0;
    for(; src != sequence.rend() && i < tmp.length(); ++src, ++i)
    {
        switch(*src)
        {
            case 'A':
            case 'a':
                tmp[i] = 'T';
                break;
            case 'T':
            case 't':
                tmp[i] = 'A';
                break;
            case 'C':
            case 'c':
                tmp[i] = 'G';
                break;
            case 'G':
            case 'g':
                tmp[i] = 'C';
                break;
            default:
                tmp[i] = *src;
        }
    }
    sequence.swap(tmp);
}

double gc(std::string &sequence)
{
    double content = 0.0, size = static_cast<double>(sequence.length());
    for (auto base = sequence.begin(); base != sequence.end(); ++base)
        if (*base == 'G' || *base == 'g' || *base == 'C' || *base == 'c') content += 1.0/size;
    return content;
}

void Fasta::open(std::string &filename)
{
    this->isOpen = true;
    this->reader.open(filename);
    if (!this->reader.is_open())
    {
        throw fileException("Unable to open reference fasta: " +filename);
    }
    std::string index_path = filename + ".fai";
    if (boost::filesystem::exists(index_path))
    {
        std::vector<std::string> contigs = bioio::read_fasta_index_contig_names(index_path);
        for (auto contig = contigs.begin(); contig != contigs.end(); ++contig) chromosomeMap(*contig);
        bioio::FastaIndex tmp_index = bioio::read_fasta_index(index_path);
        for (auto entry = tmp_index.begin(); entry != tmp_index.end(); ++entry) this->contigIndex[chromosomeMap(entry->first)] = entry->second;
    }
    else if (boost::filesystem::exists(boost::filesystem::path(filename).replace_extension(".fai")))
    {
        index_path = boost::filesystem::path(filename).replace_extension(".fai").string();
        std::vector<std::string> contigs = bioio::read_fasta_index_contig_names(index_path);
        for (auto contig = contigs.begin(); contig != contigs.end(); ++contig) chromosomeMap(*contig);
        bioio::FastaIndex tmp_index = bioio::read_fasta_index(index_path);
        for (auto entry = tmp_index.begin(); entry != tmp_index.end(); ++entry) this->contigIndex[chromosomeMap(entry->first)] = entry->second;
    }
    else throw fileException("Unable to locate fasta index: " + filename);
    if (!this->contigIndex.size()) throw fileException("No contigs found in fasta index: " + index_path);
//    this->lru.reserve(CACHE_SIZE);
}

std::string Fasta::getSeq(chrom contig, coord start, coord end)
{
    return this->getSeq(contig, start, end, false);
}

std::string Fasta::getSeq(chrom contig, coord start, coord end, bool revComp)
{
    if (!this->isOpen) return "";
    this->calls++;
    std::string output;
    for (coord i = start; i <= end; i+=PAGE_SIZE)
    {
        indexType page = this->pageForCoord(contig, i);
        if (!this->pageCache.count(page))
        {
            this->misses++;
            this->pageCache[page] = this->readSeq(contig, i);
        }
        this->updateLRU(page);
        output += this->pageCache[page];
    }
    if (!revComp) return (output.substr(0, end-start+1));
    else
    {
        output = output.substr(0, end-start+1);
        complement(output);
        return output;
    }
}

void Fasta::updateLRU(indexType page)
{
    /**/
    this->lru.remove(page);
    /*/bool isIn = false;
    auto pos = this->lru.begin();
    for (; pos != this->lru.end(); ++pos)
        if (*pos == page)
        {
            isIn = true;
            break;
        }
    if (isIn) {
        this->lru.erase(pos);
    }
     else /**/while (this->lru.size() >= CACHE_SIZE)
     {
         this->pageCache.erase(this->lru.front());
         this->lru.pop_front();
     }
    this->lru.push_back(page);
}

indexType Fasta::pageForContig(chrom contig)
{
    static std::unordered_map<chrom, indexType> pageIndex;
    if (!pageIndex.size()) pageIndex[0] = 0;
    if (pageIndex.count(contig)) return pageIndex[contig];
    if (!this->contigIndex.count(contig)) throw invalidContigException("No such contig: " + getChromosomeName(contig));
    chrom firstContig = contig;
    for (; firstContig > 0; --firstContig) if (pageIndex.count(firstContig)) break;
    indexType idx = pageIndex[firstContig];
    for (chrom i = firstContig; i < contig; ++i)
    {
        idx += ceil(static_cast<double>(this->contigIndex[i].length)/PAGE_SIZE);
        pageIndex[i+1] = idx;
        
    }
    return idx;
}

std::string Fasta::readSeq(chrom contig, coord pos)
{
    if (!this->contigIndex.count(contig)) throw invalidContigException("No such contig: " + getChromosomeName(contig));
    return (bioio::read_fasta_contig(this->reader, this->contigIndex[contig], pos - 1, PAGE_SIZE));
}

indexType Fasta::pageForCoord(chrom contig, coord pos)
{
    return this->pageForContig(contig) + floor(static_cast<double>(pos - 1)/PAGE_SIZE);
}

Fasta::~Fasta()
{
    this->reader.close();
    this->pageCache.clear();
    std::cout << this->misses << "cache misses out of" << this->calls << "requests" << std::endl;
}

