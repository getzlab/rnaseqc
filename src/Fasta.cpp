//
//  Fasta.cpp
//  IntervalTree
//
//  Created by Aaron Graubert on 5/23/18.
//  Copyright © 2018 Aaron Graubert. All rights reserved.
//

#include "Fasta.h"
#include <boost/filesystem.hpp>
#include <cmath>
#include <utility>

namespace rnaseqc {
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
    
    //Given an internal chromosome ID, get the name it corresponds to
    std::string getChromosomeName(chrom idx)
    {
        for (auto entry = chromosomes.begin(); entry != chromosomes.end(); ++entry) if (entry->second == idx) return entry->first;
        throw invalidContigException("Invalid chromosome index");
    }
    
    //Get reverse complement of a sequence
    void complement(std::string &sequence)
    {
        std::string tmp = sequence;
        auto src = sequence.rbegin();
        for(unsigned int i = 0; src != sequence.rend() && i < tmp.length(); ++src, ++i)
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
    
    //Count GC content in a sequence
    double gc(std::string &sequence)
    {
        double content = 0.0, size = static_cast<double>(sequence.length());
        for (auto base = sequence.begin(); base != sequence.end(); ++base)
            if (*base == 'G' || *base == 'g' || *base == 'C' || *base == 'c') content += 1.0/size;
        return content;
    }
    
    // Open a fasta file
    void Fasta::open(std::string &filename)
    {
        this->_open = true;
        this->reader.open(filename);
        if (!this->reader.is_open())
        {
            throw fileException("Unable to open reference fasta: " +filename);
        }
        std::string index_path = filename + ".fai";
        // Check if the index exists at filepath.fai
        if (boost::filesystem::exists(boost::filesystem::path(filename).replace_extension(".fai")))
            index_path = boost::filesystem::path(filename).replace_extension(".fai").string();
        // otherwise fail if the index doesn't exist at filepath.fasta.fai
        else if (!boost::filesystem::exists(index_path)) throw fileException("Unable to locate fasta index: " + filename);
        // Import chromosome names from index
        std::vector<std::string> contigs = bioio::read_fasta_index_contig_names(index_path);
        for (auto contig = contigs.begin(); contig != contigs.end(); ++contig) chromosomeMap(*contig);
        // Then allow bioio to parse the index
        bioio::FastaIndex tmp_index = bioio::read_fasta_index(index_path);
        for (auto entry = tmp_index.begin(); entry != tmp_index.end(); ++entry) this->contigIndex[chromosomeMap(entry->first)] = entry->second;
        if (!this->contigIndex.size()) throw fileException("No contigs found in fasta index: " + index_path);
    }
    
    //Get a forward strand sequence {contig}:{start}-{end}
    std::string Fasta::getSeq(chrom contig, coord start, coord end)
    {
        return this->getSeq(contig, start, end, Strand::Forward);
    }

    bool Fasta::isOpen() const {
        return this->_open;
    }
    
    //Get a sequence {contig}:{start}-{end}, and optionally return its reverse complement
    std::string Fasta::getSeq(chrom contig, coord start, coord end, Strand strand)
    {
        //NOTE: Coordinates must be 0-based, end-exclusive.
        if (!this->isOpen()) return "";
        std::string output;
        // Determine the coordinate for the start of the page which contains the start of this sequence
        coord pageOffset = (floor(start / PAGE_SIZE) * PAGE_SIZE);
        for (coord i = pageOffset; i < end; i+=PAGE_SIZE) //Iterate over pages until we have all the pages required
        {
            this->calls++; // Increment number of pages that were requested
            indexType page = this->pageForCoord(contig, i); // Get page index corresponding to this coordinate
            if (!this->pageCache.count(page)) // If that page isn't cached
            {
                this->misses++; // Increment number of pages that were actually read
                this->pageCache[page] = this->readSeq(contig, i); // Read page from fasta
            }
            this->updateLRU(page); // Update the cache state
            output += this->pageCache[page]; // Append this page to the output
        }
        if (start-pageOffset >= output.size())
        {
            std::cerr << "Unable to fetch sequence" << std::endl;
            std::cerr << "Target region (GTF+1):\t" << getChromosomeName(contig) << ":" << start+1 << "-" << end << std::endl;
            std::cerr << "# pages fetched:\t" << output.length() / PAGE_SIZE << std::endl;
            std::cerr << "This contig page indices:\t[" << this->pageForContig(contig) << ", " << this->pageForContig(contig+1) << ")" << std::endl;
            std::cerr << "Sequence page indices:\t[" << this->pageForCoord(contig, start) << ", " << this->pageForCoord(contig, end) << "]" << std::endl;
        }
        // Extract desired sequence from the output (which is complete pages)
        output = output.substr(start-pageOffset, end-start);
        if (strand == Strand::Reverse) complement(output);
        return output;
    }
    
    // Bump page to top of the LRU, dropping older pages as necessary
    void Fasta::updateLRU(indexType page)
    {
        this->lru.remove(page);
        while (this->lru.size() >= CACHE_SIZE)
        {
            this->pageCache.erase(this->lru.front());
            this->lru.pop_front();
        }
        this->lru.push_back(page);
    }
    
    // Get page idx for position 0 of a contig
    indexType Fasta::pageForContig(chrom contig)
    {
        static std::unordered_map<chrom, indexType> pageIndex; // Function caches the lookup table to save time
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
    
    // Read a full page starting at this position
    std::string Fasta::readSeq(chrom contig, coord pos)
    {
        if (!this->contigIndex.count(contig)) throw invalidContigException("No such contig: " + getChromosomeName(contig));
        return (bioio::read_fasta_contig(this->reader, this->contigIndex[contig], pos, PAGE_SIZE));
    }
    
    indexType Fasta::pageForCoord(chrom contig, coord pos)
    {
        return this->pageForContig(contig) + floor(static_cast<double>(pos)/PAGE_SIZE);
    }
    
    Fasta::~Fasta()
    {
        this->reader.close();
        this->pageCache.clear();
        if (this->misses) std::cerr << this->misses << " cache misses out of " << this->calls << " requests" << std::endl;
    }

    bool Fasta::hasContig(chrom contig) const {
        return this->contigIndex.count(contig);
    }
}

