//
//  BamReader.hpp
//  RNA-SeQC
//
//  Created by Aaron Graubert on 10/3/18.
//  Copyright © 2018 Aaron Graubert. All rights reserved.
//

#ifndef BamReader_h
#define BamReader_h

#include "Fasta.h"
#include <stdio.h>
#include <mutex>
#include <string>
#include <set>
#include <SeqLib/BamReader.h>
#include <SeqLib/BamHeader.h>
#include <SeqLib/BamRecord.h>
#include <htslib/cram/cram.h> // I really don't like using unofficial APIs, but not much choice here.

namespace rnaseqc {

    struct referenceHTSMismatch : public std::exception {
        std::string error;
        referenceHTSMismatch(std::string msg) : error(msg) {};
    };


    class SynchronizedReader {
        std::mutex mtx;
    protected:
        unsigned long read_count;
    public:
        SynchronizedReader() : mtx(), read_count() {
            
        }
        
        void lock()
        {
            this->mtx.lock();
        }
        
        void unlock()
        {
            this->mtx.unlock();
        }
        
        unsigned long get_count() const
        {
            return this->read_count;
        }
    };
    
    class SeqlibReader : public SynchronizedReader, public SeqLib::BamReader {
        std::string reference_path;
        std::set<chrom> valid_chroms;
        bool user_cram_reference;
    public:
        
        SeqlibReader() : reference_path(), valid_chroms(), user_cram_reference(false) {}
        
        bool next(SeqLib::BamRecord&);
        
        const SeqLib::BamHeader getHeader() const {
            return this->Header();
        }
        
        bool open(std::string filepath) {
            if (this->reference_path.length()) {
                auto htsfile = hts_open(filepath.c_str(), "r");
                hts_set_fai_filename(htsfile, this->reference_path.c_str());
                if (htsfile->format.format == htsExactFormat::cram) {
                    this->user_cram_reference = true;
                    // Cram handling is very dumb. All of this nonsense is just because htslib is incredibly opaque about reference handling
                    // Even with a user-provided reference, htslib only uses it if the MD5 matches
                    // So here we load up the file, check if it's a cram, then get a list of chromosomes that htslib decides to use
                    cram_fd *cram = static_cast<cram_fd*>(htsfile->fp.cram);
                    if (cram->refs && cram->refs->nref > 0)
                        for (unsigned int i = 0; i < cram->refs->nref; ++i)
                            if (this->reference_path == std::string(cram->refs->ref_id[i]->fn))
                                this->valid_chroms.insert(
                                    chromosomeMap(cram->refs->ref_id[i]->name)
                                );
                    this->SetCramReference(this->reference_path); // Consider moving out of if statement, if there's any meaningful use to having a reference set on a non-cram
                }
                hts_close(htsfile);
            }
            this->Open(filepath);
            return this->IsOpen();
        }
        
        void addReference(std::string filepath) {
            this->reference_path = filepath;
        }
        
        inline bool validateChromosome(const chrom c) {
            // For crams, we only validate chromosomes which matched our reference. Otherwise yes!
            return this->user_cram_reference ? this->valid_chroms.count(c) > 0 : true;
        }
        
    };
    
    typedef SeqLib::BamRecord Alignment;
}

#endif /* BamReader_h */
