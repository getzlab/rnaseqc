//
//  BamReader.hpp
//  RNA-SeQC
//
//  Created by Aaron Graubert on 10/3/18.
//  Copyright © 2018 Aaron Graubert. All rights reserved.
//

#ifndef BamReader_h
#define BamReader_h

#include <stdio.h>
#include <mutex>
#include <string>
#include <SeqLib/BamReader.h>
#include <SeqLib/BamHeader.h>
#include <SeqLib/BamRecord.h>

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

class SeqlibReader : public SynchronizedReader {
    SeqLib::BamReader bam;
public:
    SeqlibReader() : bam() {
    }
    
    bool next(SeqLib::BamRecord&);
    
    const SeqLib::BamHeader getHeader() const {
        return this->bam.Header();
    }
    
    bool open(std::string filepath) {
        this->bam.Open(filepath);
        return this->bam.IsOpen();
    }
    
    void addReference(std::string filepath) {
        this->bam.SetCramReference(filepath);
    }
    
};

typedef SeqLib::BamRecord Alignment;
#endif /* BamReader_h */
