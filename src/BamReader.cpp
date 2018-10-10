//
//  BamReader.cpp
//  RNA-SeQC
//
//  Created by Aaron Graubert on 10/3/18.
//  Copyright © 2018 Aaron Graubert. All rights reserved.
//

#include "BamReader.h"

bool SeqlibReader::next(SeqLib::BamRecord &read)
{
    // Must uncomment before adding multithreading
//    std::lock_guard<SeqlibReader> guard(*this);
    bool ok = this->bam.GetNextRecord(read);
    if (ok) this->read_count++;
    return ok;
}
