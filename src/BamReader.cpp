//
//  BamReader.cpp
//  RNA-SeQC
//
//  Created by Aaron Graubert on 10/3/18.
//  Copyright © 2018 Aaron Graubert. All rights reserved.
//

#include "BamReader.h"

namespace rnaseqc {
    bool SeqlibReader::next(SeqLib::BamRecord &read)
    {
        // Must uncomment before adding multithreading
        //    std::lock_guard<SeqlibReader> guard(*this);
        try {
            bool ok = this->GetNextRecord(read);
            if (ok) this->read_count++;
            return ok;
        }
        catch (std::runtime_error &e) {
            if (this->user_cram_reference) throw referenceHTSMismatch(std::string("HTSLib was unable to find a suitable reference while decoding a cram: ")+e.what());
            throw;
        }
        return false; // No way to get here
        
    }
}
