//
//  BED.cpp
//  IntervalTree
//
//  Created by Aaron Graubert on 7/11/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#include "BED.h"
#include <sstream>
#include <exception>
#include <stdexcept>

using std::ifstream;
using std::string;

namespace rnaseqc {
    ifstream& extractBED(ifstream &input, Feature &out)
    {
        try
        {
            string line;
            while(getline(input, line))
            {
                if(line[0] == '#') continue; //Do beds even have comment lines?
                std::istringstream tokenizer(line);
                string buffer;
                tokenizer >> buffer; //chromosome name
                out.chromosome = chromosomeMap(buffer);
                tokenizer >> buffer; //start
                out.start = std::stoull(buffer) + 1;
                tokenizer >> buffer; //stop
                out.end = std::stoull(buffer) + 1;
                out.feature_id = line; // add a dummy exon_id for mapping interval intersections later
                out.type = FeatureType::Exon;
                break;
            }
        }
        catch (std::exception &e)
        {
            throw bedException(std::string("Encountered an unknown error while parsing the BED: ") + e.what());
        }
        return input;
    }

}
