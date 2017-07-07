//
//  RefSeqReader.hpp
//  IntervalTree
//
//  Created by Aaron Graubert on 7/6/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#ifndef RefSeqReader_h
#define RefSeqReader_h

#include <fstream>
#include <string>
#include <map>

class RefSeqReader
{
    std::ifstream reader;
    std::map<std::string, unsigned long long> index; //stores feature name -> byte offset
public:
    RefSeqReader(std::string filename);
    RefSeqReader(std::string filename, std::string indexfile);
    ~RefSeqReader()
    {
        this->reader.close();
    }
    void indexFeature(std::string, std::string, unsigned long); // new feature name, host feature name, start pos
    std::string read(std::string, long); //reference feature name, string length (negative for a antisense read)
};

char invert(char);
std::string invert(std::string);

#endif /* RefSeqReader_h */
