//
//  RefSeqReader.cpp
//  IntervalTree
//
//  Created by Aaron Graubert on 7/6/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#include "RefSeqReader.h"
#include <sstream>
using std::string;

RefSeqReader::RefSeqReader(string filename) : reader(filename), index()
{
    string line;
    while(getline(this->reader, line))
    {
        if (line.size() >= 2 && line[0] == '>')
        {
            //this is a header line
            std::istringstream parser(line);
            string header;
            getline(parser, header, ' '); //read the header name
            this->index[header] = this->reader.tellg();
        }
    }
}

RefSeqReader::RefSeqReader(std::string filename, std::string indexfile) : reader(filename), index()
{
    string line;
    std::ifstream indexReader(indexfile);
    while(getline(indexReader, line))
    {
        std::istringstream parser(line);
        string header, offset;
        getline(parser, header, '\t'); //read header name
        getline(parser, offset, '\t'); //skip over length
        getline(parser, offset, '\t'); //read byte offset
        this->index[header] = std::stoull(offset);
    }
}

void RefSeqReader::indexFeature(string featureName, string referenceName, unsigned long start)
{
    this->reader.seekg(this->index.at(referenceName));
    string line;
    unsigned long long pos;
    while (start > 0)
    {
        if (start > line.size())
        {
            start -= line.size();
            pos = reader.tellg();
            getline(this->reader, line);
        }
        else
        {
            this->index[featureName] = pos + start;
            return;
        }
    }
}

string RefSeqReader::read(string featureName, long length)
{
    string output;
    this->reader.seekg(this->index.at(featureName));
    if (length > 0)
    {
        char* buffer = new char[length];
        this->reader.read(buffer, length);
        output = buffer;
        delete[] buffer;
        buffer = nullptr;
    }
    else
    {
        
    }
    return output;
}

char invert(char input);

std::string invert(std::string);

