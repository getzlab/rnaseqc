//
//  GTF.hpp
//  IntervalTree
//
//  Created by Aaron Graubert on 6/28/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#ifndef GTF_h
#define GTF_h

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>

typedef long long coord;

struct Feature {
    coord start, end;
    unsigned short chromosome;
    short strand;
    std::string type, gene_id, transcript_id, exon_id, transcript_type;
};

//template<> struct std::hash<Feature>
//{
//    size_t operator()(const Feature &obj)
//    {
//        long result = 1, prime = 97;
//        result = result * prime + obj.start;
//        result = result * prime + obj.end;
//        result = result * prime + obj.chromosome;
//        result = result * prime + obj.strand;
//        result = result * prime + std::hash<std::string>()(obj.type);
//        result = result * prime + std::hash<std::string>()(obj.gene_id);
//        return result;
//    }
//};

extern std::map<std::string, unsigned short> chromosomes;

unsigned short chromosomeMap(std::string);
std::ifstream& operator>>(std::ifstream&, Feature&);
std::map<std::string,std::string>& parseAttributes(std::string&, std::map<std::string,std::string>&);

#endif /* GTF_h */
