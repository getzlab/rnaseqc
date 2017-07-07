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

typedef long long coord;

struct Feature {
    coord start, end;
    unsigned short chromosome;
    short strand;
    std::string type, gene_id, transcript_id, transcript_type;
};

//extern std::map<unsigned short, int> chromosomeMap;

unsigned short chromosomeMap(std::string);
std::ifstream& operator>>(std::ifstream&, Feature&);
std::map<std::string,std::string>& parseAttributes(std::string&, std::map<std::string,std::string>&);

#endif /* GTF_h */
