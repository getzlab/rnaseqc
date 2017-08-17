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
#include <vector>
#include <sstream>

typedef long long coord;

struct Feature {
    coord start, end;
    unsigned short chromosome;
    short strand;
    std::string type, feature_id, gene_id, transcript_type, transcript_id;
    bool ribosomal;
};

bool operator==(const Feature &a, const Feature &b);
bool compIntervalStart(const Feature&, const Feature&);
bool compIntervalEnd(const Feature&, const Feature&);
bool intersectPoint(const Feature&, const coord);
bool intersectInterval(const Feature&, const Feature&);
int partialIntersect(const Feature&, const Feature&);

extern std::map<std::string, unsigned short> chromosomes;
extern std::map<std::string, std::string> geneNames;
extern std::map<std::string, coord> geneLengths;
extern std::map<std::string, coord> transcriptCodingLengths;
extern std::vector<std::string> geneList, exonList;

unsigned short chromosomeMap(std::string);
std::ifstream& operator>>(std::ifstream&, Feature&);
std::map<std::string,std::string>& parseAttributes(std::string&, std::map<std::string,std::string>&);

#endif /* GTF_h */
