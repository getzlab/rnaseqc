//
//  GTF.cpp
//  IntervalTree
//
//  Created by Aaron Graubert on 6/28/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#include "GTF.h"
#include <exception>
#include <stdexcept>
#include <boost/regex.hpp>

using std::ifstream;
using std::string;
using std::map;

const string EXON_NAME = "exon";
const boost::regex ribosomalPattern("(Mt_)?rRNA"); //For recognizing features which are rRNAs
map<string, string> geneNames, geneSeqs;
map<string, coord> geneLengths, geneCodingLengths, exonLengths;
std::map<std::string, std::vector<std::string>> exonsForGene;
std::vector<std::string> geneList, exonList;
map<string, unsigned int> exon_names;


ifstream& operator>>(ifstream &in, Feature &out)
{
    try{
        string line;
        while(getline(in, line))
        {
            if(line[0] == '#') continue; //not a feature line
            std::istringstream tokenizer(line);
            //get chr#
            string buffer;
            if(!getline(tokenizer, buffer, '\t')) throw gtfException("Unable to parse chromosome. Invalid GTF line: " + line);
            out.chromosome = chromosomeMap(buffer);
            //get track name
            if(!getline(tokenizer, buffer, '\t')) throw gtfException("Unable to parse track. Invalid GTF line: " + line);
            //get feature type
            if(!getline(tokenizer, buffer, '\t')) throw gtfException("Unable to parse feature type. Invalid GTF line: " + line);
            out.type = buffer;
            //get start pos
            if(!getline(tokenizer, buffer, '\t')) throw gtfException("Unable to parse start. Invalid GTF line: " + line);
            out.start = std::stoull(buffer);
            //get stop pos
            if(!getline(tokenizer, buffer, '\t')) throw gtfException("Unable to parse end. Invalid GTF line: " + line);
            out.end = std::stoull(buffer);
            //get score
            if(!getline(tokenizer, buffer, '\t')) throw gtfException("Unable to parse score. Invalid GTF line: " + line);
            //get strand
            if(!getline(tokenizer, buffer, '\t')) throw gtfException("Unable to parse strand. Invalid GTF line: " + line);
            switch(buffer[0])
            {
                case '+':
                    out.strand = Strand::Forward;
                    break;
                case '-':
                    out.strand = Strand::Reverse;
                    break;
                default:
                    out.strand = Strand::Unknown;
            }
            //get frame
            if(!getline(tokenizer, buffer, '\t')) throw gtfException("Unable to parse frame. Invalid GTF line: " + line);
            //get attributes
            if(!getline(tokenizer, buffer, '\t')) throw gtfException("Unable to parse attributes. Invalid GTF line: " + line);
            std::map<string, string> attributes;
            parseAttributes(buffer, attributes);
            if ( out.end < out.start)
                std::cout << "Bad fead feature range:" << out.start << " - " << out.end << std::endl;
            if (out.type == "gene" && attributes.find("gene_id") != attributes.end())
            {
                //Parse gene attributes
                out.feature_id = attributes["gene_id"];
                geneLengths[out.feature_id] = out.end - out.start + 1;
                geneList.push_back(attributes["gene_id"]);
            }
            if (out.type == "transcript" && attributes.find("transcript_id") != attributes.end()) out.feature_id = attributes["transcript_id"];
            if (attributes.find("gene_id") != attributes.end()) out.gene_id = attributes["gene_id"];
            if (out.type == "exon")
            {
                //Parse exon attributes
                if (attributes.find("exon_id") != attributes.end())
                {
                    out.feature_id = attributes["exon_id"];
                }
                else if (attributes.find("gene_id") != attributes.end())
                {
                    out.feature_id = attributes["gene_id"] + "_" + std::to_string(++exon_names[attributes["gene_id"]]);
                    std::cerr << "Unnamed exon: Gene: " << attributes["gene_id"] << " Position: [" << out.start << ", " << out.end <<  "] Inferred Exon Name: " << out.feature_id << std::endl;
                }
                exonList.push_back(out.feature_id);
                exonsForGene[out.gene_id].push_back(out.feature_id);
                geneCodingLengths[out.gene_id] += 1 + (out.end - out.start);
                exonLengths[out.feature_id] = 1 + (out.end - out.start);
            }
            if (attributes.find("transcript_type") != attributes.end()) out.transcript_type = attributes["transcript_type"];
            if (attributes.find("gene_name") != attributes.end()) geneNames[out.feature_id] = attributes["gene_name"];
            else if (attributes.find("gene_id") != attributes.end()) geneNames[out.feature_id] = attributes["gene_id"];
            out.ribosomal = boost::regex_match(out.transcript_type, ribosomalPattern);
            break;
        }

    }
    catch(std::invalid_argument &e)
    {
        throw gtfException(std::string("GTF is in an invalid format: ") + e.what());
    }
    catch(std::exception &e)
    {
        throw gtfException(std::string("Uncountered an unknown error while parsing GTF: ")+e.what());
    }
    return in;
}

std::map<std::string,std::string>& parseAttributes(std::string &intake, std::map<std::string,std::string> &attributes)
{
    std::istringstream tokenizer(intake);
    string buffer;
    while (getline(tokenizer, buffer, ';'))
    {
        std::istringstream splitter(buffer);
        string current;
        getline(splitter, current, '"');
        string key = current.substr(0, current.length()-1);
        while (key[0] == ' ') key = key.substr(1);
        getline(splitter, current, '"');
        attributes[key] = current;
    }
    return attributes;
}

bool operator==(const Feature &a, const Feature &b)
{
    if (a.start != b.start) return false;
    if (a.end != b.end) return false;
    if (a.chromosome != b.chromosome) return false;
    if (a.strand != b.strand) return false;
    if (a.type != b.type) return false;
    if (a.feature_id != b.feature_id) return false;
    return a.transcript_type == b.transcript_type;
}

bool compIntervalStart(const Feature &a, const Feature &b)
{
    return a.start < b.start;
}

bool compIntervalEnd(const Feature &a, const Feature &b)
{
    return a.end < b.end;
}

bool intersectPoint(const Feature &a, const coord x)
{
    return (x >= a.start) && (x <= a.end);
}

bool intersectInterval(const Feature &a, const Feature &b)
{
    return intersectPoint(a, b.start) || intersectPoint(a, b.end) || intersectPoint(b, a.start);
}

int partialIntersect(const Feature &target, const Feature &query)
{
    return intersectInterval(target, query) ? (
                                               1+std::min(target.end, query.end-1) - std::max(target.start, query.start)
                                               ) : 0;
}
