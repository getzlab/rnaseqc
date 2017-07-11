// IntervalTree.cpp : Defines the entry point for the console application.
//
#include "Interval.h"
#include "BED.h"
#include "Metrics.h"
#include <string>
#include <iostream>
#include <exception>
#include <stdio.h>
#include <set>
#include <regex>
#include <list>
#include <ctime>
#include <math.h>
#include <bamtools/api/BamReader.h>
#include <bamtools/api/BamAlignment.h>
//TODO: also intake a reference fasta for checking reference sequences
    //Generate refSeq Index if none provided
    //Use index to skip to chromosome identifiers as necessary
using namespace std;
using namespace BamTools;

const int CHIMERIC_DISTANCE = 2000000;

int usage(int, string,string);

vector<Feature>& extractBlocks(BamAlignment&, vector<Feature>&, unsigned short);

int main(int argc, char* argv[])
{
	try
	{
		if (argc < 3 || argc > 4) return usage(1, argv[0], "Must provide exactly two or three arguments");
        unsigned int doFragmentSize = 0u;
        ifstream reader(argv[1]);
		cout<<"Reading GTF Features..."<<endl;
        time_t t0, t1, t2;
        time(&t0);
        clock_t start_clock = clock();
        map<unsigned short, vector<Feature> > features;
        map<unsigned short, vector<Feature> > *bedFeatures; //for fragment calculations
        Feature line;
        regex rRNA("rRNA");
        int readLength = 0;
		while ((reader >> line))
		{
            if (line.type == "gene" || line.type == "exon")
            {
                features[line.chromosome].push_back(line);
            }
		}
        time(&t1);
        cout << "Finished processing GTF in " << difftime(t1, t0) << " seconds" << endl;
        map<unsigned short, IntervalTree*> trees;
        map<unsigned short, IntervalTree*> *bedTrees;
        for(auto bin = features.begin(); bin != features.end(); ++bin)
        {
            cout << "Building tree for chromosome " << bin->first << " with " << bin->second.size() << " features" << endl;
            trees[bin->first] = constructIntervalTree(bin->second);
        }
        if (argc == 4)
        {
            cout << "Parsing BED intervals for fragment size computations..." << endl;
            doFragmentSize = 1000000u;
            bedFeatures = new map<unsigned short, vector<Feature> >();
            bedTrees = new map<unsigned short, IntervalTree*>();
            ifstream bedReader(argv[3]);
            while (extractBED(bedReader, line)) (*bedFeatures)[line.chromosome].push_back(line);
            bedReader.close();
            for (auto bin = bedFeatures->begin(); bin != bedFeatures->end(); ++bin)
                (*bedTrees)[bin->first] = constructIntervalTree(bin->second);
        }
        BamReader bam;
		bam.Open(argv[2]);
        bam.LocateIndex(); //necessary?
		BamAlignment alignment;
        //const unsigned short CHR_CUTOFF = chromosomeMap("~~Have a wonderful day!~~");
		cout<<"Parsing bam..."<<endl;
        SamSequenceDictionary sequences = bam.GetHeader().Sequences;
        Metrics counter;
        map<string, string> fragments;
        //map<string, string> fragmentExons; //this is dumb... We can only compute fragment size using pairs on the same exon
        list<long long> fragmentSizes;
        map<string, double> geneCoverage, transcriptCoverage, exonCoverage;
        unsigned long long alignmentCount = 0ull;
        time_t report_time;
        time(&report_time);
        time(&t2);
		while (bam.GetNextAlignment(alignment))
		{
            ++alignmentCount;
            if (alignmentCount % 250000 == 0) time(&t2);
            if (difftime(t2, report_time) >= 10)
            {
                time(&report_time);
                cout << "Time elapsed: " << difftime(t2, t1) << "; Alignments processed: " << alignmentCount << endl;
            }
            if (!alignment.IsPrimaryAlignment()) counter.increment("Alternative Alignments");
            else if (alignment.IsFailedQC()) counter.increment("Failed Vendor QC");
            else if (alignment.MapQuality < 255u) counter.increment("Low quality reads");
            if (alignment.IsPrimaryAlignment() && !alignment.IsFailedQC() && alignment.MapQuality == 255u)
            {
                counter.increment("Total Reads");
                if (!alignment.IsPaired()) counter.increment("Unpaired Reads");
                if (alignment.IsDuplicate()) counter.increment("Duplicate Reads");
                if (alignment.Length > readLength) readLength = alignment.Length;
                if (alignment.IsMapped())
                {
                    counter.increment("Mapped Reads");
                    if (alignment.IsPaired())
                    {
                        unsigned int mismatches;
                        alignment.GetTag("NM", mismatches);
                        counter.increment("Mismatches", mismatches);
                        if (alignment.IsFirstMate())
                        {
                            counter.increment("End 1 Mapped Reads");
                            counter.increment("End 1 Mismatches", mismatches);
                            counter.increment("End 1 Bases", alignment.Length);
                        }
                        else
                        {
                            counter.increment("End 2 Mapped Reads");
                            counter.increment("End 2 Mismatches", mismatches);
                            counter.increment("End 2 Bases", alignment.Length);
                        }
                        counter.increment("Total Bases", alignment.Length);
                    }
                    if (alignment.IsDuplicate()) counter.increment("Mapped Duplicate Reads");
                    else counter.increment("Mapped Unique Reads");
                    //check length against max read length
                    if (alignment.IsPaired() && alignment.IsMateMapped())
                    {
                        if (alignment.IsFirstMate()) counter.increment("Total Mapped Pairs");
                        if (alignment.RefID != alignment.MateRefID || abs(alignment.Position - alignment.MatePosition) > CHIMERIC_DISTANCE)
                        {
                            counter.increment("Chimeric Pairs");
                            cout << (sequences.Begin()+alignment.RefID)->Name << "\t" << (sequences.Begin()+alignment.MateRefID)->Name << ":\t" << (alignment.RefID == alignment.MateRefID) << endl;
                        }
                    }
                    if (alignment.RefID == -1) counter.increment("Intergenic Reads");
                    //else check for rods? count intergenic if no rods, otherwise count overlap if more than one rod
                }
                
                if(alignment.RefID >= 0 && alignment.RefID < sequences.Size())
                {
                    //map the 0...N RefID to an internal chromosome number consistent with the GTF input
                    string chrName = (sequences.Begin()+alignment.RefID)->Name;
                    unsigned short chr = chromosomeMap(chrName);
                    vector<Feature> blocks;
                    extractBlocks(alignment, blocks, chr);
                    bool split = false;
                    long long lastEnd = -1;
                    for(auto block = blocks.begin(); block != blocks.end(); ++block)
                    {
                        if (lastEnd > 0 && !split) split = (block->start - lastEnd) > 100;
                        lastEnd = block->end;
                    }
                    if (split) counter.increment("Split Reads");
                    
                    if (trees.find(chr) != trees.end())
                    {
                        counter.increment("Alignment Blocks", blocks.size());
                        Interval current;
                        current.start = alignment.Position + 1; //0-based + 1 == 1-based
                        current.end = alignment.GetEndPosition(); //0-based, open == 1-based, closed
                        
                        bool exonic = false, intragenic = false, intronic = false, transcriptPlus = false, transcriptMinus = false, ribosomal = false;
                        string exonName;
                        
                        for (auto block = blocks.begin(); block != blocks.end(); ++block)
                        {
                            resultSet *results = trees[chr]->queryInterval(*block);
                            auto beg = results->begin();
                            auto end = results->end();
                            while (beg != end)
                            {
                                intragenic = true;
                                if (beg->strand == 1) transcriptPlus = true;
                                else if (beg->strand == -1) transcriptMinus = true;
                                if (beg->type == "exon")
                                {
                                    exonic = true;
                                    transcriptCoverage[beg->transcript_id] += ((double) block->end - block->start)/alignment.Length;
                                    double tmp = ((double) block->end - block->start) / alignment.Length;
                                    exonCoverage[beg->exon_id] += tmp;
                                }
                                else if (beg->type == "gene")
                                {
                                    intragenic = true;
                                    geneCoverage[beg->gene_id] += ((double) block->end - block->start)/alignment.Length;
                                }
//                                if (regex_search(beg->transcript_type, rRNA)) ribosomal = true;
                                if (beg->transcript_type == "rRNA") ribosomal = true;
                                ++beg;
                            }
                            delete results;
                        }
                        
                        
//                        vector<Interval> *results = tree->queryInterval(current);
//                        auto beg = results->begin();
//                        auto end = results->end();
//                        while (beg != end)
//                        {
//                            if (beg->chromosome == chr)
//                            {
//                                intragenic = true;
//                                if (beg->strand == 1) transcriptPlus = true;
//                                else if (beg->strand == -1) transcriptMinus = true;
//                                if (beg->type == "exon")
//                                {
//                                    exonic = true;
//                                    genes.insert(beg->gene_id);
//                                    transcripts.insert(beg->transcript_id);
//                                }
//                                else if (beg->type == "gene") intragenic = true;
//                                if (regex_search(beg->transcript_type, rRNA)) ribosomal = true;
//                            }
//                            else counter.increment("Chromosomal mismatches");
//                            ++beg;
//                        }
                        intronic = (intragenic && !exonic);
                        if (exonic)
                        {
                            counter.increment("Exonic");
//                            if (alignment.IsPaired() && okForFragmentSize)
//                            {
//                                auto fragment = fragments.find(alignment.Name);
//                                if (fragment == fragments.end()) //first time we've encountered a read in this pair
//                                {
//                                    fragments[alignment.Name] = exonName;
//                                }
//                                else if (exonName == fragments[alignment.Name])
//                                {
//                                    //This pair is useable for fragment statistics:  both pairs fully aligned to the same exon
//                                    if (alignment.IsFirstMate()) fragmentSizes.push_back(alignment.MatePosition - alignment.Position + alignment.Length);
//                                    else fragmentSizes.push_back(alignment.GetEndPosition() - alignment.MatePosition);
//                                    fragments.erase(fragment);
//                                }
//                            }
                        }
                        if (intragenic) counter.increment("Intragenic Reads");
                        if (!(exonic || intronic)) counter.increment("Intergenic Reads");
                        if (intronic) counter.increment("Intron/UTR");
                        if (ribosomal) counter.increment("rRNA");
                        
                        if (transcriptMinus ^ transcriptPlus)
                        {
                            if (alignment.IsPaired())
                            {
                                if (alignment.IsFirstMate())
                                {
                                    if (alignment.IsReverseStrand()) transcriptMinus ? counter.increment("End 1 Sense") : counter.increment("End 1 Antisense");
                                    else transcriptPlus ? counter.increment("End 1 Sense") : counter.increment("End 1 Antisense");
                                }
                                else
                                {
                                    if (alignment.IsReverseStrand()) transcriptMinus ? counter.increment("End 2 Sense") : counter.increment("End 2 Antisense");
                                    else transcriptPlus ? counter.increment("End 2 Sense") : counter.increment("End 2 Antisense");
                                }
                            }
                        }
                    }
                   if (doFragmentSize && bedTrees != nullptr && bedTrees->find(chr) != bedTrees->end())
                   {
                       for (auto block = blocks.begin(); block != blocks.end(); ++block)
                       {
                           resultSet *results = (*bedTrees)[chr]->queryInterval(*block);
                           if (results->size() == 1)
                           {
                               if (alignment.IsPaired())
                               {
                                   auto fragment = fragments.find(alignment.Name);
                                   string exonName = results->begin()->exon_id;
                                   if (fragment == fragments.end()) //first time we've encountered a read in this pair
                                   {
                                       fragments[alignment.Name] = exonName;
                                   }
                                   else if (exonName == fragments[alignment.Name])
                                   {
                                       //This pair is useable for fragment statistics:  both pairs fully aligned to the same exon
                                       if (alignment.IsFirstMate()) fragmentSizes.push_back(alignment.MatePosition - alignment.Position + alignment.Length);
                                       else fragmentSizes.push_back(alignment.GetEndPosition() - alignment.MatePosition);
                                       fragments.erase(fragment);
                                       --doFragmentSize;
                                       if (!doFragmentSize)
                                       {
                                           for (auto tree = bedTrees->begin(); tree != bedTrees->end(); ++tree) delete tree->second;
                                           delete bedTrees;
                                           delete bedFeatures;
                                       }
                                   }
                               }
                           }
                           delete results;
                       }
                   }
                }
                else counter.increment("Reads discarded by sequence identifier");
            }
            
		}
        for (auto tree = trees.begin(); tree != trees.end(); ++tree) delete tree->second;
        time(&t2);
        cout<< "Time Elapsed: " << difftime(t2, t1) << "; Alignments processed: " << alignmentCount << endl;
        cout << "Total runtime: " << difftime(t2, t0) << "; Total CPU Time: " << (clock() - start_clock)/CLOCKS_PER_SEC << endl;
        
        cout << "Generating report" << endl;
        
        ofstream checksum("checksum.txt");
        for(auto fragment = fragmentSizes.begin(); fragment != fragmentSizes.end(); ++fragment)
        {
            checksum << abs(*fragment) << endl;
        }
        checksum.close();
        
        double fragmentAvg = 0.0, fragmentStd = 0.0;
        for(auto fragment = fragmentSizes.begin(); fragment != fragmentSizes.end(); ++fragment)
        {
            fragmentAvg += (double) abs(*fragment) / fragmentSizes.size();
        }
        //fragmentAvg /= (double) fragmentSizes.size();
        for(auto fragment = fragmentSizes.begin(); fragment != fragmentSizes.end(); ++fragment)
        {
            fragmentStd += pow((double) (*fragment) - fragmentAvg, 2.0) / (double) fragmentSizes.size();
        }
        fragmentStd = pow(fragmentStd, 0.5);
        
        unsigned int genesDetected = 0, transcriptsDetected = 0 ;
        
        ofstream geneReport("geneReport.tsv");
        for(auto gene = geneCoverage.begin(); gene != geneCoverage.end(); ++gene)
        {
            geneReport << gene->first << "\t" << gene->second << endl;
            if (gene->second >= 5.0) ++genesDetected;
        }
        geneReport.close();
        
        ofstream transcriptReport("transcriptReport.tsv");
        for(auto transcript = transcriptCoverage.begin(); transcript != transcriptCoverage.end(); ++transcript)
        {
            transcriptReport << transcript->first << "\t" << transcript->second << endl;
            if (transcript->second >= 5.0) ++transcriptsDetected;
        }
        transcriptReport.close();
        
        ofstream exonReport("exonReport.tsv");
        for(auto exon = exonCoverage.begin(); exon != exonCoverage.end(); ++exon)
        {
            exonReport << exon->first << "\t" << exon->second << endl;
        }
        exonReport.close();
        
        ofstream output("report.tsv");
        output << "Sample\t" << argv[2] << endl;
        output << "Mapping Rate\t" << counter.frac("Mapped Reads", "Total Reads") << endl;
        output << "Unique Rate of Mapped\t" << counter.frac("Mapped Unique Reads", "Mapped Reads") << endl;
        output << "Duplicate Rate of Mapped\t" << counter.frac("Mapped Duplicate Reads", "Mapped Reads") << endl;
        output << "Base Mismatch\t" << counter.frac("Mismatches", "Total Bases") << endl;
        output << "End 1 Mapping Rate\t"<< 2.0 * counter.frac("End 1 Mapped Reads", "Total Reads") << endl;
        output << "End 2 Mapping Rate\t"<< 2.0 * counter.frac("End 2 Mapped Reads", "Total Reads") << endl;
        output << "End 1 Mismatch Rate\t" << counter.frac("End 1 Mismatches", "End 1 Bases") << endl;
        output << "End 2 Mismatch Rate\t" << counter.frac("End 2 Mismatches", "End 2 Bases") << endl;
        output << "Expression Profiling Efficiency\t" << counter.frac("Exonic", "Total Reads") << endl;
        output << "Exonic Rate\t" << counter.frac("Exonic", "Mapped Reads") << endl;
        output << "Intronic Rate\t" << counter.frac("Intron/UTR", "Mapped Reads") << endl;
        output << "Intergenic Rate\t" << counter.frac("Intergenic Reads", "Mapped Reads") << endl;
        output << "Intragenic Rate\t" << counter.frac("Intragenic Reads", "Mapped Reads") << endl;
        output << "rRNA Rate\t" << counter.frac("rRNA", "Mapped Reads") << endl;
        output << "End 1 Sense Rate\t" << (double) counter.get("End 1 Sense") / (counter.get("End 1 Sense") + counter.get("End 1 Antisense")) << endl;
        output << "End 2 Sense Rate\t" << (double) counter.get("End 2 Sense") / (counter.get("End 2 Sense") + counter.get("End 2 Antisense")) << endl;
        output << "Avg. Blocks per Read\t" << counter.frac("Alignment Blocks", "Mapped Reads") << endl;
        
        output << counter;
        output << "Read Length\t" << readLength << endl;
        output << "Genes Detected\t" << genesDetected << endl;
        output << "Transcripts Detected\t" << transcriptsDetected << endl;
        output << "Average Fragment Length\t" << fragmentAvg << endl;
        output << "Fragment Length Std\t" << fragmentStd << endl;
        output.close();
        
//        Note  Adrenal
//        Chimeric Pairs  0 <----
//        Fragment Length StdDev  273
//        Estimated Library Size  0
//        Fragment Length Mean  259


	}
    catch (std::length_error &e)
    {
        cout<<"Unable to parse the GFT lines"<<endl;
        cerr<<e.what()<<endl;
        return 1;
    }
    catch (std::range_error &e)
    {
        cout<<"Invalid chromosome range"<<endl;
        cerr<<e.what()<<endl;
        return 2;
    }
    catch(std::domain_error &e)
    {
        cout<<"Unable to perform string conversion"<<endl;
        cerr<<e.what()<<endl;
        return 3;
    }
	catch (...)
	{
		return usage(-1, argv[0], "Unknown error");
	}
    return 0;
}

vector<Feature>& extractBlocks(BamAlignment &alignment, vector<Feature> &blocks, unsigned short chr)
{
    auto beg = alignment.CigarData.begin();
    auto end = alignment.CigarData.end();
    coord start = alignment.Position + 1;
    while (beg != end)
    {
        CigarOp current = *(beg++);
        Feature block;
        switch(current.Type)
        {
            case 'M':
            case '=':
            case 'X':
                block.start = start;
                block.chromosome = chr;
                block.end = start + current.Length;
                blocks.push_back(block);
            case 'N':
            case 'D':
                start += current.Length;
            case 'H':
            case 'P':
            case 'S':
            case 'I':
                break;
            default:
                throw std::invalid_argument("Unrecognized Cigar Op");
        }
    }
    return blocks;
}



int usage(int ret, string name, string msg)
{
	cout << endl << "Interval Tree Test" << endl;
	cout << "Usage: $ " << name << " <path to gff file> <path to bam/sam file> [path to bed file for fragment size calculations]" << endl;
	cout << endl << msg << endl;
	return ret;
}


/*
 metrics provider:
 Give access to common resources (genes, transcripts, sequence dictionary, etc)
 Each provider is given an opportunity to parse each alignment in whichever way they see fit
 (Possibly cache common per-bam data based on the first provider that needs it)
 at the end, all providers spit out a report
 */
