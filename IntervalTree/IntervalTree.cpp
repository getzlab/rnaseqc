// IntervalTree.cpp : Defines the entry point for the console application.
//

#include "BED.h"
#include "Metrics.h"
#include <string>
#include <iostream>
#include <exception>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <set>
#include <regex>
#include <list>
#include <ctime>
#include <string>
#include <math.h>
#include <bamtools/api/BamReader.h>
#include <bamtools/api/BamAlignment.h>
#include <args.hxx>
#include <boost/filesystem.hpp>
using namespace std;
using namespace args;
using namespace BamTools;

const int CHIMERIC_DISTANCE = 2000000;
const unsigned int MAX_READ_LENGTH = 100000u;
const unsigned int FRAGMENT_SIZE_SAMPLES = 1000000u;
const unsigned int LOW_QUALITY_READS_THRESHOLD = 255u;
const unsigned int EXON_MISMATCH_THRESHOLD = 6u;
const unsigned int EXON_QUALITY_THRESHOLD = 0u;
const int SPLIT_DISTANCE = 100;
const string NM = "NM";


vector<Feature>& extractBlocks(BamAlignment&, vector<Feature>&, unsigned short);

list<Feature>* intersectBlock(Feature&, list<Feature>&);
void trimFeatures(BamAlignment&, list<Feature> &);

unsigned int fragmentSizeMetrics(unsigned int, map<unsigned short, list<Feature>>*, map<string, string>&, list<long long>&, SamSequenceDictionary&, vector<Feature>&, BamAlignment&);

void exonAlignmentMetrics(map<unsigned short, list<Feature>>&, Metrics&, SamSequenceDictionary&, map<string, double>&, map<string, double>&, vector<Feature>&, BamAlignment&);

int main(int argc, char* argv[])
{
    ArgumentParser parser("rnaSeQC");
    HelpFlag help(parser, "help", "Display this message and quit", {'h', "help"});
    Positional<string> gtfFile(parser, "gtf", "The input GTF file containing features to check the bam against");
    Positional<string> bamFile(parser, "bam", "The input SAM/BAM file containing reads to process");
    ValueFlag<string> bedFile(parser, "bed", "Optional input BED file containing non-overlapping exons used for fragment size calculations", {"bed"});
    Flag debugMode(parser, "debug", "Include values of various internal constants in the output", {'d', "debug"});
    Positional<string> outputDir(parser, "output", "Output directory");
	try
	{
        parser.ParseCLI(argc, argv);
        if (!gtfFile) throw ValidationError("No GTF file provided");
        if (!bamFile) throw ValidationError("No BAM file provided");
        if (!outputDir) throw ValidationError("No output directory provided");
        unsigned int doFragmentSize = 0u;
        ifstream reader(gtfFile.Get());
		cout<<"Reading GTF Features..."<<endl;
        time_t t0, t1, t2;
        time(&t0);
        clock_t start_clock = clock();
        
        //BEGIN GTF PARSING
        
        map<unsigned short, list<Feature>> features;
        map<unsigned short, list<Feature> > *bedFeatures; //for fragment calculations
        Feature line;
        int readLength = 0;
		while ((reader >> line))
		{
            if (line.type == "gene" || line.type == "exon")
            {
                features[line.chromosome].push_back(line);
            }
		}
        for (auto beg = features.begin(); beg != features.end(); ++beg) beg->second.sort(compIntervalStart);
        time(&t1);
        cout << "Finished processing GTF in " << difftime(t1, t0) << " seconds" << endl;
        
        //END GTF PARSING
        //BEGIN BED PARSING
        
        if (bedFile)
        {
            cout << "Parsing BED intervals for fragment size computations..." << endl;
            doFragmentSize = FRAGMENT_SIZE_SAMPLES;
            bedFeatures = new map<unsigned short, list<Feature> >();
            ifstream bedReader(bedFile.Get());
            while (extractBED(bedReader, line)) (*bedFeatures)[line.chromosome].push_back(line);
            bedReader.close();
        }
        
        //END BED PARSING
        //BEGIN BAM PARSING
        
        BamReader bam;
        const string bamFilename = bamFile.Get();
		bam.Open(bamFilename);
        bam.LocateIndex(); //necessary?
		BamAlignment alignment;
		cout<<"Parsing bam..."<<endl;
        SamSequenceDictionary sequences = bam.GetHeader().Sequences;
        Metrics counter;
        map<string, string> fragments;
        list<long long> fragmentSizes;
        map<string, double> geneCoverage, exonCoverage;
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
            
            //BEGIN READ FLAG METRICS
            
            if (!alignment.IsPrimaryAlignment()) counter.increment("Alternative Alignments");
            else if (alignment.IsFailedQC()) counter.increment("Failed Vendor QC");
            else if (alignment.MapQuality < LOW_QUALITY_READS_THRESHOLD) counter.increment("Low quality reads");
            if (alignment.IsPrimaryAlignment() && !alignment.IsFailedQC() /*&& alignment.MapQuality >= 255u*/)
            {
                counter.increment("Total Reads");
                if (!alignment.IsPaired()) counter.increment("Unpaired Reads");
                if (alignment.IsDuplicate()) counter.increment("Duplicate Reads");
                if (alignment.Length > MAX_READ_LENGTH) continue;
                if (alignment.Length > readLength) readLength = alignment.Length;
                if (alignment.IsMapped())
                {
                    counter.increment("Mapped Reads");
                    unsigned int mismatches;
                    alignment.GetTag(NM, mismatches);
                    if (alignment.IsPaired())
                    {
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
                        
                    }
                    counter.increment("Mismatches", mismatches);
                    counter.increment("Total Bases", alignment.Length);
                    if (alignment.IsDuplicate()) counter.increment("Mapped Duplicate Reads");
                    else counter.increment("Mapped Unique Reads");
                    //check length against max read length
                    if (alignment.IsPaired() && alignment.IsMateMapped() && alignment.IsProperPair())
                    {
                        if (alignment.IsFirstMate()) counter.increment("Total Mapped Pairs");
                        if (alignment.RefID != alignment.MateRefID || abs(alignment.Position - alignment.MatePosition) > CHIMERIC_DISTANCE)
                        {
                            counter.increment("Chimeric Pairs");
                            cout << (sequences.Begin()+alignment.RefID)->Name << "\t" << (sequences.Begin()+alignment.MateRefID)->Name << ":\t" << (alignment.RefID == alignment.MateRefID) << endl;
                            continue;
                        }
                    }
                    
                    //END READ FLAG METRICS
                    //BEGIN READ EXON METRICS
                    
                    if (alignment.RefID < 0 || alignment.RefID >= sequences.Size()) cout << "SEQNAME: " << alignment.Name<<endl;
                    else if(mismatches <= EXON_MISMATCH_THRESHOLD && alignment.IsProperPair() && alignment.MapQuality > EXON_QUALITY_THRESHOLD)
                    {
                        vector<Feature> blocks;
                        string chrName = (sequences.Begin()+alignment.RefID)->Name;
                        unsigned short chr = chromosomeMap(chrName);
                        extractBlocks(alignment, blocks, chr);
                        trimFeatures(alignment, features[chr]);
                        exonAlignmentMetrics(features, counter, sequences, geneCoverage, exonCoverage, blocks, alignment);
                        //END READ EXON METRICS
                        //BEGIN READ FRAGMENT METRICS
                        
                        
                        if (doFragmentSize && bedFeatures != nullptr && bedFeatures->find(chr) != bedFeatures->end())
                        {
                            doFragmentSize = fragmentSizeMetrics(doFragmentSize, bedFeatures, fragments, fragmentSizes, sequences, blocks, alignment);
                        }
                    }
                    //End of exon parsing block
                    else counter.increment("Reads excluded from exon counts");
                }
                
            }
            
		}
        
        if (!boost::filesystem::exists(outputDir.Get()))
        {
            boost::filesystem::create_directories(outputDir.Get());
        }
        time(&t2);
        cout<< "Time Elapsed: " << difftime(t2, t1) << "; Alignments processed: " << alignmentCount << endl;
        cout << "Total runtime: " << difftime(t2, t0) << "; Total CPU Time: " << (clock() - start_clock)/CLOCKS_PER_SEC << endl;
        
        cout << "Generating report" << endl;
        
        unsigned int genesDetected = 0;
        
        ofstream geneReport(outputDir.Get()+"/geneReport.tsv");
        for(auto gene = geneCoverage.begin(); gene != geneCoverage.end(); ++gene)
        {
            geneReport << gene->first << "\t" << gene->second << endl;
            if (gene->second >= 5.0) ++genesDetected;
        }
        geneReport.close();
        
        ofstream exonReport(outputDir.Get()+"/exonReport.tsv");
        for(auto exon = exonCoverage.begin(); exon != exonCoverage.end(); ++exon)
        {
            exonReport << exon->first << "\t" << exon->second << endl;
        }
        exonReport.close();
        
        ofstream output(outputDir.Get()+"/report.tsv");
        output << "Sample\t" << gtfFile.Get() << endl;
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
        output << "Intronic Rate\t" << counter.frac("Intronic", "Mapped Reads") << endl;
        output << "Intergenic Rate\t" << counter.frac("Intergenic Reads", "Mapped Reads") << endl;
        output << "Intragenic Rate\t" << counter.frac("Intragenic Reads", "Mapped Reads") << endl;
        output << "rRNA Rate\t" << counter.frac("rRNA", "Mapped Reads") << endl;
        output << "End 1 Sense Rate\t" << (double) counter.get("End 1 Sense") / (counter.get("End 1 Sense") + counter.get("End 1 Antisense")) << endl;
        output << "End 2 Sense Rate\t" << (double) counter.get("End 2 Sense") / (counter.get("End 2 Sense") + counter.get("End 2 Antisense")) << endl;
        output << "Avg. Blocks per Read\t" << counter.frac("Alignment Blocks", "Mapped Reads") << endl;
        
        output << counter;
        output << "Read Length\t" << readLength << endl;
        output << "Genes Detected\t" << genesDetected << endl;
        output << "Total Reads (incl. supplemental and failed reads)\t" << alignmentCount << endl;
        
        if (fragmentSizes.size())
        {
            ofstream checksum(outputDir+"/checksum.txt");
            for(auto fragment = fragmentSizes.begin(); fragment != fragmentSizes.end(); ++fragment)
            {
                checksum << abs(*fragment) << endl;
            }
            checksum.close();
            
            double fragmentAvg = 0.0, fragmentStd = 0.0, fragmentMed = 0.0, fragmentMedDev = 0.0;
            fragmentSizes.sort();
            auto median = fragmentSizes.begin();
            for (int midpoint = fragmentSizes.size() / 2; midpoint > 0; --midpoint) ++median;
            fragmentMed = (double) *median;
            vector<double> deviations;
            double size = (double) fragmentSizes.size();
            for(auto fragment = fragmentSizes.begin(); fragment != fragmentSizes.end(); ++fragment)
            {
                fragmentAvg += (double) abs(*fragment) / size;
                deviations.push_back(fabs((double) (*fragment) - fragmentMed));
            }
            sort(deviations.begin(), deviations.end());
            fragmentMedDev = (double) deviations[deviations.size()/2] * 1.4826;
            for(auto fragment = fragmentSizes.begin(); fragment != fragmentSizes.end(); ++fragment)
            {
                fragmentStd += pow((double) (*fragment) - fragmentAvg, 2.0) / size;
            }
            fragmentStd = pow(fragmentStd, 0.5);
            
            output << "Average Fragment Length\t" << fragmentAvg << endl;
            output << "Fragment Length Std\t" << fragmentStd << endl;
            output << "Fragment Length MAD_Std\t" << fragmentMedDev << endl;
        }
        
        if (debugMode.Get())
        {
            output << "[DEBUG]CHIMERIC_DISTANCE\t" << CHIMERIC_DISTANCE << endl;
            output << "[DEBUG]MAX_READ_LENGTH\t" << MAX_READ_LENGTH << endl;
            output << "[DEBUG]FRAGMENT_SIZE_SAMPLES\t" << FRAGMENT_SIZE_SAMPLES << endl;
            output << "[DEBUG]LOW_QUALITY_READS_THRESHOLD\t" << LOW_QUALITY_READS_THRESHOLD << endl;
            output << "[DEBUG]EXON_MISMATCH_THRESHOLD\t" << EXON_MISMATCH_THRESHOLD << endl;
            output << "[DEBUG]EXON_QUALITY_THRESHOLD\t" << EXON_QUALITY_THRESHOLD << endl;
            output << "[DEBUG]SPLIT_DISTANCE\t" << SPLIT_DISTANCE << endl;
            
        }
        
        output.close();
        
//        Note  Adrenal
//        Chimeric Pairs  0 <----
//        Fragment Length StdDev  273
//        Estimated Library Size  0
//        Fragment Length Mean  259


	}
    catch (args::Help)
    {
        cout << parser;
        return 4;
    }
    catch (args::ParseError &e)
    {
        cout << parser << endl;
        cout << "Argument parsing error: " << e.what() << endl;
        return 5;
    }
    catch (args::ValidationError &e)
    {
        cout << parser << endl;
        cout << "Argument validation error: " << e.what() << endl;
        return 6;
    }
    catch (std::invalid_argument &e)
    {
        cout << "Invalid argument type provided: " << e.what() << endl;
        return 7;
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
        cout << parser << endl;
        cout << "Unknown error" << endl;
        return -1;
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

void trimFeatures(BamAlignment &alignment, list<Feature> &features)
{
    //trim intervals upstream of this block
    while (features.size() && features.front().end < alignment.Position) features.pop_front();
}

list<Feature>* intersectBlock(Feature &block, list<Feature> &features)
{
    list<Feature> *output = new list<Feature>();
    for (auto current = features.begin(); current != features.end() && current->start <= block.end; ++current)
    {
        if (intersectInterval(block, *current)) output->push_back(*current);
    }
    return output;
}




void exonAlignmentMetrics(map<unsigned short, list<Feature>> &features, Metrics &counter, SamSequenceDictionary &sequences, map<string, double> &geneCoverage, map<string, double> &exonCoverage, vector<Feature> &blocks, BamAlignment &alignment)
{
    counter.increment("Reads used for Intron/Exon counts");
    //map the 0...N RefID to an internal chromosome number consistent with the GTF input
    string chrName = (sequences.Begin()+alignment.RefID)->Name;
    unsigned short chr = chromosomeMap(chrName);
    
    bool split = false;
    long long lastEnd = -1;
    for(auto block = blocks.begin(); block != blocks.end(); ++block)
    {
        if (lastEnd > 0 && !split) split = (block->start - lastEnd) > SPLIT_DISTANCE;
        lastEnd = block->end;
    }
    
    //No longer skipping this step since trimFeatures inadvertently ensures it's always activated
    counter.increment("Alignment Blocks", blocks.size());
    Feature current;
    current.start = alignment.Position + 1; //0-based + 1 == 1-based
    current.end = alignment.GetEndPosition(); //0-based, open == 1-based, closed
    
    bool exonic = false, intragenic = false, intronic = false, transcriptPlus = false, transcriptMinus = false, ribosomal = false;
    string exonName;
    vector<set<string> > genes;
    Collector exonCoverageCollector(&exonCoverage);
    
    for (auto block = blocks.begin(); block != blocks.end(); ++block)
    {
        genes.push_back(set<string>());
        list<Feature> *results = intersectBlock(*block, features[chr]);
        auto beg = results->begin();
        auto end = results->end();
        while (beg != end)
        {
            if (beg->strand == 1) transcriptPlus = true;
            else if (beg->strand == -1) transcriptMinus = true;
            if (beg->type == "exon")
            {
                intragenic = true;
                double tmp = ((double) block->end - block->start) / alignment.Length;
                //exonCoverage[beg->feature_id] += tmp;
                exonCoverageCollector.add(beg->gene_id, beg->feature_id, tmp);
                genes.rbegin()->insert(beg->gene_id);
            }
            else if (beg->type == "gene")
            {
                intragenic = true;
            }
            if (beg->transcript_type == "rRNA" || beg->transcript_type == "Mt_rRNA") ribosomal = true;
            ++beg;
        }
        delete results;
    }
    bool doExonMetrics = false;
    if (genes.size() > 1)
    {
        set<string> last = genes.front();
        for (int i = 1; i < genes.size(); ++i)
        {
            set<string> tmp;
            set_intersection(last.begin(), last.end(), genes[i].begin(), genes[i].end(), inserter(tmp, tmp.begin()));
            last = tmp;
        }
        for (auto gene = last.begin(); gene != last.end(); ++gene)
        {
            geneCoverage[*gene]++;
            exonCoverageCollector.collect(*gene);
            doExonMetrics = true;
        }
    }
    else if (genes.size() == 1) for (auto gene = genes[0].begin(); gene != genes[0].end(); ++gene)
    {
        geneCoverage[*gene]++;
        exonCoverageCollector.collect(*gene);
        doExonMetrics = true;
    }
    if (!exonCoverageCollector.isDirty())
    {
        counter.increment("Intronic");
    }
    else if (doExonMetrics)
    {
        counter.increment("Exonic");
        if (split) counter.increment("Split Reads");
    }
    if (intragenic) counter.increment("Intragenic Reads");
    else counter.increment("Intergenic Reads");
    if (ribosomal) counter.increment("rRNA");
    
    if (transcriptMinus ^ transcriptPlus && alignment.IsPaired())
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

unsigned int fragmentSizeMetrics(unsigned int doFragmentSize, map<unsigned short, list<Feature>> *bedFeatures, map<string, string> &fragments, list<long long> &fragmentSizes, SamSequenceDictionary &sequences, vector<Feature> &blocks, BamAlignment &alignment)
{
    string chrName = (sequences.Begin()+alignment.RefID)->Name;
    unsigned short chr = chromosomeMap(chrName);
    bool firstBlock = true, sameExon = true;
    string exonName = "";
    trimFeatures(alignment, (*bedFeatures)[chr]);
    for (auto block = blocks.begin(); sameExon && block != blocks.end(); ++block)
    {
        list<Feature> *results = intersectBlock(*block, (*bedFeatures)[chr]);
        if (results->size() == 1)
        {
            if (alignment.IsPaired())
            {
                if (firstBlock) exonName = results->begin()->feature_id;
                else if (exonName != results->begin()->feature_id)
                {
                    sameExon = false;
                    break;
                }
                
            }
        }
        else sameExon = false;
        delete results;
        firstBlock = false;
    }
    if (sameExon && exonName.size())
    {
        auto fragment = fragments.find(alignment.Name);
        if (fragment == fragments.end()) //first time we've encountered a read in this pair
        {
            fragments[alignment.Name] = exonName;
        }
        else if (exonName == fragments[alignment.Name])
        {
            //This pair is useable for fragment statistics:  both pairs fully aligned to the same exon
            fragmentSizes.push_back(abs(alignment.InsertSize));
            fragments.erase(fragment);
            --doFragmentSize;
            if (!doFragmentSize)
            {
                cout << "Completed taking fragment size samples" << endl;
                delete bedFeatures;
            }
        }
    }
    return doFragmentSize;
}
