// IntervalTree.cpp : Defines the entry point for the console application.

//Include headers
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
//TODO: ensure full intersection of read and gene (possibly block and exon as well)

template <typename T>
bool isIn(set<T> &s, T member)
{
    return s.size() ? s.find(member) != s.end() : false;
}

const string NM = "NM";
bool debugging = false;
bool dumbMode = false;

//Utility functions
unsigned int extractBlocks(BamAlignment&, vector<Feature>&, unsigned short);
list<Feature>* intersectBlock(Feature&, list<Feature>&);
void trimFeatures(BamAlignment&, list<Feature> &);

//Metrics functions
unsigned int fragmentSizeMetrics(unsigned int, map<unsigned short, list<Feature>>*, map<string, string>&, list<long long>&, SamSequenceDictionary&, vector<Feature>&, BamAlignment&);

void exonAlignmentMetrics(unsigned int, map<unsigned short, list<Feature>>&, Metrics&, SamSequenceDictionary&, map<string, double>&, map<string, double>&, vector<Feature>&, BamAlignment&, unsigned int);

void legacyExonAlignmentMetrics(unsigned int, map<unsigned short, list<Feature>>&, Metrics&, SamSequenceDictionary&, map<string, double>&, map<string, double>&, vector<Feature>&, BamAlignment&, unsigned int);

int main(int argc, char* argv[])
{
    //Set up command line syntax
    ArgumentParser parser("rnaSeQC");
    HelpFlag help(parser, "help", "Display this message and quit", {'h', "help"});
    Positional<string> gtfFile(parser, "gtf", "The input GTF file containing features to check the bam against");
    Positional<string> bamFile(parser, "bam", "The input SAM/BAM file containing reads to process");
    ValueFlag<string> bedFile(parser, "bed", "Optional input BED file containing non-overlapping exons used for fragment size calculations", {"bed"});
    ValueFlag<int> chimericDistance(parser, "dist", "Set the maximum accepted distance between read mates.  Mates beyond this distance will be counted as chimeric pairs", {"chimeric-distance"});
    ValueFlag<unsigned int> maxReadLength(parser, "length", "Set the maximum accepted length.  Reads longer than this threshold are discarded", {"read-length"});
    ValueFlag<unsigned int> fragmentSamples(parser, "samples", "Set the number of samples to take when computing fragment sizes.  Requires the --bed argument", {"samples"});
    ValueFlag<unsigned int> lowQualityThreshold(parser, "quality", "Set the lower bound on read quality. Reads below this number are counted as low quality BUT ARE STILL USED IN COUNTS.  See --exon-quality to discard reads based on quality", {"low-quality"});
    ValueFlag<unsigned int> exonQualityThreshold(parser,"quality", "Set the lower bound on read quality for exon coverage counting.  Reads below this number are excluded from coverage metrics", {"exon-quality"});
    ValueFlag<unsigned int> exonMismatchThreshold(parser, "mismatches", "Set the maximum number of allowed mismatches between a read and the reference sequence.  Reads with more than this number of mismatches are excluded from coverage metrics", {"exon-mismatch"});
    ValueFlag<int> splitDistance(parser, "dist", "Set the maximum distance between aligned blocks of a read.  Reads with aligned blocks separated by more than this distance are counted as split reads, BUT ARE STILL USED IN COUNTS", {"split-distance"});
    Flag debugMode(parser, "debug", "Include values of various internal constants in the output", {'d', "debug"});
    Flag LegacyMode(parser, "legacy", "Use legacy exon metrics.  Matches output of RNA-SeQC 1.1.6", {"legacy"});
    Flag dumbModeQ(parser, "dumb", "Use dumb mode", {"dumb"});
    Positional<string> outputDir(parser, "output", "Output directory");
	try
	{
        //parse and validate the command line arguments
        parser.ParseCLI(argc, argv);
        if (!gtfFile) throw ValidationError("No GTF file provided");
        if (!bamFile) throw ValidationError("No BAM file provided");
        if (!outputDir) throw ValidationError("No output directory provided");
        
        const int CHIMERIC_DISTANCE = chimericDistance ? chimericDistance.Get() : 2000000;
        const unsigned int MAX_READ_LENGTH = maxReadLength ? maxReadLength.Get() : 100000u;
        const unsigned int FRAGMENT_SIZE_SAMPLES = fragmentSamples ? fragmentSamples.Get() : 1000000u;
        const unsigned int LOW_QUALITY_READS_THRESHOLD = lowQualityThreshold ? lowQualityThreshold.Get() : 255u;
        const unsigned int EXON_MISMATCH_THRESHOLD = 6u;
        const unsigned int EXON_QUALITY_THRESHOLD = exonQualityThreshold ? exonQualityThreshold.Get() : 0u;
        const int SPLIT_DISTANCE = 100;
        debugging = debugMode.Get();
        dumbMode = dumbModeQ.Get();
        
        //Define variables
        Metrics counter; //main tracker for various metrics
        int readLength = 0; //longest read encountered so far
        time_t t0, t1, t2; //various timestamps to record execution time
        clock_t start_clock = clock(); //timer used to compute CPU time
        BamReader bam;
        const string bamFilename = bamFile.Get();
        SamSequenceDictionary sequences; //for chromosome lookup
        map<string, double> geneCoverage, exonCoverage; //counters for read coverage of genes and exons
        unsigned long long alignmentCount = 0ull; //count of how many alignments we've seen so far
        
        //fragment size variables
        unsigned int doFragmentSize = 0u; //count of remaining fragment size samples to record
        map<unsigned short, list<Feature>> features; //map of chr -> genes/exons; parsed from GTF
        map<unsigned short, list<Feature> > *bedFeatures; //similar map, but parsed from BED for fragment sizes only
        map<string, string> fragments; //Map of alignment name -> exonID to ensure mates map to the same exon for
        list<long long> fragmentSizes; //list of fragment size samples taken so far
        
        //Parse the GTF and extract features
        {
            Feature line; //current feature being read from the gtf
            ifstream reader(gtfFile.Get());
            
            cout<<"Reading GTF Features..."<<endl;
            time(&t0);
            while ((reader >> line))
            {
                //Just keep genes and exons.  We don't care about transcripts or any other feature types
                if (line.type == "gene" || line.type == "exon")
                {
                    features[line.chromosome].push_back(line);
                }
            }
        }
		//ensure that the features are sorted.  This MUST be true for the exon alignment metrics
        for (auto beg = features.begin(); beg != features.end(); ++beg) beg->second.sort(compIntervalStart);
        time(&t1); //record the time taken to parse the GTF
        cout << "Finished processing GTF in " << difftime(t1, t0) << " seconds" << endl;
        
        if (bedFile) //If we were given a BED file, parse it for fragment size calculations
        {
             Feature line; //current feature being read from the bed
            cout << "Parsing BED intervals for fragment size computations..." << endl;
            doFragmentSize = FRAGMENT_SIZE_SAMPLES;
            bedFeatures = new map<unsigned short, list<Feature> >();
            ifstream bedReader(bedFile.Get());
            //extract each line of the bed and insert it into the bedFeatures map
            while (extractBED(bedReader, line)) (*bedFeatures)[line.chromosome].push_back(line);
            bedReader.close();
        }
        
        //Begin parsing the bam.  Each alignment is run through various sets of metrics
        {
            BamAlignment alignment; //current bam alignment
            
            time_t report_time; //used to ensure that stdout isn't spammed if the program runs super fast
            bam.Open(bamFilename);
            bam.LocateIndex(); //load in the index, if found.  Slightly improves IO perf
            sequences = bam.GetHeader().Sequences; //read the sequence dictionary from the header
            cout<<"Parsing bam..."<<endl;
            time(&report_time);
            time(&t2);
            while (bam.GetNextAlignment(alignment))
            {
                //try to print an update to stdout every 250,000 reads, but no more than once every 10 seconds
                ++alignmentCount;
                if (alignmentCount % 250000 == 0) time(&t2);
                if (difftime(t2, report_time) >= 10)
                {
                    time(&report_time);
                    cout << "Time elapsed: " << difftime(t2, t1) << "; Alignments processed: " << alignmentCount << endl;
                }
                //count metrics based on basic read data
                if (!alignment.IsPrimaryAlignment()) counter.increment("Alternative Alignments");
                else if (alignment.IsFailedQC()) counter.increment("Failed Vendor QC");
                else if (alignment.MapQuality < LOW_QUALITY_READS_THRESHOLD) counter.increment("Low quality reads");
                if (alignment.IsPrimaryAlignment() && !alignment.IsFailedQC() /*&& alignment.MapQuality >= 255u*/)
                {
                    counter.increment("Total Reads");
                    if (!alignment.IsPaired()) counter.increment("Unpaired Reads");
                    if (alignment.IsDuplicate()) counter.increment("Duplicate Reads");
//                    if (alignment.Length > MAX_READ_LENGTH) continue;
//                    if (alignment.Length > readLength) readLength = alignment.Length;
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
                        unsigned int alignmentSize = alignment.GetEndPosition() - alignment.Position + 1;
                        if (alignmentSize > MAX_READ_LENGTH) continue;
                        if (alignmentSize > readLength) readLength = alignment.Length;
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
                        
                        //now record intron/exon metrics by intersecting filtered reads with the list of features
                        if (alignment.RefID < 0 || alignment.RefID >= sequences.Size())
                        {
                            //The read had an unrecognized RefID (one not defined in the bam's header)
                            cout << "Unrecognized RefID on alignment: " << alignment.Name<<endl;
                        }
                        else if(mismatches <= EXON_MISMATCH_THRESHOLD && alignment.IsProperPair() && alignment.MapQuality >= EXON_QUALITY_THRESHOLD)
                        {
                            vector<Feature> blocks;
                            string chrName = (sequences.Begin()+alignment.RefID)->Name;
                            unsigned short chr = chromosomeMap(chrName); //parse out a chromosome shorthand
                            unsigned int length = extractBlocks(alignment, blocks, chr); //extract each cigar block from the alignment
                            trimFeatures(alignment, features[chr]); //drop features that appear before this read
                            
                            //run the read through exon metrics
                            //exonAlignmentMetrics(SPLIT_DISTANCE, features, counter, sequences, geneCoverage, exonCoverage, blocks, alignment, length);
                            legacyExonAlignmentMetrics(SPLIT_DISTANCE, features, counter, sequences, geneCoverage, exonCoverage, blocks, alignment, length);
                            
                            //if fragment size calculations were requested, we still have samples to take, and the chromosome exists within the provided bed
                            if (alignment.IsPaired() && doFragmentSize && bedFeatures != nullptr && bedFeatures->find(chr) != bedFeatures->end())
                            {
                                doFragmentSize = fragmentSizeMetrics(doFragmentSize, bedFeatures, fragments, fragmentSizes, sequences, blocks, alignment);
                            }
                        }
                        else counter.increment("Reads excluded from exon counts");
                    }
                    
                }
                
            } //end of bam alignment loop
        } //end of bam alignment scope
		
        //use boost to ensure that the output directory exists before the metrics are dumped to it
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
        geneReport << fixed;
        //iterate over every gene with coverage reported.  If it had at leat 5 reads, also count it as 'detected'
        for(auto gene = geneCoverage.begin(); gene != geneCoverage.end(); ++gene)
        {
            geneReport << gene->first << "\t" << gene->second << endl;
            if (gene->second >= 5.0) ++genesDetected;
        }
        geneReport.close();
        
        ofstream exonReport(outputDir.Get()+"/exonReport.tsv");
        exonReport << fixed;
        //iterate over every exon with coverage reported
        for(auto exon = exonCoverage.begin(); exon != exonCoverage.end(); ++exon)
        {
            exonReport << exon->first << "\t" << exon->second << endl;
        }
        exonReport.close();
        
        ofstream output(outputDir.Get()+"/report.tsv");
        //output rates and other fractions to the report
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
        //automatically dump the raw counts of all metrics to the file
        output << counter;
        //append metrics that were manually tracked
        output << "Read Length\t" << readLength << endl;
        output << "Genes Detected\t" << genesDetected << endl;
        output << "Total Reads (incl. supplemental and failed reads)\t" << alignmentCount << endl;
        
        if (fragmentSizes.size())
        {
            //If any fragment size samples were taken, also generate a fragment size report
            fragmentSizes.sort();
            
            
            
            double fragmentAvg = 0.0, fragmentStd = 0.0, fragmentMed = 0.0, fragmentMedDev = 0.0;
            //I can't believe I have to do this, but G++98's list.size() has to WALK THE ENTIRE LIST FOR SOME REASON
            double size = (double) fragmentSizes.size();
            vector<double> deviations; //list of recorded deviations from the median
            auto median = fragmentSizes.begin(); //reference the median value.  We have to walk the list to get here
            for (int midpoint = fragmentSizes.size() / 2; midpoint > 0; --midpoint) ++median;
            fragmentMed = (double) *median; //save the median value
            ofstream fragmentList(outputDir.Get()+"/fragmentSizes.txt"); //raw list of each fragment size recorded
            for(auto fragment = fragmentSizes.begin(); fragment != fragmentSizes.end(); ++fragment)
            {
                fragmentList << abs(*fragment) << endl; //record the fragment size into the output list
                fragmentAvg += (double) abs(*fragment) / size; //add this fragment's size to the mean
                deviations.push_back(fabs((double) (*fragment) - fragmentMed)); //record this fragment's deviation
            }
            fragmentList.close();
            sort(deviations.begin(), deviations.end()); //for the next line to work, we have to sort
            //now compute the median absolute deviation, an estimator for standard deviation
            fragmentMedDev = (double) deviations[deviations.size()/2] * 1.4826;
            //we have to iterate again now for the standard deviation calculation, now that we know the mean
            for(auto fragment = fragmentSizes.begin(); fragment != fragmentSizes.end(); ++fragment)
            {
                fragmentStd += pow((double) (*fragment) - fragmentAvg, 2.0) / size;
            }
            fragmentStd = pow(fragmentStd, 0.5); //compute the standard deviation
            
            output << "Average Fragment Length\t" << fragmentAvg << endl;
            output << "Fragment Length Std\t" << fragmentStd << endl;
            output << "Fragment Length MAD_Std\t" << fragmentMedDev << endl;
        }
        
        if (debugMode.Get())
        {
            //output the constants used
            output << "[DEBUG]CHIMERIC_DISTANCE\t" << CHIMERIC_DISTANCE << endl;
            output << "[DEBUG]MAX_READ_LENGTH\t" << MAX_READ_LENGTH << endl;
            output << "[DEBUG]FRAGMENT_SIZE_SAMPLES\t" << FRAGMENT_SIZE_SAMPLES << endl;
            output << "[DEBUG]LOW_QUALITY_READS_THRESHOLD\t" << LOW_QUALITY_READS_THRESHOLD << endl;
            output << "[DEBUG]EXON_MISMATCH_THRESHOLD\t" << EXON_MISMATCH_THRESHOLD << endl;
            output << "[DEBUG]EXON_QUALITY_THRESHOLD\t" << EXON_QUALITY_THRESHOLD << endl;
            output << "[DEBUG]SPLIT_DISTANCE\t" << SPLIT_DISTANCE << endl;
            
        }
        
        output.close();

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

unsigned int extractBlocks(BamAlignment &alignment, vector<Feature> &blocks, unsigned short chr)
{
    //parse the cigar string and populate the provided vector with each block of the read
    auto beg = alignment.CigarData.begin();
    auto end = alignment.CigarData.end();
    coord start = alignment.Position + 1;
    unsigned int alignedSize = 0;
    while (beg != end)
    {
        CigarOp current = *(beg++);
        Feature block;
        switch(current.Type)
        {
            case 'M':
            case '=':
            case 'X':
                //M, =, and X blocks are aligned, so push back this block
                block.start = start;
                block.chromosome = chr;
                block.end = start + current.Length; //1-based, closed
                blocks.push_back(block);
                alignedSize += current.Length;
            case 'N':
            case 'D':
                //M, =, X, N, and D blocks all advance the start position of the next block
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
    return alignedSize;
}

void trimFeatures(BamAlignment &alignment, list<Feature> &features)
{
    //trim intervals upstream of this block
    //Since alignments are sorted, if an alignment occurs beyond any features, these features can be dropped
    while (features.size() && features.front().end < alignment.Position) features.pop_front();
}

list<Feature>* intersectBlock(Feature &block, list<Feature> &features)
{
    list<Feature> *output = new list<Feature>();
    //since we've trimmed the beginning of the features, we start from the new beginning here
    //There should be little overhead (at most ~1 gene worth of exons on either end of the block)
    for (auto current = features.begin(); current != features.end() && current->start <= block.end; ++current)
    {
        //check that the current feature actually intersects the current block of the current alignment
        if (intersectInterval(block, *current)) output->push_back(*current);
    }
    return output;
}

void legacyExonAlignmentMetrics(unsigned int SPLIT_DISTANCE, map<unsigned short, list<Feature>> &features, Metrics &counter, SamSequenceDictionary &sequences, map<string, double> &geneCoverage, map<string, double> &exonCoverage, vector<Feature> &blocks, BamAlignment &alignment, unsigned int length)
{
    string chrName = (sequences.Begin()+alignment.RefID)->Name;
    unsigned short chr = chromosomeMap(chrName); //generate the chromosome shorthand name
    
    //check for split reads by iterating over all the blocks of this read
    bool split = false;
    long long lastEnd = -1; // used for split read detection
    for(auto block = blocks.begin(); block != blocks.end(); ++block)
    {
        if (lastEnd > 0 && !split) split = (block->start - lastEnd) > SPLIT_DISTANCE;
        lastEnd = block->end;
    }
    
    counter.increment("Alignment Blocks", blocks.size());
    counter.increment("Reads used for Intron/Exon counts");
    //Bamtools uses 0-based indexing because it's the 'norm' in computer science, even though bams are 1-based
    Feature current; //current block of the alignment (used while iterating)
    current.start = alignment.Position+1; //0-based + 1 == 1-based
    current.end = alignment.GetEndPosition(); //0-based, open == 1-based, closed
    
    vector<set<string> > genes; //each set is the set of genes intersected by the current block (one set per block)
    Collector exonCoverageCollector(&exonCoverage); //Collects coverage counts for later (counts may be discarded)
    bool intragenic = false, transcriptPlus = false, transcriptMinus = false, ribosomal = false, doExonMetrics = false; //various booleans for keeping track of the alignment
    
    for (auto block = blocks.begin(); block != blocks.end(); ++block)
    {
        bool blockWasIntragenic = false; //legacy
        genes.push_back(set<string>()); //create a new set for this block
        list<Feature> *results = intersectBlock(*block, features[chr]); //grab the list of intersecting features
        for (auto result = results->begin(); result != results->end(); ++result)
        {
            if (result->strand == 1) transcriptPlus = true;
            else if (result->strand == -1) transcriptMinus = true;
            if (result->type == "exon")
            {
 
                int intersectionSize = partialIntersect(*result, *block);
                //check that this block fully overlaps the feature
                //(if any bases of the block don't overlap, then the read is discarded)
                /*/ //toggle here to switch between genes and exon
#define gene_mode
                double tmp = (double) intersectionSize / length;
                exonCoverageCollector.add(result->gene_id, result->feature_id, tmp);
                /*/
                if (intersectionSize == block->end - block->start)
                {
                    genes.rbegin()->insert(result->gene_id);
                    blockWasIntragenic = true; //legacy
                    double tmp = (double) intersectionSize / length;
                    exonCoverageCollector.add(result->gene_id, result->feature_id, tmp);
                }
                /**/
                //probably look into ensuring that the intersection size spans this whole block
                
            }
            else if (result->type == "gene")
            {
                intragenic = true;
#ifdef gene_mode
                int intersectionSize = partialIntersect(*result, *block);
                if (intersectionSize == block->end - block->start)
                {
                    //store the exon split dosage coverage in the collector for now
                    genes.rbegin()->insert(result->gene_id);
                    blockWasIntragenic = true; //legacy
                }
#endif
                //we don't record the gene name here because in terms of gene coverage and detection, we only care about exons (apparently not)
            }
            if (result->transcript_type == "rRNA" || result->transcript_type == "Mt_rRNA") ribosomal = true;
        }
        delete results; //clean up dynamic allocation
        //legacy as fuck
        if (split && !blockWasIntragenic && genes.size() == 1 && blocks.size() > 1)
        {
            genes.pop_back();
            split = false;
        }
    }
    set<string> seen;
    /*doExonMetrics = true;
    //toggle the following line to switch between all blocks or any blocks mode
    if (split)
    {
        //legacy code handles split reads differently
        for (auto run = genes.begin(); run != genes.end(); ++run) if(!run->size()) doExonMetrics = false;
        if (doExonMetrics)
        {
            doExonMetrics = false;
            for (auto run = genes.begin(); run != genes.end(); ++run)
            {
                for (auto gene = run->begin(); gene != run->end(); ++gene)
                {
                    if (!isIn(seen, *gene))
                    {
                        //probably just easier to record the intragenic state of each alignment block
                        //cout << alignment.Name << " -> " << *gene << endl;
                        if (exonCoverageCollector.queryGene(*gene))
                        {
                            geneCoverage[*gene]++;
                            //cout << "+" << endl;
                        }
                        exonCoverageCollector.collect(*gene);
                        //cout << endl;
                        doExonMetrics = true;
                        seen.insert(*gene);
                    }
                }
            }
        }
    }
    else
    {
        for (auto gene = genes.rbegin()->begin(); gene != genes.rbegin()->end(); ++gene)
        {
            //probably just easier to record the intragenic state of each alignment block
            //cout << alignment.Name << " -> " << *gene << endl;
            if (exonCoverageCollector.queryGene(*gene))
            {
                geneCoverage[*gene]++;
                //cout << "+" << endl;
            }
            exonCoverageCollector.collectSingle(*gene);
            //cout << endl;
            doExonMetrics = true;
        }
    }*/
    for (auto gene = genes.rbegin()->begin(); gene != genes.rbegin()->end(); ++gene)
    {
        //probably just easier to record the intragenic state of each alignment block
        //cout << alignment.Name << " -> " << *gene << endl;
        if (exonCoverageCollector.queryGene(*gene))
        {
            geneCoverage[*gene]++;
            //cout << "+" << endl;
            //cout << "Dumb: " << alignment.Name << " " << *gene << endl;
        }
        exonCoverageCollector.collectSingle(alignment.Name, *gene);
        //cout << endl;
        doExonMetrics = true;
    }
    
    if (!exonCoverageCollector.isDirty()) //a.k.a: No exons were detected at all on any block of the read
    {
        if (intragenic)
        {
            counter.increment("Intronic");
            counter.increment("Intragenic Reads");
        }
        else counter.increment("Intergenic Reads");
    }
    else if (doExonMetrics) //if exons were detected and at least one exon ended up being collected, we count this as exonic
    {
        counter.increment("Exonic");
        counter.increment("Intragenic Reads");
        if (split) counter.increment("Split Reads");
    }
    else
    {
        //It's unclear how to properly classify these reads
        //They had exon coverage, but aligned to multiple genes
        //Any exon and gene coverage they had was discarded and not recorded
        counter.increment("Intron/Exon disqualified reads");
    }
    if (ribosomal) counter.increment("rRNA");
    //also record strandedness counts
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

void exonAlignmentMetrics(unsigned int SPLIT_DISTANCE, map<unsigned short, list<Feature>> &features, Metrics &counter, SamSequenceDictionary &sequences, map<string, double> &geneCoverage, map<string, double> &exonCoverage, vector<Feature> &blocks, BamAlignment &alignment, unsigned int length)
{
    string chrName = (sequences.Begin()+alignment.RefID)->Name;
    unsigned short chr = chromosomeMap(chrName); //generate the chromosome shorthand name
    
    //check for split reads by iterating over all the blocks of this read
    bool split = false;
    long long lastEnd = -1; // used for split read detection
    for(auto block = blocks.begin(); block != blocks.end(); ++block)
    {
        if (lastEnd > 0 && !split) split = (block->start - lastEnd) > SPLIT_DISTANCE;
        lastEnd = block->end;
    }
    
    counter.increment("Alignment Blocks", blocks.size());
    counter.increment("Reads used for Intron/Exon counts");
    //Bamtools uses 0-based indexing because it's the 'norm' in computer science, even though bams are 1-based
    Feature current; //current block of the alignment (used while iterating)
    current.start = alignment.Position+1; //0-based + 1 == 1-based
    current.end = alignment.GetEndPosition(); //0-based, open == 1-based, closed
    
    vector<set<string> > genes; //each set is the set of genes intersected by the current block (one set per block)
    Collector exonCoverageCollector(&exonCoverage); //Collects coverage counts for later (counts may be discarded)
    bool intragenic = false, transcriptPlus = false, transcriptMinus = false, ribosomal = false, doExonMetrics = false; //various booleans for keeping track of the alignment
    
    bool legacyFirstBlock = false;
    //toggle next line for legacy block
    //legacyFirstBlock = true;
    for (auto block = blocks.begin(); block != blocks.end(); ++block)
    {
        bool blockWasIntragenic = false; //legacy
        genes.push_back(set<string>()); //create a new set for this block
        list<Feature> *results = intersectBlock(*block, features[chr]); //grab the list of intersecting features
        for (auto result = results->begin(); result != results->end(); ++result)
        {
            if (result->strand == 1) transcriptPlus = true;
            else if (result->strand == -1) transcriptMinus = true;
            if (result->type == "exon")
            {
                int intersectionSize = partialIntersect(*result, *block);
                //check that this block fully overlaps the feature
                //(if any bases of the block don't overlap, then the read is discarded)

                
                /**/ //Toggle here to switch between exon/gene requirements
                double tmp = (double) intersectionSize / length;
                exonCoverageCollector.add(result->gene_id, result->feature_id, tmp);
#define Silly_toggle
                /*/
                if (intersectionSize == block->end - block->start)
                {
                    //store the exon split dosage coverage in the collector for now
                    genes.rbegin()->insert(result->gene_id);
                    double tmp = (double) intersectionSize / length;
                    exonCoverageCollector.add(result->gene_id, result->feature_id, tmp);
                    blockWasIntragenic = true; //legacy

                }
                /**/
                
            }
            else if (result->type == "gene")
            {
                intragenic = true;
#ifdef Silly_toggle
                int intersectionSize = partialIntersect(*result, *block);
                if (intersectionSize == block->end - block->start)
                {
                    //store the exon split dosage coverage in the collector for now
                    genes.rbegin()->insert(result->gene_id);
                    blockWasIntragenic = true; //legacy
                }
#endif
                //we don't record the gene name here because in terms of gene coverage and detection, we only care about exons (apparently not)
            }
            if (result->transcript_type == "rRNA" || result->transcript_type == "Mt_rRNA") ribosomal = true;
        }
        delete results; //clean up dynamic allocation
        //legacy as fuck
        if ((split || legacyFirstBlock) && !blockWasIntragenic && genes.size() == 1 && blocks.size() > 1)
        {
            genes.pop_back();
            if (split) split = false;
        }
        legacyFirstBlock = false;
    }
    /**/ //toggle here to switch between any gene or common gene settings
    set<string> seen;
    doExonMetrics = true;
    //toggle the following line to switch between all blocks or any blocks mode
    for (auto run = genes.begin(); run != genes.end(); ++run) if(!run->size()) doExonMetrics = false;
    if (doExonMetrics)
    {
        doExonMetrics = false;
        for (auto run = genes.begin(); run != genes.end(); ++run)
        {
            for (auto gene = run->begin(); gene != run->end(); ++gene)
            {
                if (!isIn(seen, *gene))
                {
                    //probably just easier to record the intragenic state of each alignment block
                    if (exonCoverageCollector.queryGene(*gene)) geneCoverage[*gene]++;
                    //cout << alignment.Name << " -> " << *gene << endl;
                    exonCoverageCollector.collect(*gene);
                    //cout << endl;
                    doExonMetrics = true;
                    seen.insert(*gene);
                }
            }
        }
    }
    
    /*/
    
    if (genes.size() > 1)
    {
        //if there was more than one block, iterate through each block's set of genes and intersect them
        //In the end, we only care about genes that are common to each block
        //In theory, there's only one gene per block (in most cases) but I won't limit us on that assumption
        set<string> last = genes.front();
        for (int i = 1; i < genes.size(); ++i)
        {
            set<string> tmp;
            set_intersection(last.begin(), last.end(), genes[i].begin(), genes[i].end(), inserter(tmp, tmp.begin()));
            last = tmp;
        }
        //after the intersection, iterate over the remaining genes and record their coverage
        for (auto gene = last.begin(); gene != last.end(); ++gene)
        {
            if (exonCoverageCollector.queryGene(*gene)) geneCoverage[*gene]++;
//            cout << alignment.Name << endl;
            exonCoverageCollector.collect(*gene); //collect and keep exon coverage for this gene
            doExonMetrics = true;
        }
    }
    //if there was only one block, we just do the same coverage recording as above
    else if (genes.size() == 1) for (auto gene = genes[0].begin(); gene != genes[0].end(); ++gene)
    {
        if (exonCoverageCollector.queryGene(*gene)) geneCoverage[*gene]++;
//        cout << alignment.Name << endl;
        exonCoverageCollector.collect(*gene);
        doExonMetrics = true;
    }/**/
    if (!exonCoverageCollector.isDirty()) //a.k.a: No exons were detected at all on any block of the read
    {
        if (intragenic)
        {
            counter.increment("Intronic");
            counter.increment("Intragenic Reads");
        }
        else counter.increment("Intergenic Reads");
    }
    else if (doExonMetrics) //if exons were detected and at least one exon ended up being collected, we count this as exonic
    {
        counter.increment("Exonic");
        counter.increment("Intragenic Reads");
        if (split) counter.increment("Split Reads");
    }
    else
    {
        //It's unclear how to properly classify these reads
        //They had exon coverage, but aligned to multiple genes
        //Any exon and gene coverage they had was discarded and not recorded
        counter.increment("Intron/Exon disqualified reads");
    }
    if (ribosomal) counter.increment("rRNA");
    //also record strandedness counts
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
    unsigned short chr = chromosomeMap(chrName); //generate the chromosome shorthand referemce
    bool firstBlock = true, sameExon = true; //for keeping track of the alignment state
    string exonName = ""; // the name of the intersected exon from the bed
    
    trimFeatures(alignment, (*bedFeatures)[chr]); //trim out the features to speed up intersections
    for (auto block = blocks.begin(); sameExon && block != blocks.end(); ++block)
    {
        //for each block, intersect it with the bed file features
        list<Feature> *results = intersectBlock(*block, (*bedFeatures)[chr]);
        if (results->size() == 1) //if the block intersected more than one exon, it's immediately disqualified
        {
            if (firstBlock) exonName = results->begin()->feature_id; //record the exon name on the first pass
            else if (exonName != results->begin()->feature_id) //ensure the same exon name on subsequent passes
            {
                sameExon = false;
                break;
            }
        }
        else sameExon = false;
        delete results; //clean up dynamic allocation
        firstBlock = false;
    }
    if (sameExon && exonName.size()) //if all blocks intersected the same exon, take a fragment size sample
    {
        //both mates in a pair have to intersected the same exon in order for the pair to qualify for the sample
        auto fragment = fragments.find(alignment.Name);
        if (fragment == fragments.end()) //first time we've encountered a read in this pair
        {
            fragments[alignment.Name] = exonName;
        }
        else if (exonName == fragments[alignment.Name]) //second time we've encountered a read in this pair
        {
            //This pair is useable for fragment statistics:  both pairs fully aligned to the same exon
            fragmentSizes.push_back(abs(alignment.InsertSize));
            fragments.erase(fragment);
            --doFragmentSize;
            if (!doFragmentSize)
            {
                cout << "Completed taking fragment size samples" << endl;
                delete bedFeatures; //after taking all the samples we need, clean up the dynamic allocation
            }
        }
    }
    //return the remaining count of fragment samples to take
    return doFragmentSize;
}
