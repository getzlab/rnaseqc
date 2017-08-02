// IntervalTree.cpp : Defines the entry point for the console application.

//Include headers
#include "BED.h"
#include "Metrics.h"
#include "Expression.h"
#include <string>
#include <iostream>
#include <exception>
#include <iterator>
#include <stdio.h>
#include <set>
#include <regex>
#include <list>
#include <ctime>
#include <math.h>
#include <bamtools/api/BamReader.h>
#include <bamtools/api/BamAlignment.h>
#include <args.hxx>
#include <boost/filesystem.hpp>
using namespace std;
using namespace args;
using namespace BamTools;

const string NM = "NM";
bool debugging = false;

int main(int argc, char* argv[])
{
    //Set up command line syntax
    ArgumentParser parser("rnaSeQC");
    HelpFlag help(parser, "help", "Display this message and quit", {'h', "help"});
    Positional<string> gtfFile(parser, "gtf", "The input GTF file containing features to check the bam against");
    Positional<string> bamFile(parser, "bam", "The input SAM/BAM file containing reads to process");
    ValueFlag<string> bedFile(parser, "BEDFILE", "Optional input BED file containing non-overlapping exons used for fragment size calculations", {"bed"});
    ValueFlag<int> chimericDistance(parser, "DISTANCE", "Set the maximum accepted distance between read mates.  Mates beyond this distance will be counted as chimeric pairs. Default: 2Mb", {"chimeric-distance"});
    ValueFlag<unsigned int> maxReadLength(parser, "LENGTH", "Set the maximum accepted length.  Reads longer than this threshold are discarded. Default: 100Kb", {"read-length"});
    ValueFlag<unsigned int> fragmentSamples(parser, "SAMPLES", "Set the number of samples to take when computing fragment sizes.  Requires the --bed argument. Default: 1,000,000", {"samples"});
    ValueFlag<unsigned int> lowQualityThreshold(parser, "QUALITY", "Set the lower bound on read quality. Reads below this number are counted as low quality BUT ARE STILL USED IN COUNTS. See --exon-quality to discard reads based on quality. Default: 255", {"low-quality"});
    ValueFlag<unsigned int> exonQualityThreshold(parser,"QUALITY", "Set the lower bound on read quality for exon coverage counting. Reads below this number are excluded from coverage metrics. Default: 255", {"exon-quality"});
    ValueFlag<unsigned int> exonMismatchThreshold(parser, "MISMATCHES", "Set the maximum number of allowed mismatches between a read and the reference sequence. Reads with more than this number of mismatches are excluded from coverage metrics. Default: 6", {"exon-mismatch"});
    ValueFlag<int> splitDistance(parser, "DISTANCE", "Set the maximum distance between aligned blocks of a read.  Reads with aligned blocks separated by more than this distance are counted as split reads, BUT ARE STILL USED IN COUNTS. Default: 100bp", {"split-distance"});
    Flag debugMode(parser, "debug", "Include values of various internal constants in the output", {'d', "debug"});
    Flag LegacyMode(parser, "legacy", "Use legacy gene counting rules.  Gene counts match output of RNA-SeQC 1.1.6", {"legacy"});
    CounterFlag verbosity(parser, "verbose", "Give some feedback about what's going on.  Supply this argument twice for progress updates while parsing the bam", {'v', "verbose"});
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
        const unsigned int EXON_MISMATCH_THRESHOLD = exonMismatchThreshold ? exonMismatchThreshold.Get() : 6u;
        const unsigned int EXON_QUALITY_THRESHOLD = exonQualityThreshold ? exonQualityThreshold.Get() : 255u;
        const int SPLIT_DISTANCE = 100;
        const int VERBOSITY = verbosity ? verbosity.Get() : 0;
        debugging = debugMode.Get();
        
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
            
            if (VERBOSITY) cout<<"Reading GTF Features..."<<endl;
            time(&t0);
            while ((reader >> line))
            {
                if(LegacyMode.Get() && line.end == line.start)
                {
                    //legacy code excludes single base exons
                    if (VERBOSITY > 1) cout<<"Legacy excluded exon: " << line.feature_id << endl;
                    continue;
                }
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
        if (VERBOSITY) cout << "Finished processing GTF in " << difftime(t1, t0) << " seconds" << endl;
        
        if (bedFile) //If we were given a BED file, parse it for fragment size calculations
        {
             Feature line; //current feature being read from the bed
            if (VERBOSITY) cout << "Parsing BED intervals for fragment size computations..." << endl;
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
            if (VERBOSITY) cout<<"Parsing bam..."<<endl;
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
                    if (VERBOSITY > 1) cout << "Time elapsed: " << difftime(t2, t1) << "; Alignments processed: " << alignmentCount << endl;
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
                                if (alignment.IsDuplicate())counter.increment("Duplicate Fragments");
                                else counter.increment("Unique Fragments");
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
                        if (alignment.IsDuplicate())counter.increment("Mapped Duplicate Reads");
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
                                if (VERBOSITY > 1) cout << "Chimeric pair: " << (sequences.Begin()+alignment.RefID)->Name << "\t" << (sequences.Begin()+alignment.MateRefID)->Name << ":\t" << (alignment.RefID == alignment.MateRefID) << endl;
                                continue;
                            }
                        }
                        
                        //now record intron/exon metrics by intersecting filtered reads with the list of features
                        if (alignment.RefID < 0 || alignment.RefID >= sequences.Size())
                        {
                            //The read had an unrecognized RefID (one not defined in the bam's header)
                            if (VERBOSITY) cout << "Unrecognized RefID on alignment: " << alignment.Name<<endl;
                        }
                        else if(mismatches <= EXON_MISMATCH_THRESHOLD && alignment.IsProperPair() && alignment.MapQuality >= EXON_QUALITY_THRESHOLD)
                        {
                            vector<Feature> blocks;
                            string chrName = (sequences.Begin()+alignment.RefID)->Name;
                            unsigned short chr = chromosomeMap(chrName); //parse out a chromosome shorthand
                            //extract each cigar block from the alignment
                            unsigned int length = extractBlocks(alignment, blocks, chr);
                            trimFeatures(alignment, features[chr]); //drop features that appear before this read
                            
                            //run the read through exon metrics
                            //exonAlignmentMetrics(SPLIT_DISTANCE, features, counter, sequences, geneCoverage, exonCoverage, blocks, alignment, length);
                            if (LegacyMode.Get()) legacyExonAlignmentMetrics(SPLIT_DISTANCE, features, counter, sequences, geneCoverage, exonCoverage, blocks, alignment, length);
                            else exonAlignmentMetrics(SPLIT_DISTANCE, features, counter, sequences, geneCoverage, exonCoverage, blocks, alignment, length);
                            
                            //if fragment size calculations were requested, we still have samples to take, and the chromosome exists within the provided bed
                            if (alignment.IsPaired() && doFragmentSize && bedFeatures != nullptr && bedFeatures->find(chr) != bedFeatures->end())
                            {
                                doFragmentSize = fragmentSizeMetrics(doFragmentSize, bedFeatures, fragments, fragmentSizes, sequences, blocks, alignment);
                                if (VERBOSITY > 1 && !doFragmentSize) cout << "Completed taking fragment size samples" << endl;
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
        if (VERBOSITY)
        {
            cout<< "Time Elapsed: " << difftime(t2, t1) << "; Alignments processed: " << alignmentCount << endl;
            cout << "Total runtime: " << difftime(t2, t0) << "; Total CPU Time: " << (clock() - start_clock)/CLOCKS_PER_SEC << endl;
            cout << "Estimating library complexity..." << endl;
        }
        double duplicates = (double) counter.get("Duplicate Fragments");
        double unique = (double) counter.get("Unique Fragments");
        double numReads = duplicates + unique;
        unsigned int minReads = 0u, minError = UINT_MAX;
        if (duplicates > 0)
        {
            //If there are no duplicates, the estimate is useless, so skip it
            for (double x = unique; x < 1e9; ++x)
            {
                double estimate = x * (1.0 - exp(-1.0 * numReads / x)); //lander-waterman
                unsigned int error = (unsigned int) fabs(estimate - unique);
                if (error < minError)
                {
                    minError = error;
                    minReads = (unsigned int) x;
                }
            }
        }

        if (VERBOSITY) cout << "Generating report" << endl;
        
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
        output << "Disqualification Rate\t" << counter.frac("Intron/Exon disqualified reads", "Mapped Reads") << endl;
        output << "rRNA Rate\t" << counter.frac("rRNA", "Mapped Reads") << endl;
        output << "End 1 Sense Rate\t" << (double) counter.get("End 1 Sense") / (counter.get("End 1 Sense") + counter.get("End 1 Antisense")) << endl;
        output << "End 2 Sense Rate\t" << (double) counter.get("End 2 Sense") / (counter.get("End 2 Sense") + counter.get("End 2 Antisense")) << endl;
        output << "Avg. Blocks per Read\t" << counter.frac("Alignment Blocks", "Mapped Reads") << endl;
        //automatically dump the raw counts of all metrics to the file
        output << counter;
        //append metrics that were manually tracked
        output << "Read Length\t" << readLength << endl;
        output << "Genes Detected\t" << genesDetected << endl;
        output << "Estimated Library Complexity\t" << minReads << endl;
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
