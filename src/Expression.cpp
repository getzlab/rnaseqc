//
//  Expression.cpp
//  IntervalTree
//
//  Created by Aaron Graubert on 8/2/17.
//  Copyright © 2017 Aaron Graubert. All rights reserved.
//

#include "Expression.h"
#include <algorithm>

using namespace BamTools;
using std::vector;
using std::list;
using std::map;
using std::string;
using std::set;
using std::cout;
using std::endl;
using std::pair;


//this actually is the legacy version, but it works out the same and makes alignment size math a little easier
unsigned int extractBlocks(BamAlignment &alignment, vector<Feature> &blocks, chrom chr, bool legacy)
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
                block.strand = alignment.IsReverseStrand() ? -1 : 1;
                blocks.push_back(block);
                alignedSize += current.Length;
            case 'N':
            case 'D':
                //M, =, X, N, and D blocks all advance the start position of the next block
                start += current.Length;
            case 'H':
            case 'P':
//            case 'S':
            case 'I':
                break;
            case 'S':
                if (legacy) alignedSize += current.Length;
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
    while (!features.empty() && features.front().end < alignment.Position) features.pop_front();
}

void trimFeatures(BamAlignment &alignment, list<Feature> &features, BaseCoverage &coverage)
{
    //trim intervals upstream of this block
    //Since alignments are sorted, if an alignment occurs beyond any features, these features can be dropped
    while (!features.empty() && features.front().end < alignment.Position)
    {
        if (features.front().type == "gene") coverage.compute(features.front()); //Once this gene leaves the search window, compute coverage
        features.pop_front();
    }
}

void dropFeatures(std::list<Feature> &features, BaseCoverage &coverage)
{
    for (auto feat = features.begin(); feat != features.end(); ++feat) if (feat->type == "gene") coverage.compute(*feat);
    features.clear();
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

void legacyExonAlignmentMetrics(unsigned int SPLIT_DISTANCE, map<chrom, list<Feature>> &features, Metrics &counter, SamSequenceDictionary &sequences, vector<Feature> &blocks, BamAlignment &alignment, unsigned int length, unsigned short stranded, BaseCoverage &baseCoverage)
{
    string chrName = (sequences.Begin()+alignment.RefID)->Name;
    chrom chr = chromosomeMap(chrName); //generate the chromosome shorthand name
    const bool dbg = false;
    if (dbg) cout << "~" << alignment.Name;

    //check for split reads by iterating over all the blocks of this read
    bool split = false;
    long long lastEnd = -1; // used for split read detection
    for(auto block = blocks.begin(); block != blocks.end(); ++block)
    {
        if (lastEnd > 0 && !split) split = (block->start - lastEnd) > SPLIT_DISTANCE - 1;
        lastEnd = block->end;
    }


    counter.increment("Alignment Blocks", blocks.size());
    counter.increment("Reads used for Intron/Exon counts");
    //Bamtools uses 0-based indexing because it's the 'norm' in computer science, even though sams are 1-based
    Feature current; //current block of the alignment (used while iterating)
    current.start = alignment.Position+1; //0-based + 1 == 1-based
    current.end = alignment.GetEndPosition(); //0-based, open == 1-based, closed

    list<Feature> *results = intersectBlock(current, features[chr]);

    vector<set<string> > genes; //each set is the set of genes intersected by the current block (one set per block)
    bool intragenic = false, transcriptPlus = false, transcriptMinus = false, ribosomal = false, doExonMetrics = false, exonic = false; //various booleans for keeping track of the alignment
    bool legacyNotSplit = false; //Legacy bug to override a read being split
    for (auto result = results->begin(); result != results->end(); ++result)
    {
        Feature exon;
        bool legacyFoundExon = false, legacyFoundGene = false;
        map<string, float> legacySplitDosage;
        legacyNotSplit = false;
        if (result->type == "gene")
        {
            intragenic = true;
            if (result->strand == 1) transcriptPlus = true;
            else if (result->strand == -1) transcriptMinus = true;
            for (auto block = blocks.begin(); block != blocks.end(); ++block)
            {
                //if (block->start == block->end) continue;
                if (stranded)
                {
                    bool target = alignment.IsReverseStrand() ^ alignment.IsFirstMate();
                    if (stranded == 1) target = !target;
                    if (result->strand != (target ? 1 : -1))
                    {
                        continue;
                    }
                }

//                bool legacyMatchIntron = false;
                bool firstexon = false;
                legacyFoundExon = false;
//                if (intersectInterval(*result, *block))
                {
                    legacyFoundGene = true;
//                    bias.checkBias(*result, *block); //bias checking now handled in baseCoverage
                    for (auto ex = results->begin(); ex != results->end() && !firstexon ; ++ex)
                    {
                        if (ex->type == "exon" && ex->gene_id == result->gene_id && intersectInterval(*ex, *block)  )
                        {
                            if (result->ribosomal) ribosomal = true;
                            if (partialIntersect(*ex, *block) == (block->end - block->start))
                            {
                                exon = *ex;
                                exonic = true;
                                firstexon = true;
                                legacyFoundExon=true; //should this be part of the loop condition?  look into overlapsIntronp
                                baseCoverage.add(*ex, block->start, block->end);
                            }
                        }

                    }
                    if (split && !legacyNotSplit)
                    {
                        if (legacyFoundExon)
                        {
                            legacySplitDosage[exon.feature_id] += (float) (block->end - block->start) / (float) alignment.Length;//length;
                        }
                        else legacyNotSplit = true;
                    }
                }
            }
            //record gene to collect, probs want to keep a list of single and partial collections
            if (legacyFoundExon)
            {
                if (split && !legacyNotSplit)
                {
                    for (auto coverage = legacySplitDosage.begin(); coverage != legacySplitDosage.end(); ++coverage)
                    {
                        exonCounts[coverage->first] += coverage->second;
                        if (dbg) cout << "\t" << coverage->first << " " << coverage->second;
                    }
                }
                else
                {
                    //If read was not detected as split or the legacy bug changed it to unsplit, only record last exon
                    exonCounts[exon.feature_id] += 1.0;
                    if (dbg) cout << "\t" << exon.feature_id << " 1.0";
                }
                geneCounts[exon.gene_id] += 1.0;
                baseCoverage.commit(exon.gene_id);
                doExonMetrics = true;
            }
        }
        //if (split && !legacyNotSplit && !legacyFoundGene) legacyNotSplit = true;

    }
    delete results;
    if (dbg) cout << endl;

    if (!exonic) //a.k.a: No exons were detected at all on any block of the read
    {
        if (intragenic)
        {
            counter.increment("Intronic Reads");
            counter.increment("Intragenic Reads");
        }
        else counter.increment("Intergenic Reads");
    }
    else if (doExonMetrics) //if exons were detected and at least one exon ended up being collected, we count this as exonic
    {
        counter.increment("Exonic Reads");
        counter.increment("Intragenic Reads");
        if (split && !legacyNotSplit) counter.increment("Split Reads");
    }
    else
    {
        //It's unclear how to properly classify these reads
        //They had exon coverage, but aligned to multiple genes
        //Any exon and gene coverage they had was discarded and not recorded
        counter.increment("Intron/Exon Disqualified Reads");
    }
    if (ribosomal) counter.increment("rRNA Reads");
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
    baseCoverage.reset();
}


void exonAlignmentMetrics(unsigned int SPLIT_DISTANCE, map<chrom, list<Feature>> &features, Metrics &counter, SamSequenceDictionary &sequences, vector<Feature> &blocks, BamAlignment &alignment, unsigned int length, unsigned short stranded, BaseCoverage &baseCoverage)
{
    string chrName = (sequences.Begin()+alignment.RefID)->Name;
    chrom chr = chromosomeMap(chrName); //generate the chromosome shorthand name

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
    Collector exonCoverageCollector(&exonCounts); //Collects coverage counts for later (counts may be discarded)
    bool intragenic = false, transcriptPlus = false, transcriptMinus = false, ribosomal = false, doExonMetrics = false, exonic = false; //various booleans for keeping track of the alignment

    for (auto block = blocks.begin(); block != blocks.end(); ++block)
    {
        genes.push_back(set<string>()); //create a new set for this block
        list<Feature> *results = intersectBlock(*block, features[chr]); //grab the list of intersecting features
        for (auto result = results->begin(); result != results->end(); ++result)
        {
            if (stranded)
            {
                bool target = alignment.IsReverseStrand() ^ alignment.IsFirstMate();
                if (stranded == 1) target = !target;
                if (result->strand != (target ? 1 : -1))
                {

                    continue;
                }
            }
            //if (stranded && result->strand != (alignment.IsReverseStrand() ^ (stranded == 1 ? alignment.IsFirstMate() : alignment.IsSecondMate()) ? -1 : 1)) continue;
            if (result->strand == 1) transcriptPlus = true;
            else if (result->strand == -1) transcriptMinus = true;
            //else...what, exactly?
            if (result->type == "exon")
            {
                exonic = true;
                int intersectionSize = partialIntersect(*result, *block);
                //check that this block fully overlaps the feature
                //(if any bases of the block don't overlap, then the read is discarded)
                if (intersectionSize == block->end - block->start)
                {
                    //store the exon split dosage coverage in the collector for now
                    genes.rbegin()->insert(result->gene_id);
                    double tmp = static_cast<double>(intersectionSize) / length;
                    exonCoverageCollector.add(result->gene_id, result->feature_id, tmp);
                    baseCoverage.add(*result, block->start, block->end); //provisionally add per-base coverage to this gene

                }

            }
            else if (result->type == "gene")
            {
                intragenic = true;
                //we don't record the gene name here because in terms of gene coverage and detection, we only care about exons

                //bias checking now handled in baseCoverage
//                bias.checkBias(*result, *block);
            }
            if (result->ribosomal) ribosomal = true;
        }
        delete results; //clean up dynamic allocation
    }

    if (genes.size() >= 1)
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
            if (exonCoverageCollector.queryGene(*gene)) geneCounts[*gene]++;
            exonCoverageCollector.collect(*gene); //collect and keep exon coverage for this gene
            baseCoverage.commit(*gene); //keep the per-base coverage recorded on this gene
            doExonMetrics = true;
        }
    }

    if (!exonic) //a.k.a: No exons were detected at all on any block of the read
    {
        if (intragenic)
        {
            counter.increment("Intronic Reads");
            counter.increment("Intragenic Reads");
        }
        else counter.increment("Intergenic Reads");
    }
    else if (doExonMetrics) //if exons were detected and at least one exon ended up being collected, we count this as exonic
    {
        counter.increment("Exonic Reads");
        counter.increment("Intragenic Reads");
        if (split) counter.increment("Split Reads");
    }
    else
    {
        //It's unclear how to properly classify these reads
        //They had exon coverage, but aligned to multiple genes
        //Any exon and gene coverage they had was discarded and not recorded
        counter.increment("Intron/Exon Disqualified Reads");
    }
    if (ribosomal) counter.increment("rRNA Reads");
    //also record strandedness counts
    //TODO: check standing metrics.  Counts are probably off because of null intron/exon calls
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
    baseCoverage.reset();
}

unsigned int fragmentSizeMetrics(unsigned int doFragmentSize, map<chrom, list<Feature>> *bedFeatures, map<string, string> &fragments, map<long long, unsigned long> &fragmentSizes, SamSequenceDictionary &sequences, vector<Feature> &blocks, BamAlignment &alignment)
{
    string chrName = (sequences.Begin()+alignment.RefID)->Name;
    chrom chr = chromosomeMap(chrName); //generate the chromosome shorthand referemce
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
                delete results;
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
            fragmentSizes[abs(alignment.InsertSize)] += 1;
            fragments.erase(fragment);
            --doFragmentSize;
            if (!doFragmentSize)
            {
                delete bedFeatures; //after taking all the samples we need, clean up the dynamic allocation
                bedFeatures = nullptr;
            }
        }
    }
    //return the remaining count of fragment samples to take
    return doFragmentSize;
}


std::string buildSequence(BamTools::BamAlignment &alignment)
{
    if (alignment.HasTag("MD"))
    {
        std::string readSeq = alignment.AlignedBases, mdSeq;
        alignment.GetTag("MD", mdSeq);
        //Let's assume the CIGAR string matches the MD tag (for now)
        
        
    }
    return "";
}

/*unsigned int extractBlocks(BamAlignment &alignment, vector<Feature> &blocks, chrom chr)
{
    unsigned int alignedSize = 0;
    auto beg = alignment.CigarData.begin();
    auto end = alignment.CigarData.end();
    coord start = alignment.Position + 1;
    Feature block;
    bool filled = false;
    while (beg != end)
    {
        CigarOp current = *(beg++);
        switch(current.Type)
        {
            case 'M':
            case '=':
            case 'X':
                //M, =, and X blocks are aligned, so push back this block
                if (!filled) block.start = start;
                block.chromosome = chr;
                block.end = start + current.Length; //1-based, closed
                //blocks.push_back(block);
                alignedSize += current.Length;
                filled = true;
                start += current.Length;
                break;
            case 'N':
                //N blocks push back the previous M block
                //The goal of splitting this up is that MNM strings produce 2 blocks, but MDM would produce one long block
                if (filled)
                {
                    blocks.push_back(block);
                    filled = false;
                }
            case 'D':
                //D blocks only advance the start position of the next block
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
    if (filled)
    {
        blocks.push_back(block);
    }
    return alignedSize;
}*/
