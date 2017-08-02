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

//this actually is the legacy version, but it works out the same and makes alignment size math a little easier
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
                
                if (intersectionSize == block->end - block->start)
                {
                    genes.rbegin()->insert(result->gene_id);
                    blockWasIntragenic = true; //legacy
                    double tmp = (double) intersectionSize / length;
                    exonCoverageCollector.add(result->gene_id, result->feature_id, tmp);
                }
                
            }
            else if (result->type == "gene")
            {
                intragenic = true;
                //we don't record the gene name here because in terms of gene coverage and detection, we only care about exons (apparently not)
            }
            if (result->transcript_type == "rRNA" || result->transcript_type == "Mt_rRNA") ribosomal = true;
        }
        delete results; //clean up dynamic allocation
        //legacy split read handling:
        if (split && !blockWasIntragenic && genes.size() == 1 && blocks.size() > 1)
        {
            genes.pop_back();
            split = false;
        }
    }
    //legacy code only examines 'final' (closer to the end of the contig) segment of each read
    for (auto gene = genes.rbegin()->begin(); gene != genes.rbegin()->end(); ++gene)
    {
        if (exonCoverageCollector.queryGene(*gene))
        {
            geneCoverage[*gene]++;
        }
        if(split) exonCoverageCollector.collect(*gene);
        else exonCoverageCollector.collectSingle(*gene);
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
    
    for (auto block = blocks.begin(); block != blocks.end(); ++block)
    {
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
                if (intersectionSize == block->end - block->start)
                {
                    //store the exon split dosage coverage in the collector for now
                    genes.rbegin()->insert(result->gene_id);
                    double tmp = (double) intersectionSize / length;
                    exonCoverageCollector.add(result->gene_id, result->feature_id, tmp);
                    
                }
                
            }
            else if (result->type == "gene")
            {
                intragenic = true;
                //we don't record the gene name here because in terms of gene coverage and detection, we only care about exons (apparently not)
            }
            if (result->transcript_type == "rRNA" || result->transcript_type == "Mt_rRNA") ribosomal = true;
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
            if (exonCoverageCollector.queryGene(*gene)) geneCoverage[*gene]++;
            exonCoverageCollector.collect(*gene); //collect and keep exon coverage for this gene
            doExonMetrics = true;
        }
    }
    //if there was only one block, we just do the same coverage recording as above
    /*else if (genes.size() == 1) for (auto gene = genes[0].begin(); gene != genes[0].end(); ++gene)
     {
     if (exonCoverageCollector.queryGene(*gene)) geneCoverage[*gene]++;
     exonCoverageCollector.collect(*gene);
     doExonMetrics = true;
     }*/
    
    /*if (exonCoverageCollector.sum() > 1.0)
    {
        cout << alignment.Name << " had excess exon coverage: " << exonCoverageCollector.sum() << endl;
    }*/
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
                delete bedFeatures; //after taking all the samples we need, clean up the dynamic allocation
                bedFeatures = nullptr;
            }
        }
    }
    //return the remaining count of fragment samples to take
    return doFragmentSize;
}

/*unsigned int extractBlocks(BamAlignment &alignment, vector<Feature> &blocks, unsigned short chr)
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
