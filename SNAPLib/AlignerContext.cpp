/*++


Module Name:

    AlignerContext.cpp

Abstract:

    Common parameters for running single & paired alignment.

Authors:

    Ravi Pandya, May, 2012

Environment:

    User mode service.

Revision History:

    Integrated from SingleAligner.cpp & PairedAligner.cpp

--*/

#include "stdafx.h"
#include "Compat.h"
#include "options.h"
#include "AlignerOptions.h"
#include "AlignerContext.h"
#include "AlignerStats.h"
#include "BaseAligner.h"
#include "FileFormat.h"
#include "exit.h"
#include "PairedAligner.h"
#include "Error.h"
#include "Util.h"
#include "CommandProcessor.h"
#include "T2TrRNA.h"

using std::max;
using std::min;

#if INSTRUMENTATION_FOR_PAPER
_int64 g_alignmentTimeByHitCountsOfEachSeed[MAX_HIT_SIZE_LOG_2 + 1][MAX_HIT_SIZE_LOG_2 + 1];  // In the paired-end aligner, if you have seeds A and B with hit set sizes |A| and |B| then the total time in ns gets added into g_alignmentTimeByHitCountsOfEachSeed[log2(|A|)][log2(|B|)]
_int64 g_alignmentCountByHitCountsOfEachSeed[MAX_HIT_SIZE_LOG_2 + 1][MAX_HIT_SIZE_LOG_2 + 1];  // Same as above, but just add one per time.
_int64 g_scoreCountByHitCountsOfEachSeed[MAX_HIT_SIZE_LOG_2 + 1][MAX_HIT_SIZE_LOG_2 + 1];
_int64 g_setIntersectionSizeByHitCountsOfEachSeed[MAX_HIT_SIZE_LOG_2 + 1][MAX_HIT_SIZE_LOG_2 + 1];
_int64 g_100xtotalRatioOfSetIntersectionSizeToSmallerSeedHitCountByCountsOfEachSeed[MAX_HIT_SIZE_LOG_2 + 1][MAX_HIT_SIZE_LOG_2 + 1];
_int64 g_totalSizeOfSmallerHitSet = 0;
_int64 g_totalSizeOfSetIntersection = 0;
_int64 g_alignmentsWithMoreThanOneCandidateWhereTheBestCandidateIsScoredFirst = 0;
_int64 g_alignmentsWithMoreThanOneCandidate = 0;
#endif // INSTRUMENTATION_FOR_PAPER

//
// Save the index & index directory globally so that we don't need to reload them on multiple runs.
//
GenomeIndex *g_index = NULL;
char *g_indexDirectory = NULL;

AlignerContext::AlignerContext(int i_argc, const char **i_argv, const char *i_version, AlignerExtension* i_extension)
    :
    index(NULL),
    writerSupplier(NULL),
    options(NULL),
    stats(NULL),
    extension(i_extension != NULL ? i_extension : new AlignerExtension()),
    readWriter(NULL),
    argc(i_argc),
    argv(i_argv),
    version(i_version),
    perfFile(NULL),
    statFile(NULL)
{
}

AlignerContext::~AlignerContext()
{
    delete extension;
    if (NULL != perfFile)
        fclose(perfFile);
    if (NULL != statFile)
        fclose(statFile);
    delete stats;
}


void AlignerContext::runAlignment(int argc, const char **argv, const char *version, unsigned *argsConsumed)
{
    options = parseOptions(argc, argv, version, argsConsumed, isPaired());

    if (NULL == options) { // Didn't parse correctly
        *argsConsumed = argc;
        return;
    }

#ifdef _MSC_VER
    useTimingBarrier = options->useTimingBarrier;
#endif

#if INSTRUMENTATION_FOR_PAPER
    for (int i = 0; i < MAX_HIT_SIZE_LOG_2 + 1; i++) {
        for (int j = 0; j < MAX_HIT_SIZE_LOG_2 + 1; j++) {
            g_alignmentTimeByHitCountsOfEachSeed[i][j] = 0;
            g_alignmentCountByHitCountsOfEachSeed[i][j] = 0;
            g_scoreCountByHitCountsOfEachSeed[i][j] = 0;
            g_setIntersectionSizeByHitCountsOfEachSeed[i][j] = 0;
            g_100xtotalRatioOfSetIntersectionSizeToSmallerSeedHitCountByCountsOfEachSeed[i][j] = 0;
        } // j
    } // i
#endif // INSTRUMENTATION_FOR_PAPER


    if (!initialize()) {
        return;
    }
    extension->initialize();

    if (! extension->skipAlignment()) {
        WriteStatusMessage("Aligning.\n");

        beginIteration();

        runTask();

        finishIteration();

        printStats();

        nextIteration();    // This probably should get rolled into something else; it's really cleanup code, not "next iteration"
    }

    extension->finishAlignment();
    PrintBigAllocProfile();
    PrintWaitProfile();

#if INSTRUMENTATION_FOR_PAPER
    FILE* outputFile = fopen("SNAPInstrumentation.txt", "w");

    if (outputFile == NULL) {
        fprintf(stderr, "Unable to open instrumentation output file\n");
    } else {
        fprintf(outputFile, "%lld alignments have more than one candidate scored.  Of those, %lld had the best candidate scored first\n", g_alignmentsWithMoreThanOneCandidate, g_alignmentsWithMoreThanOneCandidateWhereTheBestCandidateIsScoredFirst);
        fprintf(outputFile, "Total size of set intersection %lld.  Total size of smaller set %lld\n", g_totalSizeOfSetIntersection, g_totalSizeOfSmallerHitSet);

        fprintf(outputFile, "Alignment count by hit counts of each seed\nHits 1/2");
        for (int i = 0; i < MAX_HIT_SIZE_LOG_2 + 1; i++) {
            fprintf(outputFile, "\t%d", 1 << i);
        }
        fprintf(outputFile, "\n");

        for (int i = 0; i < MAX_HIT_SIZE_LOG_2 + 1; i++) {
            fprintf(outputFile, "%d", 1 << i);
            for (int j = 0; j < MAX_HIT_SIZE_LOG_2 + 1; j++) {
                fprintf(outputFile, "\t%lld", g_alignmentCountByHitCountsOfEachSeed[i][j]);
            } // j
            fprintf(outputFile, "\n");
        } // i

        fprintf(outputFile, "\nAlignment time by hit counts of each seed (ns)\nHits 1/2");
        for (int i = 0; i < MAX_HIT_SIZE_LOG_2 + 1; i++) {
            fprintf(outputFile, "\t%d", 1 << i);
        }
        fprintf(outputFile, "\n");

        for (int i = 0; i < MAX_HIT_SIZE_LOG_2 + 1; i++) {
            fprintf(outputFile, "%d", 1 << i);
            for (int j = 0; j < MAX_HIT_SIZE_LOG_2 + 1; j++) {
                fprintf(outputFile, "\t%lld", g_alignmentTimeByHitCountsOfEachSeed[i][j]);
            } // j
            fprintf(outputFile, "\n");
        } // i

        fprintf(outputFile, "\nLocations scored by hit counts of each seed\nHits 1/2");
        for (int i = 0; i < MAX_HIT_SIZE_LOG_2 + 1; i++) {
            fprintf(outputFile, "\t%d", 1 << i);
        }
        fprintf(outputFile, "\n");

        for (int i = 0; i < MAX_HIT_SIZE_LOG_2 + 1; i++) {
            fprintf(outputFile, "%d", 1 << i);
            for (int j = 0; j < MAX_HIT_SIZE_LOG_2 + 1; j++) {
                fprintf(outputFile, "\t%lld", g_scoreCountByHitCountsOfEachSeed[i][j]);
            } // j
            fprintf(outputFile, "\n");
        } // i

        fprintf(outputFile, "\nSet intersection size by hit counts of each seed\nHits 1/2");
        for (int i = 0; i < MAX_HIT_SIZE_LOG_2 + 1; i++) {
            fprintf(outputFile, "\t%d", 1 << i);
        }
        fprintf(outputFile, "\n");

        for (int i = 0; i < MAX_HIT_SIZE_LOG_2 + 1; i++) {
            fprintf(outputFile, "%d", 1 << i);
            for (int j = 0; j < MAX_HIT_SIZE_LOG_2 + 1; j++) {
                fprintf(outputFile, "\t%lld", g_setIntersectionSizeByHitCountsOfEachSeed[i][j]);
            } // j
            fprintf(outputFile, "\n");
        } // i

        fprintf(outputFile, "\n100x Total of Ratio of Set Intersection size to size of smaller seed hit count\nHits 1/2");
        for (int i = 0; i < MAX_HIT_SIZE_LOG_2 + 1; i++) {
            fprintf(outputFile, "\t%d", 1 << i);
        }
        fprintf(outputFile, "\n");

        for (int i = 0; i < MAX_HIT_SIZE_LOG_2 + 1; i++) {
            fprintf(outputFile, "%d", 1 << i);
            for (int j = 0; j < MAX_HIT_SIZE_LOG_2 + 1; j++) {
                fprintf(outputFile, "\t%lld", g_100xtotalRatioOfSetIntersectionSizeToSmallerSeedHitCountByCountsOfEachSeed[i][j]);
            } // j
            fprintf(outputFile, "\n");
        } // i

        fclose(outputFile);
    }
#endif // INSTRUMENTATION_FOR_PAPER

}

    void
AlignerContext::initializeThread()
{
    stats = newStats(); // separate copy per thread
    stats->extra = extension->extraStats();
    readWriter = writerSupplier != NULL ? writerSupplier->getWriter() : NULL;
    extension = extension->copy();
}

    void
AlignerContext::runThread()
{
    extension->beginThread();
    runIterationThread();
    if (readWriter != NULL) {
        readWriter->close();
        delete readWriter;
    }
    extension->finishThread();
}

    void
AlignerContext::finishThread(AlignerContext* common)
{
    common->stats->add(stats);
    delete stats;
    stats = NULL;
    delete extension;
    extension = NULL;
}

    std::unordered_set<_int64>
AlignerContext::rrnaPosSet()
{
    std::unordered_set<_int64> newSet;
    for (unsigned i = 0; i < T2T_RRNA_NROW; ++i)
        for (_int64 j = T2T_RRNA_RANGE[i][0]; j <= T2T_RRNA_RANGE[i][1]; ++j)
            newSet.insert(j);
    return newSet;
}

    std::unordered_map<_int64, _int8>
AlignerContext::icPosMap()
{
    _int64 i = 0, j = 0, k = 0;
    std::unordered_map<_int64, _int8> newMap;
    int numContigs = index->getGenome()->getNumContigs();
    for (i = 0; i < numContigs; i++)
    {
        const Genome::Contig* contig = index->getGenome()->getContigByOriginalContigNumber(OriginalContigNum(i));
        const _int64 ctgLength = index->getGenome()->getContigByOriginalContigNumber(OriginalContigNum(i))->length;
		if (*contig->name == 'I' && *(contig->name + 1) == 'C')
		{
			_int64 pos = contig->beginningLocation;
			_int64 len = contig->length - index->getGenome()->getChromosomePadding();
			fprintf(stderr, "%d\t%" PRId64 "\t%" PRId64 "\n", j + 1, pos, len);
			for (k = 0; k < len; ++k)
				newMap.insert(std::make_pair(j, pos + k));
			++j;
		}
	}
    return newMap;
}

    std::unordered_set<_int64>
AlignerContext::hskPosSet()
{
    std::unordered_set<_int64> newSet;
    for (unsigned i = 0; i < T2T_HSK_NROW; ++i)
        for (_int64 j = T2T_HSK_RANGE[i][0]; j <= T2T_HSK_RANGE[i][1]; ++j)
            newSet.insert(j);
    return newSet;
}

    std::unordered_set<_int64>
AlignerContext::hskCovSet()
{
    std::unordered_set<_int64> newSet;
    return newSet;
}

    bool
AlignerContext::initialize()
{
    if (g_indexDirectory == NULL || strcmp(g_indexDirectory, options->indexDir) != 0) {
        delete g_index;
        g_index = NULL;
        delete g_indexDirectory;
        g_indexDirectory = new char [strlen(options->indexDir) + 1];
        strcpy(g_indexDirectory, options->indexDir);

        if (strcmp(options->indexDir, "-") != 0) {
            WriteStatusMessage("Loading index from directory... ");
 
            fflush(stdout);
            _int64 loadStart = timeInMillis();
            index = GenomeIndex::loadFromDirectory((char*) options->indexDir, options->mapIndex, options->prefetchIndex);
            if (index == NULL) {
                WriteErrorMessage("Index load failed, aborting.\n");
                soft_exit(1);
            }
            g_index = index;

            const int basesBufferSize = 30;
            char basesBuffer[basesBufferSize];

            _int64 loadTime = timeInMillis() - loadStart;
             WriteStatusMessage("%llds.  %s bases, seed size %d.\n",
                    loadTime / 1000, FormatUIntWithCommas(index->getGenome()->getCountOfBases(), basesBuffer, basesBufferSize), index->getSeedLength());

             if (index->getMajorVersion() < 5 || (index->getMajorVersion() == 5 && index->getMinorVersion() == 0)) {
                 WriteErrorMessage("WARNING: The version of the index you're using was built with an earlier version of SNAP and will result in Ns in the reference NOT matching Ns in reads.\n         If you do not want this behavior, rebuild the index.\n");
             }
         } else {
            WriteStatusMessage("no alignment, input/output only\n");
        }
    } else {
        index = g_index;
    }
    isT2T = index->getGenome()->getCountOfBases() >= T2T_GENOME_BASES ? true : false;
    // check chr sizes in order
    if (isT2T)
    {
        int numContigs = index->getGenome()->getNumContigs();
        if (numContigs < T2T_CHROMOSOME_NROW)
            isT2T = false;
        else
        {
            for (int i = 0; i < T2T_CHROMOSOME_NROW; ++i) // check chr order
            {
                const _int64 ctgLength = index->getGenome()->getContigByOriginalContigNumber(OriginalContigNum(i))->length;
                if (ctgLength - index->getGenome()->getChromosomePadding() != T2T_CHROMOSOME_SIZES[i])
                {
                    isT2T = false;
                    break;
                }
            }
        }
    }
    rrnapos = rrnaPosSet();
    hskpos = hskPosSet();
    icpos = icPosMap();
    hasIC = icpos.size() ? true : false;
    maxHits_ = options->maxHits;
    maxDist_ = options->maxDist;
    maxDistForIndels_ = options->maxDistForIndels;
    extraSearchDepth = options->extraSearchDepth;
    disabledOptimizations = options->disabledOptimizations;
    ignoreAlignmentAdjustmentForOm = options->ignoreAlignmentAdjustmentsForOm;
    altAwareness = options->altAwareness;
    emitALTAlignments = options->emitALTAlignments;
    maxSecondaryAlignmentAdditionalEditDistance = options->maxSecondaryAlignmentAdditionalEditDistance;
    maxSecondaryAlignments = options->maxSecondaryAlignments;
    maxSecondaryAlignmentsPerContig = options->maxSecondaryAlignmentsPerContig;
    maxScoreGapToPreferNonALTAlignment = options->maxScoreGapToPreferNonALTAlignment;
    useAffineGap = options->useAffineGap;
    matchReward = options->matchReward;
    subPenalty = options->subPenalty;
    gapOpenPenalty = options->gapOpenPenalty;
    gapExtendPenalty = options->gapExtendPenalty;
    fivePrimeEndBonus = options->fivePrimeEndBonus;
    threePrimeEndBonus = options->threePrimeEndBonus;
    useSoftClipping = options->useSoftClipping;

    if (maxSecondaryAlignmentAdditionalEditDistance < 0 && (maxSecondaryAlignments < 1000000 || maxSecondaryAlignmentsPerContig > 0)) {
        WriteErrorMessage("You set -omax and/or -mpc without setting -om.  They're meaningful only in the context of -om, so you probably didn't really mean to do that.\n");
        soft_exit(1);
    }
    minReadLength = options->minReadLength;

    if (index != NULL && (int)minReadLength < index->getSeedLength()) {
        WriteErrorMessage("The min read length (%d) must be at least the seed length (%d), or there's no hope of aligning reads that short.\n", minReadLength, index->getSeedLength());
        soft_exit(1);
    }

    if (options->perfFileName != NULL) {
        perfFile = fopen(options->perfFileName, "a");
        if (NULL == perfFile) {
            WriteErrorMessage("Unable to open perf file '%s'\n", options->perfFileName);
            soft_exit(1);
        }
    }

    if (options->statFileName != NULL) {
        statFile = fopen(options->statFileName, "w");
        if (NULL == statFile) {
            WriteErrorMessage("Unable to open stats file '%s'\n", options->statFileName);
            soft_exit(1);
        }
    }

    DataSupplier::ThreadCount = options->numThreads;

    return true;
}


    void
AlignerContext::beginIteration()
{
    writerSupplier = NULL;
    alignStart = timeInMillis();
    clipping = options->clipping;
    totalThreads = options->numThreads;
    bindToProcessors = options->bindToProcessors;
    maxDist = maxDist_;
    maxDistForIndels = maxDistForIndels_;
    maxHits = maxHits_;
    numSeedsFromCommandLine = options->numSeedsFromCommandLine;
    seedCoverage = options->seedCoverage;
    minWeightToCheck = options->minWeightToCheck;
    if (stats != NULL) {
        delete stats;
    }
    stats = newStats();
    stats->extra = extension->extraStats();
    extension->beginIteration();
    
    memset(&readerContext, 0, sizeof(readerContext));
    readerContext.clipping = options->clipping;
    readerContext.defaultReadGroup = options->defaultReadGroup;
    readerContext.genome = index != NULL ? index->getGenome() : NULL;
    readerContext.ignoreSecondaryAlignments = options->ignoreSecondaryAlignments;
    readerContext.ignoreSupplementaryAlignments = options->ignoreSecondaryAlignments;   // Maybe we should split them out
    readerContext.preserveFASTQComments = options->preserveFASTQComments;
    DataSupplier::ExpansionFactor = options->expansionFactor;

    typeSpecificBeginIteration();

    if (UnknownFileType != options->outputFile.fileType) {
        const FileFormat* format;
        if (SAMFile == options->outputFile.fileType) {
            format = FileFormat::SAM[options->useM];
        } else if (BAMFile == options->outputFile.fileType) {
            format = FileFormat::BAM[options->useM];
        } else {
            //
            // This shouldn't happen, because the command line parser should catch it.  Perhaps you've added a new output file format and just
            // forgoten to add it here.
            //
            WriteErrorMessage("AlignerContext::beginIteration(): unknown file type %d for '%s'\n", options->outputFile.fileType, options->outputFile.fileName);
            soft_exit(1);
        }
        format->setupReaderContext(options, &readerContext);

        writerSupplier = format->getWriterSupplier(options, readerContext.genome);
        ReadWriter* headerWriter = writerSupplier->getWriter();
        headerWriter->writeHeader(readerContext, options->sortOutput, argc, argv, version, options->rgLineContents, options->outputFile.omitSQLines);
        headerWriter->close();
        delete headerWriter;
    }
}

    void
AlignerContext::finishIteration()
{
    extension->finishIteration();

    if (NULL != writerSupplier) {
        if (options->dropIndexBeforeSort) {
            g_index->dropIndex();

        }

        writerSupplier->close();  // This is where the sort happens
        delete writerSupplier;
        writerSupplier = NULL;

        if (options->dropIndexBeforeSort) {
            //
            // Since we dropped the index part of the index, now we need to delete it completely.  We can't do it earlier because
            // sort needs the contigs from the genome.
            //
            delete g_index;
            g_index = NULL;
            index = NULL;
            delete g_indexDirectory;
            g_indexDirectory = NULL;
        }
    }

    alignTime = /*timeInMillis() - alignStart -- use the time from ParallelTask.h, that may exclude memory allocation time*/ time;
}

    bool
AlignerContext::nextIteration()
{
    //
    // This thing is a vestage of when we used to allow parameter ranges.
    //
    typeSpecificNextIteration();
    return false;
}

//
// Take an integer and a percentage, and turn it into a string of the form "number (percentage%)<padding>" where
// number has commas and the whole thing is padded out with spaces to a specific length.
//
char *numPctAndPad(char *buffer, _uint64 num, double pct, size_t desiredWidth, size_t bufferLen)
{
    _ASSERT(desiredWidth + 1 < bufferLen); // < to leave room for trailing null.

    FormatUIntWithCommas(num, buffer, bufferLen);
    const size_t percentageBufferSize = 100; // Plenty big enough for any value
    char percentageBuffer[percentageBufferSize];

    snprintf(percentageBuffer, percentageBufferSize, " (%.02f%%)", pct);
    if (strlen(percentageBuffer) + strlen(buffer) >= bufferLen || desiredWidth >= bufferLen) { // >= accounts for terminating null
        WriteErrorMessage("numPctAndPad: overflowed output buffer\n");
        buffer[0] = '\0';
        return buffer;
    }

    strcat(buffer, percentageBuffer);
    for (size_t x = strlen(buffer); x < desiredWidth; x++) {
        strcat(buffer, " ");
    }

    return buffer;
}

char *pctAndPad(char * buffer, double pct, size_t desiredWidth, size_t bufferLen, bool useDecimal, bool printPercentSign = true)
{
    _ASSERT(desiredWidth + 1 < bufferLen);

    const size_t percentageBufferSize = 100; // Plenty big enough for any value
    char percentageBuffer[percentageBufferSize];

    if (useDecimal) {
        snprintf(percentageBuffer, percentageBufferSize, "%.02f%s", pct, (printPercentSign ? "%" : ""));
    } else {
        snprintf(percentageBuffer, percentageBufferSize, "%d%s",  (unsigned)((100.0 * pct) + .5), (printPercentSign ? "%" : ""));
        snprintf(percentageBuffer, percentageBufferSize, "%d%s",  (unsigned)((100.0 * pct) + .5), (printPercentSign ? "%" : ""));
    }

    if (strlen(percentageBuffer) + 1 > bufferLen) {
        WriteErrorMessage("pctAndPad: buffer too small\n");
        buffer[0] = '\0';
        return buffer;
    }

    strcpy(buffer, percentageBuffer);
    for (size_t x = strlen(buffer); x < desiredWidth; x++) {
        strcat(buffer, " ");
    }

    return buffer;
}

    void
AlignerContext::printStats()
{
    if (NULL == statFile)
    {
        WriteStatusMessage("Total Reads    rRNA Reads           Aligned, MAPQ >= %2d    Aligned, MAPQ < %2d     Unaligned              Too Short/Too Many Ns  %s%s%sReads/s   Time in Aligner (s)%s%s\n", MAPQ_LIMIT_FOR_SINGLE_HIT, MAPQ_LIMIT_FOR_SINGLE_HIT,
            (stats->filtered > 0) ? "Filtered               " : "",
            (stats->extraAlignments) ? "Extra Alignments  " : "",
            isPaired() ? "%Pairs    " : "   ",
            options->profile ? (!options->sortOutput ? " Read Align Write(& compress)" : " Read Align Write") : "",
            (isPaired() && options->profileAffineGap) ? " %AgSingle %AgUsedSingle AG/Edit" : ""
            );
    }

    const size_t strBufLen = 50; // Way more than enough for 64 bit numbers with commas
    char tooShort[strBufLen];
    char rrna[strBufLen];
    char single[strBufLen];
    char multi[strBufLen];
    char unaligned[strBufLen];
    char numReads[strBufLen];
    char readsPerSecond[strBufLen];
    char alignTimeString[strBufLen];

    char filtered[strBufLen];
    char extraAlignments[strBufLen];
    char pctPairs[strBufLen];
    char pctRead[strBufLen];
    char pctAlign[strBufLen];
    char pctWrite[strBufLen];
    char pctAg[strBufLen];    
    char pctAg2[strBufLen];    
    char agRatio[strBufLen];
    _int64 totalTime = stats->millisReading + stats->millisAligning + stats->millisWriting;

    /*
                        total
                        |  rrna
                        |  |  single
                        |  |  |  multi
                        |  |  |  |  unaligned        reads/s
                        |  |  |  |  |  too short     |  time
                        |  |  |  |  |  |  filtered   |  | %Read
                        |  |  |  |  |  |  | extra    |  | | %Align
                        |  |  |  |  |  |  | | pairs  |  | | | %Write
                        |  |  |  |  |  |  | | |      |  | | | | Ag
                        |  |  |  |  |  |  | | |      |  | | | | | AgUsed
                        v  v  v  v  v  v  v v v      v  v v v v v v v AG/Edit
    */
    if (NULL == statFile)
    {
        WriteStatusMessage("%s %s %s %s %s %s %s%s%s   %-9s %s%s%s%s%s%s%s\n",
            FormatUIntWithCommas(stats->totalReads, numReads, strBufLen, 14),
            numPctAndPad(rrna, stats->rrnaReads, 100.0 * stats->rrnaReads / stats->totalReads, 20, strBufLen),
            numPctAndPad(single, stats->singleHits, 100.0 * stats->singleHits / stats->totalReads, 22, strBufLen),
            numPctAndPad(multi, stats->multiHits, 100.0 * stats->multiHits / stats->totalReads, 22, strBufLen),
            numPctAndPad(unaligned, stats->notFound, 100.0 * stats->notFound / stats->totalReads, 22, strBufLen),
            numPctAndPad(tooShort, stats->uselessReads , 100.0 * stats->uselessReads / max(stats->totalReads, (_int64)1), 22, strBufLen),
            (stats->filtered > 0) ? numPctAndPad(filtered, stats->filtered, 100.0 * stats->filtered / stats->totalReads, 23, strBufLen) : "",
            (stats->extraAlignments > 0) ? FormatUIntWithCommas(stats->extraAlignments, extraAlignments, strBufLen, 18) : "",
            isPaired() ? pctAndPad(pctPairs,  100.0 * stats->alignedAsPairs / stats->totalReads, 7, strBufLen, true) : "",
            FormatUIntWithCommas((_uint64)(1000 * stats->totalReads / max(alignTime, (_int64)1)), readsPerSecond, strBufLen), // Aligntime is in ms
            FormatUIntWithCommas((alignTime + 500) / 1000, alignTimeString, strBufLen, 20),
            options->profile ? pctAndPad(pctRead, (double)stats->millisReading / (double)totalTime, 5, strBufLen, false) : "",
            options->profile ? pctAndPad(pctAlign, (double)stats->millisAligning / (double)totalTime, 6, strBufLen, false) : "",
            options->profile ? pctAndPad(pctWrite, (double)stats->millisWriting / (double)totalTime, 6, strBufLen, false) : "",
            (isPaired() && options->profileAffineGap) ? pctAndPad(pctAg, (double)stats->agForcedSingleEndAlignment / (double)stats->totalReads, 10, strBufLen, true) : "",
            (isPaired() && options->profileAffineGap) ? pctAndPad(pctAg2, (double)stats->agUsedSingleEndAlignment / (double)stats->totalReads, 14, strBufLen, true) : "",
            options->profileAffineGap ? pctAndPad(agRatio, (double)stats->affineGapCalls / (double)stats->lvCalls * 100, 8, strBufLen, true, true) : ""
        );
        // hsk rds, cov, dep
        WriteStatusMessage("HSK_READS_COV_AVGDEP\t%" PRId64 "\t%.2f\t%.2f\n",
                              stats->hskReads,
                              100.0 * stats->hskCov.size() / T2T_HSK_SIZE,
                              1.0 * stats->hskBases / T2T_HSK_SIZE);
    }
    // sort IC by reads number and index
    std::vector<_int64> ic;
    std::unordered_map<_int8, _int64>::iterator it;
    for (it = stats->icReads.begin(); it != stats->icReads.end(); ++it)
        if (it->second)
            ic.insert(ic.end(), (it->second << 32) | (UINT32_MAX - it->first));
    sort(ic.begin(), ic.end(), std::greater<_int64>());
    std::vector<_int64>::iterator x;
    if (NULL == statFile)
    {
        for (x = ic.begin(); x != ic.end(); ++x)
            WriteStatusMessage("IC%d:%" PRId64 ";", UINT32_MAX - (_int32)*x, *x>>32);
        if (ic.size())
            WriteStatusMessage("\n");
    }
    if (NULL != perfFile) {
        fprintf(perfFile, "maxHits\tmaxDist\t%% reads not useless\t%% reads single hit\t%% reads multi hit\t%% reads not found\tLV calls\taffine gap calls\t%% aligned as pairs\ttotal reads\treads/s\n");

        char lvBuf[strBufLen];
        char agBuf[strBufLen];
        char totalReadsBuf[strBufLen];
        char timePerReadBuf[strBufLen];
        fprintf(perfFile, "%d\t%d\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%s\t%s\t%0.2f%%\t%s\t%s\n",
                maxHits_, maxDist_, 
                100.0 * (stats->totalReads - stats->uselessReads) / max(stats->totalReads, (_int64) 1),
                100.0 * stats->singleHits / stats->totalReads,
                100.0 * stats->multiHits / stats->totalReads,
                100.0 * stats->notFound / stats->totalReads,
                FormatUIntWithCommas(stats->lvCalls, lvBuf, strBufLen),
                FormatUIntWithCommas(stats->affineGapCalls, agBuf, strBufLen),
                100.0 * stats->alignedAsPairs / stats->totalReads,
                FormatUIntWithCommas(stats->totalReads, totalReadsBuf, strBufLen),
                FormatUIntWithCommas((1000 * (stats->totalReads - stats->uselessReads)) / max(alignTime, (_int64)1), timePerReadBuf, strBufLen));

        fprintf(perfFile,"\n");
    }
    // output stats to json via -ss
    if (NULL != statFile) {
        fputs("{\n", statFile);
        fputs("\t\"alignment\": {\n", statFile);
        fprintf(statFile, "\t\t\"total_reads\": %" PRId64 ",\n", stats->totalReads);
        fprintf(statFile, "\t\t\"rrna_reads\": %" PRId64 ",\n", stats->rrnaReads);
        fprintf(statFile, "\t\t\"rrna_pct\": %.3f,\n", stats->totalReads ? 100.0 * stats->rrnaReads / stats->totalReads : 0);
        fprintf(statFile, "\t\t\"aligned_mq10+\": %" PRId64 ",\n", stats->singleHits);
        fprintf(statFile, "\t\t\"aligned_mq10-\": %" PRId64 ",\n", stats->multiHits);
        fprintf(statFile, "\t\t\"unaligned\": %" PRId64 ",\n", stats->notFound);
        fprintf(statFile, "\t\t\"too_short_or_too_much_N\": %" PRId64 ",\n", stats->uselessReads);
        fprintf(statFile, "\t\t\"filtered\": %" PRId64 "\n", stats->filtered);
        fputs("\t},\n", statFile);
        fputs("\t\"hsk\": {\n", statFile);
        fprintf(statFile, "\t\t\"reads\": %" PRId64 ",\n", stats->hskReads);
        fprintf(statFile, "\t\t\"coverage\": %.3f,\n", 100.0 * stats->hskCov.size() / T2T_HSK_SIZE);
        fprintf(statFile, "\t\t\"depth\": %.3f\n", 1.0 * stats->hskBases / T2T_HSK_SIZE);
        if (!ic.size())
            fputs("\t}\n", statFile);
        else
        {
            fputs("\t},\n", statFile);
            fputs("\t\"ic\": {\n", statFile);
            if (ic.size() == 1)
            {
                x = ic.begin();
                fprintf(statFile, "\t\t\"%d\": %" PRId64 "\n", UINT32_MAX - (_int32)*x, *x>>32);
            }
            else
            {
                for (x = ic.begin(); x != ic.end() - 1; ++x)
                    fprintf(statFile, "\t\t\"%d\": %" PRId64 ",\n", UINT32_MAX - (_int32)*x, *x>>32);
                fprintf(statFile, "\t\t\"%d\": %" PRId64 "\n", UINT32_MAX - (_int32)*x, *x>>32);
            }
            fputs("\t}\n", statFile);
        }
        fputc('}', statFile);
    }


#if TIME_HISTOGRAM
    if (stats->backwardsTimeStamps != 0) {
        //
        // I'm pretty sure this is fixed, so don't print unless we actually see it it.
        //
        WriteStatusMessage("\n%lld time stamps were negative, for a total of %lld ns.  They are otherwise ignored.\n", stats->backwardsTimeStamps, stats->totalBackwardsTimeStamps);
    }

    WriteStatusMessage("\nPer-read alignment time histogram:\nlog2(ns)\tcount\ttotalTime(ns)\ttimePerRead\tReads\tTime\n");
    _int64 totalReads = 0;
    _int64 totalTimeX = 0;  // totalTime is already used
    for (int i = 0; i < 31; i++) {
        totalReads += stats->countByTimeBucket[i];
        totalTimeX += stats->nanosByTimeBucket[i];
    }

    _int64 readsSoFar = 0;
    _int64 timeSoFar = 0;
    for (int i = 0; i < 31; i++) {
        readsSoFar += stats->countByTimeBucket[i];
        timeSoFar += stats->nanosByTimeBucket[i];
        _int64 timePerRead;
        if (stats->countByTimeBucket[i] == 0) {
            timePerRead = 0;
        } else {
            timePerRead = stats->nanosByTimeBucket[i] / stats->countByTimeBucket[i];
        }

        WriteStatusMessage("%d\t%lld\t%lld\t%lld\t%f\t%f\n", i, stats->countByTimeBucket[i], stats->nanosByTimeBucket[i], timePerRead, (double)readsSoFar / (double)totalReads, (double)timeSoFar/(double)totalTimeX);
    }

    WriteStatusMessage("\nPer-read alignment count and time by MAPQ\nMAPQ\tcount\ttotalTime(ns)\ttimePerRead\tReads\tTime\tTime Per Read (right scale)\n");
    totalReads = stats->countOfUnaligned;
    totalTimeX = stats->timeOfUnaligned;
    for (int i = 0; i <= 70; i++) {
        totalReads += stats->countByMAPQ[i];
        totalTimeX += stats->timeByMAPQ[i];
    }
    readsSoFar = stats->countOfUnaligned;
    timeSoFar = stats->timeOfUnaligned;


    WriteStatusMessage("*\t%lld\t%lld\t%lld\t%f\t%f\t%lld\n", stats->countOfUnaligned, stats->timeOfUnaligned, stats->timeOfUnaligned / stats->countOfUnaligned, 
        (double)readsSoFar / (double)totalReads, (double)timeSoFar / (double)totalTimeX, stats->timeOfUnaligned / stats->countOfUnaligned / 1000000);
    for (int i = 0; i <= 70; i++) {
        readsSoFar += stats->countByMAPQ[i];
        timeSoFar += stats->timeByMAPQ[i];
        _int64 timePerRead;
        if (stats->countByMAPQ[i] == 0) {
            timePerRead = 0;
        } else {
            timePerRead = stats->timeByMAPQ[i] / stats->countByMAPQ[i];
        }
        WriteStatusMessage("%d\t%lld\t%lld\t%lld\t%f\t%f\t%lld\n", i, stats->countByMAPQ[i], stats->timeByMAPQ[i], timePerRead, (double)readsSoFar / (double)totalReads, (double)timeSoFar / (double)totalTimeX, timePerRead / 1000000);
    }

    WriteStatusMessage("\nPer-read alignment count and time by final edit distance\nEditDistance\tcount\ttotalTime(ns)\ttimePerRead\tReads\tTime\tTime per read (right scale)\n");
    totalReads = stats->countOfUnaligned;
    totalTimeX = stats->timeOfUnaligned;
    for (int i = 0; i < 31; i++) {
        totalReads += stats->countByNM[i];
        totalTimeX += stats->timeByNM[i];
    }
    readsSoFar = 0;
    timeSoFar = 0;
    for (int i = 0; i < 31; i++) {
        readsSoFar += stats->countByNM[i];
        timeSoFar += stats->timeByNM[i];
        _int64 timePerRead;
        if (stats->countByNM[i] == 0) {
            timePerRead = 0;
        } else {
            timePerRead = stats->timeByNM[i] / stats->countByNM[i];
        }
        WriteStatusMessage("%d\t%lld\t%lld\t%lld\t%f\t%f\t%lld\n", i, stats->countByNM[i], stats->timeByNM[i], timePerRead, (double)readsSoFar / (double)totalReads, (double)timeSoFar / (double)totalTimeX, timePerRead / 1000000);
    }
    readsSoFar += stats->countOfUnaligned;
    timeSoFar += stats->timeOfUnaligned;
    WriteStatusMessage("*\t%lld\t%lld\t%lld\t%f\t%f\t%lld\n", stats->countOfUnaligned, stats->timeOfUnaligned, stats->timeOfUnaligned / stats->countOfUnaligned, 
        (double)readsSoFar / (double)totalReads, (double)timeSoFar / (double)totalTimeX, stats->timeOfUnaligned / stats->countOfUnaligned / 1000000);

#endif // TIME_HISTOGRAM


    stats->printHistograms(stdout);

#ifdef  TIME_STRING_DISTANCE
    WriteStatusMessage("%llds, %lld calls in BSD noneClose, not -1\n",  stats->nanosTimeInBSD[0][1]/1000000000, stats->BSDCounts[0][1]);
    WriteStatusMessage("%llds, %lld calls in BSD noneClose, -1\n",      stats->nanosTimeInBSD[0][0]/1000000000, stats->BSDCounts[0][0]);
    WriteStatusMessage("%llds, %lld calls in BSD close, not -1\n",      stats->nanosTimeInBSD[1][1]/1000000000, stats->BSDCounts[1][1]);
    WriteStatusMessage("%llds, %lld calls in BSD close, -1\n",          stats->nanosTimeInBSD[1][0]/1000000000, stats->BSDCounts[1][0]);
    WriteStatusMessage("%llds, %lld calls in Hamming\n",                stats->hammingNanos/1000000000,         stats->hammingCount);
#endif  // TIME_STRING_DISTANCE

    extension->printStats();
}



        AlignerOptions*
AlignerContext::parseOptions(
    int i_argc,
    const char **i_argv,
    const char *i_version,
    unsigned *argsConsumed,
    bool      paired)
{
    argc = i_argc;
    argv = i_argv;
    version = i_version;

    g_suppressStatusMessages = false;   // This is a global, so it would carry over between runs (either in daemon mode or comma syntax).  Reset it here to get the expected behavior.
    g_suppressErrorMessages = false;    // ditto

    AlignerOptions *options;

    if (paired) {
        options = new PairedAlignerOptions("snap-aligner paired <index-dir> <inputFile(s)> [<options>] where <input file(s)> is a list of files to process.\n");
    } else {
        options = new AlignerOptions("snap-aligner single <index-dir> <inputFile(s)> [<options>] where <input file(s)> is a list of files to process.\n");
    }

    options->extra = extension->extraOptions();
    if (argc < 3) {
        WriteErrorMessage("Too few parameters\n");
        options->usage();
        delete options;
        return NULL;
    }

    options->indexDir = argv[1];
    struct InputList {
        SNAPFile    input;
        InputList*  next;
    } *inputList = NULL;

    //
    // Now build the input array and parse options.
    //

    bool inputFromStdio = false;

    int i;
    int nInputs = 0;
    for (i = 2; i < argc; i++) { // Starting at 2 skips single/paired and the index

        if (',' == argv[i][0]  && '\0' == argv[i][1]) {
            i++;    // Consume the comma
            break;
        }

        int argsConsumed;
        SNAPFile input;
        if (SNAPFile::generateFromCommandLine(argv+i, argc-i, &argsConsumed, &input, paired, true)) {
            if (input.isStdio) {
                if (CommandPipe != NULL) {
                    WriteErrorMessage("You may not use stdin/stdout in daemon mode\n");
                    delete options;
                    return NULL;
                }

                if (inputFromStdio) {
                    WriteErrorMessage("You specified stdin ('-') specified for more than one input, which isn't permitted.\n");
                    delete options;
                    return NULL;
                } else {
                    inputFromStdio = true;
                }
            }

            InputList *listEntry = new InputList;
            listEntry->input = input;
            listEntry->next = inputList;
            inputList = listEntry;      // Yes, this puts them in backwards.  a) We reverse them at the end and b) it doesn't matter anyway

            nInputs++;
            i += argsConsumed - 1;
            continue;
        }

        bool done;
        int oldI = i;
        if (!options->parse(argv, argc, i, &done)) {
            WriteErrorMessage("Didn't understand options starting at %s\n", argv[oldI]);
            options->usage();
            delete options;
            return NULL;
        }

        if (done) {
            i++;    // For the ',' arg
            break;
        }
    }

    if (0 == nInputs) {
        WriteErrorMessage("No input files specified.\n");
        delete options;
        return NULL;
    }

    if (options->maxDist + options->extraSearchDepth >= MAX_K) {
        WriteErrorMessage("You specified too large of a maximum edit distance combined with extra search depth.  The must add up to less than %d.\n", MAX_K);
        WriteErrorMessage("Either reduce their sum, or change MAX_K in LandauVishkin.h and recompile.\n");
        delete options;
        return NULL;
    }

    if (options->maxSecondaryAlignmentAdditionalEditDistance > (int)options->extraSearchDepth) {
        WriteErrorMessage("You can't have the max edit distance for secondary alignments (-om) be bigger than the max search depth (-D)\n");
        delete options;
        return NULL;
    }

    options->nInputs = nInputs;
    options->inputs = new SNAPFile[nInputs];
    for (int j = nInputs - 1; j >= 0; j --) {
        // The loop runs backwards so that we reverse the reversing that we did when we built it.  Not that it matters anyway.
        _ASSERT(NULL != inputList);
        options->inputs[j] = inputList->input;
        InputList *dying = inputList;
        inputList = inputList->next;
        delete dying;
    }
    _ASSERT(NULL == inputList);

    *argsConsumed = i;
    return options;
}
