/*++

Module Name:

    BaseAligner.h

Abstract:

    Header for SNAP genome aligner

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once

#include "AlignmentResult.h"
#include "LandauVishkin.h"
#include "AffineGap.h"
#include "AffineGapVectorized.h"
#include "BigAlloc.h"
#include "ProbabilityDistance.h"
#include "AlignerStats.h"
#include "directions.h"
#include "GenomeIndex.h"
#include "AlignmentAdjuster.h"
#include "AlignerOptions.h"

extern bool doAlignerPrefetch;

class BaseAligner {
public:

    BaseAligner(
        GenomeIndex             *i_genomeIndex, 
        unsigned                 i_maxHitsToConsider, 
        unsigned                 i_maxK,
        unsigned                 i_maxReadSize,
        unsigned                 i_maxSeedsToUse,
        double                   i_maxSeedCoverage,
        unsigned                 i_minWeightToCheck,
        unsigned                 i_extraSearchDepth,
        DisabledOptimizations    i_disabledOptimizations,
        bool                     i_useAffineGap,           
        bool                     i_ignoreAlignmentAdjustmentsForOm,
        bool                     i_altAwareness,
        bool                     i_emitALTAlignments,
        int                      i_maxScoreGapToPreferNonAltAlignment,
        int                      i_maxSecondaryAlignmentsPerContig,
        LandauVishkin<1>        *i_landauVishkin = NULL,
        LandauVishkin<-1>       *i_reverseLandauVishkin = NULL,
        unsigned                 i_matchReward = 1,
        unsigned                 i_subPenalty = 4,
        unsigned                 i_gapOpenPenalty = 6,
        unsigned                 i_gapExtendPenalty = 1,
        unsigned                 i_fivePrimeEndBonus = 10,
        unsigned                 i_threePrimeEndBonus = 5,
        AlignerStats            *i_stats = NULL,
        BigAllocator            *allocator = NULL);

    virtual ~BaseAligner();

        bool
    AlignRead(
        Read                    *read,
        SingleAlignmentResult   *primaryResult,
        SingleAlignmentResult   *firstALTResult,                        // This only gets filled in if there's a good ALT result that's not the primary result.
        int                      maxEditDistanceForSecondaryResults,
        _int64                   secondaryResultBufferSize,
        _int64                  *nSecondaryResults,
        _int64                   maxSecondaryResults,         // The most secondary results to return; always return the best ones
        SingleAlignmentResult   *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by AlignRead()
        _int64                   maxCandidatesForAffineGapBufferSize,
        _int64                  *nCandidatesForAffineGap,
        SingleAlignmentResult   *candidatesForAffineGap, // Alignment candidates that need to be rescored using affine gap
        bool                     useHamming  = false           // run Hamming distance based scoring instead of Landau-Vishkin
    );      // Retun value is true if there was enough room in the secondary alignment buffer for everything that was found.


        bool
    alignAffineGap(
        Read* read,
        SingleAlignmentResult* result,
        SingleAlignmentResult* firstALTResult,
        _int64 nCandidatesForAffineGap,
        SingleAlignmentResult* candidatesForAffineGap // Alignment candidates that need to be rescored using affine gap
    );
    //
    // Statistics gathering.
    //

    _int64 getNHashTableLookups() const {return nHashTableLookups;}
    _int64 getLocationsScoredWithLandauVishkin() const { return nLocationsScoredWithLandauVishkin; }
    _int64 getLocationsScoredWithAffineGap() const { return nLocationsScoredWithAffineGap; }
    _int64 getNHitsIgnoredBecauseOfTooHighPopularity() const {return nHitsIgnoredBecauseOfTooHighPopularity;}
    _int64 getNReadsIgnoredBecauseOfTooManyNs() const {return nReadsIgnoredBecauseOfTooManyNs;}
    _int64 getNIndelsMerged() const {return nIndelsMerged;}
    void addIgnoredReads(_int64 newlyIgnoredReads) {nReadsIgnoredBecauseOfTooManyNs += newlyIgnoredReads;}

    const char *getRCTranslationTable() const {return rcTranslationTable;}

    inline int getMaxK() const {return maxK;}

    inline void setMaxK(int maxK_) {maxK = maxK_;}

    inline void setReadId(int readId_) {readId = readId_;}

    const char *getName() const {return "Base Aligner";}

    inline bool checkedAllSeeds() {return popularSeedsSkipped == 0;}

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

    void *operator new(size_t size, BigAllocator *allocator) {_ASSERT(size == sizeof(BaseAligner)); return allocator->allocate(size);}
    void operator delete(void *ptr, BigAllocator *allocator) {/* do nothing.  Memory gets cleaned up when the allocator is deleted.*/}
 
    inline bool getExplorePopularSeeds() {return explorePopularSeeds;}
    inline void setExplorePopularSeeds(bool newValue) {explorePopularSeeds = newValue;}

    inline bool getStopOnFirstHit() {return stopOnFirstHit;}
    inline void setStopOnFirstHit(bool newValue) {stopOnFirstHit = newValue;}

    static size_t getBigAllocatorReservation(GenomeIndex *index, bool ownLandauVishkin, unsigned maxHitsToConsider, unsigned maxReadSize, unsigned seedLen, 
        unsigned numSeedsFromCommandLine, double seedCoverage, int maxSecondaryAlignmentsPerContig, unsigned extraSearchDepth);

    static const unsigned UnusedScoreValue = 0xffff;


private:

    inline int scoreLimit(bool forALT);

    bool hadBigAllocator;

    LandauVishkin<> *landauVishkin;
    LandauVishkin<-1> *reverseLandauVishkin;
    bool ownLandauVishkin;
	bool altAwareness;
    bool emitALTAlignments;
    int maxScoreGapToPreferNonAltAlignment;

    // AffineGap<> *affineGap;
    // AffineGap<-1> *reverseAffineGap;

    AffineGapVectorized<> *affineGap;
    AffineGapVectorized<-1> *reverseAffineGap;

    // Affine gap scoring parameters
    int matchReward;
    int subPenalty;
    int gapOpenPenalty;
    int gapExtendPenalty;

    ProbabilityDistance *probDistance;

    AlignmentAdjuster alignmentAdjuster;

    // Maximum distance to merge candidates that differ in indels over.
#ifdef LONG_READS
    static const unsigned maxMergeDist = 64; // Must be even and <= 64
#else
    static const unsigned maxMergeDist = 48; // Must be even and <= 64
#endif
    char rcTranslationTable[256];

    _int64 nHashTableLookups;
    _int64 nLocationsScoredWithLandauVishkin;
    _int64 nLocationsScoredWithAffineGap;
    _int64 nHitsIgnoredBecauseOfTooHighPopularity;
    _int64 nReadsIgnoredBecauseOfTooManyNs;
    _int64 nIndelsMerged;

    //
    // A bitvector indexed by offset in the read indicating whether this seed is used.
    // This is here to avoid doing a memory allocation in the aligner.
    //
    BYTE *seedUsed;
    BYTE *seedUsedAsAllocated;  // Use this for deleting seedUsed.

    inline bool IsSeedUsed(unsigned indexInRead) const {
        return (seedUsed[indexInRead / 8] & (1 << (indexInRead % 8))) != 0;
    }

    inline void SetSeedUsed(unsigned indexInRead) {
        seedUsed[indexInRead / 8] |= (1 << (indexInRead % 8));
    }

    struct Candidate {
        Candidate() {init();}
        void init();

        unsigned        score;
        int             seedOffset;
        double          matchProbability;
        GenomeLocation  origGenomeLocation;
    };

    static const unsigned hashTableElementSize = maxMergeDist;   // The code depends on this, don't change it

    void decomposeGenomeLocation(GenomeLocation genomeLocation, _uint64 *highOrder, _uint64 *lowOrder)
    {
        *lowOrder = (_uint64)GenomeLocationAsInt64(genomeLocation) % hashTableElementSize;
        if (NULL != highOrder) {
            *highOrder = (_uint64)GenomeLocationAsInt64(genomeLocation) - *lowOrder;
        }
    }

    struct HashTableElement {
        HashTableElement();
        void init();

        //
        // Doubly linked list for the weight buckets.
        //
        HashTableElement    *weightNext;
        HashTableElement    *weightPrev;

        //
        // Singly linked list for the hash table buckets.
        //
        HashTableElement    *next;

        _uint64              candidatesUsed;    // Really candidates we still need to score
        _uint64              candidatesScored;

        GenomeLocation       baseGenomeLocation;
        unsigned             weight;
        unsigned             lowestPossibleScore;
        unsigned             bestScore;
        int                  bestAGScore;
        GenomeLocation       bestScoreGenomeLocation; // adjusted location after scoring
        GenomeLocation       bestScoreOrigGenomeLocation; // location before scoring
        Direction            direction;
        bool                 allExtantCandidatesScored;
        double               matchProbabilityForBestScore;
        bool                 usedAffineGapScoring;
        int                  basesClippedBefore;
        int                  basesClippedAfter;
        int                  agScore;
        int                  seedOffset;

        Candidate            candidates[hashTableElementSize];
    };

    struct ScoreSet {
        ScoreSet();
        void init();
        void init(SingleAlignmentResult* result);

        int bestScore;
        GenomeLocation bestScoreGenomeLocation;
        GenomeLocation bestScoreOrigGenomeLocation; // location before scoring
        Direction bestScoreDirection;
        bool bestScoreUsedAffineGapScoring;
        int bestScoreBasesClippedBefore;
        int bestScoreBasesClippedAfter;
        int bestScoreAGScore;
        int bestScoreSeedOffset;
        double bestScoreMatchProbability;


        double probabilityOfAllCandidates;
        double probabilityOfBestCandidate;

        void updateProbabilitiesForNearbyMatch(double probabilityOfMatchBeingReplaced);   // For the "nearby match" code
        void updateProbabilitiesForNewMatch(double newProbability, double matchProbabilityOfNearbyMatch);
        inline void updateProbabilityOfAllMatches(double oldProbability) {
            probabilityOfAllCandidates = __max(0, probabilityOfAllCandidates - oldProbability);
        }
        inline void updateProbabilityOfBestMatch(double newProbability) {
            probabilityOfBestCandidate = newProbability;
            probabilityOfAllCandidates += newProbability;
        }
        void updateBestScore(
            GenomeLocation genomeLocation,
            GenomeLocation origGenomeLocation,
            unsigned score, 
            bool useAffineGap, 
            int agScore, 
            double matchProbability, 
            unsigned int& lvScoresAfterBestFound, 
            BaseAligner::HashTableElement *elementToScore, 
            SingleAlignmentResult *secondaryResults, 
            _int64* nSecondaryResults, 
            _int64 secondaryResultBufferSize,
            bool anyNearbyCandidatesAlreadyScored,
            int maxEditDistanceForSecondaryResults, 
            bool *overflowedSecondaryBuffer,
            _int64 maxCandidatesForAffineGapBufferSize,
            _int64* nCandidatesForAffineGap,
            SingleAlignmentResult* candidatesForAffineGap,
            unsigned extraSearchDepth);

        bool updateBestScore(SingleAlignmentResult* result) {
            probabilityOfAllCandidates += result->matchProbability;
            if (result->agScore > bestScoreAGScore || (result->agScore == bestScoreAGScore && result->matchProbability > bestScoreMatchProbability)) {
                bestScore = result->score;
                bestScoreAGScore = result->agScore;
                bestScoreMatchProbability = result->matchProbability;
                bestScoreGenomeLocation = result->location;
                bestScoreOrigGenomeLocation = result->origLocation;
                bestScoreDirection = result->direction;
                bestScoreUsedAffineGapScoring = result->usedAffineGapScoring;
                bestScoreBasesClippedBefore = result->basesClippedBefore;
                bestScoreBasesClippedAfter = result->basesClippedAfter;
                bestScoreSeedOffset = result->seedOffset;
                return true;
            } else {
                return false;
            }
        }

        void fillInSingleAlignmentLocation(SingleAlignmentResult* result, int popularSeedsSkipped);
        void fillInSingleAlignmentResult(SingleAlignmentResult* result, int popularSeedsSkipped);
    };

    //
    // Clearing out all of the pointers in the hash tables is expensive relative to running
    // an alignment, because usually the table is much bigger than the number of entries in it.
    // So, we avoid that expense by simply not clearing out the table at all.  Instead, along with
    // the pointers we keep an epoch number.  There's a corresponding epoch number in the
    // BaseAligner object, and if the two differ then the hash table bucket is empty.  We increment
    // the epoch number in the BaseAligner at the beginning of each alignment, thus effectively
    // clearing the hash table from the last run.
    //
    struct HashTableAnchor {
        HashTableElement *element;
        _int64            epoch;
    };

    _int64 hashTableEpoch;

    unsigned nUsedHashTableElements;
    unsigned hashTableElementPoolSize;
    HashTableElement *hashTableElementPool;

    const HashTableElement emptyHashTableElement;

    unsigned candidateHashTablesSize;
    HashTableAnchor *candidateHashTable[NUM_DIRECTIONS];
    
    HashTableElement *weightLists;
    unsigned highestUsedWeightList;
    unsigned wrapCount;
    unsigned nAddedToHashTable;

    static inline _uint64 hash(_uint64 key) {
        key = key * 131;    // Believe it or not, we spend a long time computing the hash, so we're better off with more table entries and a dopey function.
        return key;
    }


    // MAPQ parameters, currently not set to match Mason.  Using #define because VC won't allow "static const double".
#define SNP_PROB  0.001
#define GAP_OPEN_PROB  0.001
#define GAP_EXTEND_PROB  0.5

    //
    // Storage that's used during a call to AlignRead, but that's also needed by the
    // score function.  Since BaseAligner is single threaded, it's easier just to make
    // them member variables than to pass them around.
    //
    unsigned lowestPossibleScoreOfAnyUnseenLocation[NUM_DIRECTIONS];
    unsigned currRoundLowestPossibleScoreOfAnyUnseenLocation[NUM_DIRECTIONS];
    unsigned mostSeedsContainingAnyParticularBase[NUM_DIRECTIONS];
    unsigned nSeedsApplied[NUM_DIRECTIONS];
    ScoreSet scoresForAllAlignments;
    ScoreSet scoresForNonAltAlignments;
    //unsigned scoreLimit;
    unsigned minScoreThreshold; // used in affine gap to elide scoring of missed seed hits
    unsigned lvScoresAfterBestFound;
    unsigned affineGapScores;
    unsigned affineGapScoresAfterBestFound;

    int firstPassSeedsNotSkipped[NUM_DIRECTIONS];
    unsigned highestWeightListChecked;

    double totalProbabilityByDepth[AlignerStats::maxMaxHits];

        bool
    score(
        bool                     forceResult,
        Read                    *read[NUM_DIRECTIONS],
        SingleAlignmentResult   *primaryResult,
        SingleAlignmentResult   *firstALTResult,                        // This only gets filled in if there's a good ALT result that's not the primary result.
        int                      maxEditDistanceForSecondaryResults,
        _int64                   secondaryResultBufferSize,
        _int64                  *nSecondaryResults,
        SingleAlignmentResult   *secondaryResults,
        bool                    *overflowedSecondaryResultsBuffer,
        _int64                   maxCandidatesForAffineGapBufferSize,
        _int64                  *nCandidatesForAffineGap,
        SingleAlignmentResult   *candidatesForAffineGap,
        bool                     useHamming = false);

    void scoreLocationWithAffineGap(
        Read*                read[NUM_DIRECTIONS],
        Direction            direction,
        GenomeLocation       genomeLocation,
        unsigned             seedOffset,
        int                  scoreLimit,
        int                 *score,
        double              *matchProbability,
        int                 *genomeLocationOffset,
        int                 *basesClippedBefore,
        int                 *basesClippedAfter,
        int                 *agScore
    );

    void clearCandidates();

    bool findElement(GenomeLocation genomeLocation, Direction direction, HashTableElement **hashTableElement);
    void findCandidate(GenomeLocation genomeLocation, Direction direction, Candidate **candidate, HashTableElement **hashTableElement);
    void allocateNewCandidate(GenomeLocation genomeLoation, Direction direction, unsigned lowestPossibleScore, int seedOffset, Candidate **candidate, HashTableElement **hashTableElement);
    void incrementWeight(HashTableElement *element);
    void prefetchHashTableBucket(GenomeLocation genomeLocation, Direction direction);

    const Genome            *genome;
    GenomeIndex             *genomeIndex;
    unsigned                 seedLen;
    unsigned                 maxHitsToConsider;
    unsigned                 maxK;
    unsigned                 maxReadSize;
    unsigned                 maxSeedsToUseFromCommandLine; // Max number of seeds to look up in the hash table
    double                   maxSeedCoverage;  // Max seeds to used expressed as readSize/seedSize this is mutually exclusive with maxSeedsToUseFromCommandLine
    unsigned                 minWeightToCheck;
    unsigned                 extraSearchDepth;
    unsigned                 numWeightLists;
    DisabledOptimizations    disabledOptimizations;
    bool                     useAffineGap;
    bool                     ignoreAlignmentAdjustmentsForOm;
    bool                     doesGenomeIndexHave64BitLocations;
    int                      maxSecondaryAlignmentsPerContig;

    struct HitsPerContigCounts {
        _int64  epoch;          // Used hashTableEpoch, for the same reason
        int     hits;
    };

    HitsPerContigCounts *hitsPerContigCounts;   // How many alignments are we reporting for each contig.  Used to implement -mpc, otheriwse unallocated.

    char *rcReadData;
    char *rcReadQuality;
    char *reversedRead[NUM_DIRECTIONS];

    unsigned nTable[256];

    int readId;
    
    // How many overly popular (> maxHits) seeds we skipped this run
    unsigned popularSeedsSkipped;

    bool explorePopularSeeds; // Whether we should explore the first maxHits hits even for overly
                              // popular seeds (useful for filtering reads that come from a database
                              // with many very similar sequences).

    bool stopOnFirstHit;      // Whether to stop the first time a location matches with less than
                              // maxK edit distance (useful when using SNAP for filtering only).

    AlignerStats *stats;

    unsigned *hitCountByExtraSearchDepth;   // How many hits at each depth bigger than the current best edit distance.
                                            // So if the current best hit has edit distance 2, then hitCountByExtraSearchDepth[0] would
                                            // be the count of hits at edit distance 2, while hitCountByExtraSearchDepth[2] would be the count
                                            // of hits at edit distance 4.

    void finalizeSecondaryResults(
        Read                    *read,
        SingleAlignmentResult   *primaryResult,
        _int64                  *nSecondaryResults,                     // in/out
        SingleAlignmentResult   *secondaryResults,
        _int64                   maxSecondaryResults,
        int                      maxEditDistanceForSecondaryResults,
        int                      bestScore);
};
