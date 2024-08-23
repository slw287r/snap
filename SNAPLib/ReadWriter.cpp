/*++

Module Name:

    ReadWriter.cpp

Abstract:

    General file writer.

Environment:

    User mode service.

    Not thread safe.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "Read.h"
#include "SAM.h"
#include "Tables.h"
#include "RangeSplitter.h"
#include "ParallelTask.h"
#include "Util.h"
#include "ReadSupplierQueue.h"
#include "FileFormat.h"
#include "exit.h"
#include "Error.h"
#include "Genome.h"

class SimpleReadWriter : public ReadWriter
{
public:
    SimpleReadWriter(const FileFormat* i_format, DataWriter* i_writer, const Genome* i_genome, bool i_killIfTooSlow, bool i_emitInternalScore, 
                     char *i_internalScoreTag, bool i_ignoreAlignmentAdjustmentsForOm, int i_matchReward, int i_subPenalty, int i_gapOpenPenalty, 
                     int i_gapExtendPenalty, bool i_attachAlignmentTime)
        : format(i_format), writer(i_writer), genome(i_genome), killIfTooSlow(i_killIfTooSlow), lastTooSlowCheck(0), 
          emitInternalScore(i_emitInternalScore), ignoreAlignmentAdjustmentsForOm(i_ignoreAlignmentAdjustmentsForOm),
          attachAlignmentTime(i_attachAlignmentTime)
    {
        if (emitInternalScore) {
            if (strlen(i_internalScoreTag) != 2) {
                WriteErrorMessage("SimpleReadWriter: bogus internal score tag\n");
                soft_exit(1);
            }
            strcpy(internalScoreTag, i_internalScoreTag);
        } else  {
            internalScoreTag[0] = '\0';
        }

        agc.init(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty);
    }

    virtual ~SimpleReadWriter()
    {
        delete writer;
    }

	virtual bool writeHeader(const ReaderContext& context, bool sorted, int argc, const char **argv, const char *version, const char *rgLine, bool omitSQLines);

    virtual bool writeReads(const ReaderContext& context, Read *read, SingleAlignmentResult *results, _int64 nResults, bool firstIsPrimary, bool useAffineGap = false);

    virtual bool writePairs(const ReaderContext& context, Read **reads /* array of size 2 */, PairedAlignmentResult *result, _int64 nResults, 
        SingleAlignmentResult **singleResults /* array of size 2*/, _int64 *nSingleResults /* array of size 2*/, bool firstIsPrimary, bool useAffineGap = false);

    virtual void close();

private:
    const FileFormat* format;
    DataWriter* writer;
    const Genome* genome;
    LandauVishkinWithCigar lvc;
    AffineGapVectorizedWithCigar agc;

    void checkIfTooSlow();
    bool killIfTooSlow;
    _int64 lastTooSlowCheck;
    _int64 writesSinceLastTooSlowCheck;

    bool emitInternalScore;
    char internalScoreTag[3];
    bool ignoreAlignmentAdjustmentsForOm;
    bool attachAlignmentTime;
};

    bool
SimpleReadWriter::writeHeader(
    const ReaderContext& context,
    bool sorted,
    int argc,
    const char **argv,
    const char *version,
    const char *rgLine,
	bool omitSQLines)
{
    char* buffer;
    size_t size;
    size_t used;

    char *localBuffer = NULL;

	writer->inHeader(true);
    if (! writer->getBuffer(&buffer, &size)) {
        return false;
    }

    char *writerBuffer = buffer;
    size_t writerBufferSize = size;

	while (!format->writeHeader(context, buffer, size, &used, sorted, argc, argv, version, rgLine, omitSQLines)) {
        delete[] localBuffer;
        size = 2 * size;
        localBuffer = new char[size];
        buffer = localBuffer;
    }

    if (NULL == localBuffer) {
        _ASSERT(writerBuffer == buffer);
        writer->advance((unsigned)used, 0);
        writer->nextBatch();
    } else {
        size_t bytesRemainingToWrite = used;
        size_t bytesWritten = 0;
        while (bytesRemainingToWrite > 0) {
            size_t bytesToWrite = __min(bytesRemainingToWrite, writerBufferSize);
            memcpy(writerBuffer, localBuffer + bytesWritten, bytesToWrite);
            writer->advance(bytesToWrite);
            writer->nextBatch();
            if (!writer->getBuffer(&writerBuffer, &writerBufferSize)) {
                return false;
            }
            bytesWritten += bytesToWrite;
            bytesRemainingToWrite -= bytesToWrite;
        }

        delete[] localBuffer;
    }

	writer->inHeader(false);
    return true;
}

    void
SimpleReadWriter::checkIfTooSlow()
{
    const _int64 tooSlowCheckPeriod = 5 * 60 * 1000;    // 5 min in ms
    const _int64 tooSlowCheckMinReadsPerCheckPeriod = 5 * 60 * 1000;    // One read/ms (or 1000 reads/s, but just on this thread).

    if (killIfTooSlow) {
        _int64 now = timeInMillis();
        if (lastTooSlowCheck + tooSlowCheckPeriod <= now) {
            if (lastTooSlowCheck != 0 && writesSinceLastTooSlowCheck < tooSlowCheckMinReadsPerCheckPeriod) {
                WriteErrorMessage("Only wrote %lld writes during a %lld minute check period; we're probably out of memory and are giving up because of -kts\n", writesSinceLastTooSlowCheck, tooSlowCheckPeriod / (60 * 1000));
                soft_exit(1);
            }

            lastTooSlowCheck = now;
            writesSinceLastTooSlowCheck = 0;
        }

        writesSinceLastTooSlowCheck++;
    } // if (killIfTooSlow)
}


    bool
SimpleReadWriter::writeReads(
    const ReaderContext& context,
    Read *read,
    SingleAlignmentResult *results,
    _int64 nResults,
    bool firstIsPrimary,
    bool useAffineGap
	)
{
    char* buffer;
    size_t size;
    size_t used;
    bool result = false;

    checkIfTooSlow();

    for (_int64 i = 0; i < nResults; i++) {
        if (results[i].status == NotFound) {
            results[i].location = InvalidGenomeLocation;
        }
    }

    //
    // We need to keep track of the offsets of all of the alignments in the output buffer so we can commit them.  However,
    // we want to avoid dynamic memory allocation as much as possible.  So, we have a static buffer on the stack that's big enough
    // for the great majority of cases, and then allocate dynamically if that's too small.  Makes for annoying, but efficient
    // code.
    //

    const int staticUsedBufferSize = 2000;
    size_t staticUsedBuffer[staticUsedBufferSize];

    GenomeLocation staticFinalLocationsBuffer[staticUsedBufferSize];

    size_t *usedBuffer;
    GenomeLocation *finalLocations;
    if (nResults <= staticUsedBufferSize) {
        usedBuffer = staticUsedBuffer;
        finalLocations = staticFinalLocationsBuffer;
    } else {
        usedBuffer = new size_t[nResults];
        finalLocations = new GenomeLocation[nResults];
    }

    for (int pass = 0; pass < 2; pass++) { // Make two passes, one with whatever buffer space is left and one with a clean buffer.
        bool blewBuffer = false;

        if (!writer->getBuffer(&buffer, &size)) {
            goto done;
        }

        used = 0;

        for (_int64 whichResult = 0; whichResult < nResults; whichResult++) {
            int addFrontClipping = 0;
            read->setAdditionalFrontClipping(results[whichResult].clippingForReadAdjustment);

            int cumulativeAddFrontClipping = 0, cumulativeAddBackClipping = 0;
            finalLocations[whichResult] = results[whichResult].location;

            unsigned nAdjustments = 0;
            // stopOnFirstHit (via -f) leaves score to the default value 0xffff and hence > 0
            if (useAffineGap && (results[whichResult].usedAffineGapScoring || results[whichResult].score > 0)) {
                while (!format->writeRead(context, &agc, buffer + used, size - used, &usedBuffer[whichResult], read->getIdLength(), read, results[whichResult].status,
                    results[whichResult].mapq, finalLocations[whichResult], results[whichResult].direction, (whichResult > 0) || !firstIsPrimary, results[whichResult].supplementary, &addFrontClipping,
                    results[whichResult].score, results[whichResult].scorePriorToClipping, emitInternalScore, internalScoreTag, attachAlignmentTime, results[whichResult].alignmentTimeInNanoseconds,
                    results[whichResult].basesClippedBefore,
                    results[whichResult].basesClippedAfter)) {

                    _ASSERT(0 == addFrontClipping || ignoreAlignmentAdjustmentsForOm); // Because of the alignment adjuster.

                    nAdjustments++;

                    if (0 == addFrontClipping) {
                        blewBuffer = true;
                        break;
                    }

                    // redo if read modified (e.g. to add soft clipping, or move alignment for a leading I.
                    const Genome::Contig *originalContig = results[whichResult].status == NotFound ? NULL
                        : genome->getContigAtLocation(results[whichResult].location);
                    const Genome::Contig *newContig = results[whichResult].status == NotFound ? NULL
                        : genome->getContigAtLocation(results[whichResult].location + addFrontClipping);
                    if (newContig == NULL || newContig != originalContig || finalLocations[whichResult] + addFrontClipping > newContig->beginningLocation + newContig->length - genome->getChromosomePadding() ||
                        nAdjustments > read->getDataLength()) {
                        //
                        // Altering this would push us over a contig boundary, or we're stuck in a loop.  Just give up on the read.
                        //
                        results[whichResult].status = NotFound;
                        results[whichResult].location = InvalidGenomeLocation;
                        results[whichResult].score = -1;
                        results[whichResult].direction = FORWARD;
                        finalLocations[whichResult] = InvalidGenomeLocation;
                    } else {
                        if (addFrontClipping < 0) { // Insertion (soft-clip)
                            cumulativeAddFrontClipping += addFrontClipping;
                            if (results[whichResult].direction == FORWARD) {
                                read->setAdditionalFrontClipping(-cumulativeAddFrontClipping);
                            } else {
                                read->setAdditionalBackClipping(-cumulativeAddFrontClipping);
                            }
                        } else { // Deletion
                            finalLocations[whichResult] = results[whichResult].location + addFrontClipping;
                        }
                    }
                } // while formatting doesn't work
            } else {
                while (!format->writeRead(context, &lvc, buffer + used, size - used, &usedBuffer[whichResult], read->getIdLength(), read, results[whichResult].status,
                    results[whichResult].mapq, finalLocations[whichResult], results[whichResult].direction, (whichResult > 0) || !firstIsPrimary, results[whichResult].supplementary, &addFrontClipping,
                    results[whichResult].scorePriorToClipping, emitInternalScore, internalScoreTag, attachAlignmentTime, results[whichResult].alignmentTimeInNanoseconds)) {

                    _ASSERT(0 == addFrontClipping || ignoreAlignmentAdjustmentsForOm); // Because of the alignment adjuster.

                    nAdjustments++;

                    if (0 == addFrontClipping) {
                        blewBuffer = true;
                        break;
                    }

                    // redo if read modified (e.g. to add soft clipping, or move alignment for a leading I.
                    const Genome::Contig *originalContig = results[whichResult].status == NotFound ? NULL
                        : genome->getContigAtLocation(results[whichResult].location);
                    const Genome::Contig *newContig = results[whichResult].status == NotFound ? NULL
                        : genome->getContigAtLocation(results[whichResult].location + addFrontClipping);
                    if (newContig == NULL || newContig != originalContig || finalLocations[whichResult] + addFrontClipping > originalContig->beginningLocation + originalContig->length - genome->getChromosomePadding() ||
                        nAdjustments > read->getDataLength()) {
                        //
                        // Altering this would push us over a contig boundary, or we're stuck in a loop.  Just give up on the read.
                        //
                        results[whichResult].status = NotFound;
                        results[whichResult].location = InvalidGenomeLocation;
                        results[whichResult].score = -1;
                        results[whichResult].direction = FORWARD;
                        finalLocations[whichResult] = InvalidGenomeLocation;
                    } else {
                        if (addFrontClipping > 0) {
                            cumulativeAddFrontClipping += addFrontClipping;
                            read->setAdditionalFrontClipping(cumulativeAddFrontClipping);
                        }
                        finalLocations[whichResult] += addFrontClipping;
                    }
                } // while formatting doesn't work
            }
            if (blewBuffer) {
                break;
            }

            used += usedBuffer[whichResult];
            _ASSERT(used <= size);

            if (used > 0xffffffff) {
                 WriteErrorMessage("SimpleReadWriter:writeReads: used too big\n");
                 soft_exit(1);
            }
        } // for each result.

        if (!blewBuffer) {
            //
            // Everything worked OK.
            //
            for (int whichResult = 0; whichResult < nResults; whichResult++) {
                writer->advance((unsigned)usedBuffer[whichResult], finalLocations[whichResult]);
            }
            result = true;
            goto done;
        }

        if (pass == 1) {
            WriteErrorMessage("Failed to write into fresh buffer; trying providing the -wbs switch with a larger value\n");
            soft_exit(1);
        }

        if (!writer->nextBatch()) {
            goto done;
        }
    } // for each pass (i.e., not empty, empty buffer)
    
done:
    if (usedBuffer != staticUsedBuffer) {
        delete[] usedBuffer;
        usedBuffer = NULL;

        delete[] finalLocations;
        finalLocations = NULL;
    }

    read->setAdditionalFrontClipping(0);

    return result;
}

    bool
SimpleReadWriter::writePairs(
    const ReaderContext& context, 
    Read **reads /* array of size NUM_READS_PER_PAIR */, 
    PairedAlignmentResult *result, 
    _int64 nResults,
    SingleAlignmentResult **singleResults /* array of size NUM_READS_PER_PAIR*/, 
    _int64 *nSingleResults /* array of size NUM_READS_PER_PAIR*/, 
    bool firstIsPrimary,
    bool useAffineGap)
{
    bool retVal = false;
    //
    // We need to write all alignments for the pair into the same buffer, so that a write from
    // some other thread doesn't separate them.  We make two passes, trying to write into the 
    // existing buffer, and then into a clean one.  If that doesn't work, abort the alignment
    // run and ask for a bigger write buffer.
    //
    const int staticUsedBufferSize = 2000;
    size_t staticUsedBuffer[NUM_READS_PER_PAIR][staticUsedBufferSize];
    GenomeLocation staticLocationBuffer[NUM_READS_PER_PAIR][staticUsedBufferSize];

    _int64 NullNSingleResults[NUM_READS_PER_PAIR] = { 0,0 };
    if (NULL == nSingleResults) {
        nSingleResults = NullNSingleResults;
    }

    checkIfTooSlow();

    GenomeLocation *finalLocations[NUM_READS_PER_PAIR];
    size_t *usedBuffer[NUM_READS_PER_PAIR];
    if (nResults + nSingleResults[0] <= staticUsedBufferSize && nResults + nSingleResults[1] <= staticUsedBufferSize) {
        usedBuffer[0] = staticUsedBuffer[0];
        usedBuffer[1] = staticUsedBuffer[1];
        finalLocations[0] = staticLocationBuffer[0];
        finalLocations[1] = staticLocationBuffer[1];
    } else {
        usedBuffer[0] = new size_t[nResults * NUM_READS_PER_PAIR + nSingleResults[0] + nSingleResults[1]];
        usedBuffer[1] = usedBuffer[0] + nResults + nSingleResults[0];
        finalLocations[0] = new GenomeLocation[nResults * NUM_READS_PER_PAIR + nSingleResults[0] + nSingleResults[1]];
        finalLocations[1] = finalLocations[0] + nResults + nSingleResults[0];
    }


    //
    // For paired reads, we need to have the same QNAME for both of them, and it needs to be unique among all other
    // reads in the dataset.  For now, all we do is see if the read names end in /1 and /2, and if so truncate them.
    //
    size_t idLengths[NUM_READS_PER_PAIR];
    idLengths[0] = reads[0]->getIdLength();
    idLengths[1] = reads[1]->getIdLength();
    if (idLengths[0] == idLengths[1] && idLengths[0] > 2 && reads[0]->getId()[idLengths[0]-2] == '/' && reads[1]->getId()[idLengths[0]-2] == '/') {
        char lastChar0, lastChar1;
        lastChar0 = reads[0]->getId()[idLengths[0] - 1];
        lastChar1 = reads[1]->getId()[idLengths[1] - 1];
        if ((lastChar0 == '1' || lastChar0 == '2') && (lastChar1 == '1' || lastChar1 == '2') && 
            lastChar0 != lastChar1) {
                idLengths[0] -= 2;
                idLengths[1] -= 2;
        }
    }

    for (int pass = 0; pass < 2; pass++) {

        char* buffer;
        size_t size;
        size_t used = 0;

        bool fitInBuffer = true;

        if (!writer->getBuffer(&buffer, &size)) {
            goto done;
        }

        //
        // Write all of the pair alignments into the buffer.
        //
        for (int whichAlignmentPair = 0; whichAlignmentPair < nResults; whichAlignmentPair++) {
            reads[0]->setAdditionalFrontClipping(result[whichAlignmentPair].clippingForReadAdjustment[0]);
            reads[1]->setAdditionalFrontClipping(result[whichAlignmentPair].clippingForReadAdjustment[1]);

            GenomeLocation locations[2];
            locations[0] = result[whichAlignmentPair].status[0] != NotFound ? result[whichAlignmentPair].location[0] : InvalidGenomeLocation;
            locations[1] = result[whichAlignmentPair].status[1] != NotFound ? result[whichAlignmentPair].location[1] : InvalidGenomeLocation;

            bool secondReadLocationChanged;
            int cumulativePositiveAddFrontClipping[NUM_READS_PER_PAIR] = { 0, 0 };
            bool writeOrderChanged;

            do {
                size_t tentativeUsed = 0;
                secondReadLocationChanged = false;
                writeOrderChanged = false;


                int writeOrder[2];  // The order in which we write the reads, which is just numerical by genome location.  SO writeOrder[0] gets written first, and writeOrder[1] second.

                if (locations[0] <= locations[1]) {
                    writeOrder[0] = 0;
                    writeOrder[1] = 1;
                } else {
                    writeOrder[0] = 1;
                    writeOrder[1] = 0;
                }

                size_t spaceUsed[2] = { usedBuffer[0][whichAlignmentPair], usedBuffer[1][whichAlignmentPair] };
                bool outOfSpace = false;
                if (!format->writePairs(context, &lvc, &agc, useAffineGap, buffer + used + tentativeUsed, size - used - tentativeUsed, 
                    spaceUsed, idLengths, reads, locations, &result[whichAlignmentPair], whichAlignmentPair != 0 || !firstIsPrimary,
                    emitInternalScore, internalScoreTag, attachAlignmentTime, writeOrder, cumulativePositiveAddFrontClipping, 
                    &secondReadLocationChanged, &outOfSpace)) {
                    
                    if (outOfSpace) {
                        //
                        // We failed because we ran out of buffer.
                        //
                        goto blownBuffer;
                    }
                }

                // Check if the write order is correct, if not redo
                int newWriteOrder[2];
                if (locations[0] <= locations[1]) {
                    newWriteOrder[0] = 0;
                    newWriteOrder[1] = 1;
                } else {
                    newWriteOrder[0] = 1;
                    newWriteOrder[1] = 0;
                }

                if (writeOrder[0] != newWriteOrder[0] || writeOrder[1] != newWriteOrder[1]) {
                    writeOrder[0] = newWriteOrder[0];
                    writeOrder[1] = newWriteOrder[1];
                    writeOrderChanged = true;
                }

                usedBuffer[0][whichAlignmentPair] = spaceUsed[0];
                usedBuffer[1][whichAlignmentPair] = spaceUsed[1];

            } while (secondReadLocationChanged || writeOrderChanged);
            used += usedBuffer[0][whichAlignmentPair] + usedBuffer[1][whichAlignmentPair];

            //
            // Both reads are written into the buffer.  Save the final locations we used for when we commit.
            //
            for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
                finalLocations[whichRead][whichAlignmentPair] = locations[whichRead];
            }
        } // for each pair.

        //
        // Now write the single alignments.
        //
        for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            for (int whichAlignment = 0; whichAlignment < nSingleResults[whichRead]; whichAlignment++) {
                int addFrontClipping;
                reads[whichRead]->setAdditionalFrontClipping(singleResults[whichRead]->clippingForReadAdjustment);

                GenomeLocation location = singleResults[whichRead][whichAlignment].status != NotFound ? singleResults[whichRead][whichAlignment].location : InvalidGenomeLocation;
                int cumulativePositiveAddFrontClipping = 0;

                if (useAffineGap && (singleResults[whichRead][whichAlignment].usedAffineGapScoring || singleResults[whichRead][whichAlignment].score > 0)) {
                    while (!format->writeRead(context, &agc, buffer + used, size - used, &usedBuffer[whichRead][nResults + whichAlignment], reads[whichRead]->getIdLength(),
                        reads[whichRead], singleResults[whichRead][whichAlignment].status, singleResults[whichRead][whichAlignment].mapq, location, singleResults[whichRead][whichAlignment].direction,
                        true, singleResults[whichRead][whichAlignment].supplementary, &addFrontClipping, singleResults[whichRead][whichAlignment].score,
                        singleResults[whichRead][whichAlignment].scorePriorToClipping, emitInternalScore, internalScoreTag, attachAlignmentTime, singleResults[whichRead][whichAlignment].alignmentTimeInNanoseconds,
                        singleResults[whichRead][whichAlignment].basesClippedBefore, singleResults[whichRead][whichAlignment].basesClippedAfter)) {

                        if (0 == addFrontClipping) {
                            goto blownBuffer;
                        }

                        const Genome::Contig *originalContig = genome->getContigAtLocation(location);
                        const Genome::Contig *newContig = genome->getContigAtLocation(location + addFrontClipping);
                        if (newContig != originalContig || NULL == newContig || location + addFrontClipping > originalContig->beginningLocation + originalContig->length - genome->getChromosomePadding()) {
                            //
                            // Altering this would push us over a contig boundary.  Just give up on the read.
                            //
                            singleResults[whichRead][whichAlignment].status = NotFound;
                            location = InvalidGenomeLocation;
                            singleResults[whichRead][whichAlignment].score = -1;
                            singleResults[whichRead][whichAlignment].direction = FORWARD;
                        } else {
                            if (addFrontClipping < 0) { // Insertion (soft-clip)
                                cumulativePositiveAddFrontClipping += addFrontClipping;
                                if (singleResults[whichRead][whichAlignment].direction == FORWARD) {
                                    reads[whichRead]->setAdditionalFrontClipping(-cumulativePositiveAddFrontClipping);
                                }
                                else {
                                    reads[whichRead]->setAdditionalBackClipping(-cumulativePositiveAddFrontClipping);
                                }
                            }
                            else { // Deletion
                                location += addFrontClipping;
                            }
                        }
                    }
                } else {
                    while (!format->writeRead(context, &lvc, buffer + used, size - used, &usedBuffer[whichRead][nResults + whichAlignment], reads[whichRead]->getIdLength(),
                        reads[whichRead], singleResults[whichRead][whichAlignment].status, singleResults[whichRead][whichAlignment].mapq, location, singleResults[whichRead][whichAlignment].direction,
                        true, singleResults[whichRead][whichAlignment].supplementary, &addFrontClipping, singleResults[whichRead][whichAlignment].scorePriorToClipping, emitInternalScore, internalScoreTag,
                        attachAlignmentTime, singleResults[whichRead][whichAlignment].alignmentTimeInNanoseconds)) {

                        if (0 == addFrontClipping) {
                            goto blownBuffer;
                        }

                        const Genome::Contig *originalContig = genome->getContigAtLocation(location);
                        const Genome::Contig *newContig = genome->getContigAtLocation(location + addFrontClipping);
                        if (newContig != originalContig || NULL == newContig || location + addFrontClipping > originalContig->beginningLocation + originalContig->length - genome->getChromosomePadding()) {
                            //
                            // Altering this would push us over a contig boundary.  Just give up on the read.
                            //
                            singleResults[whichRead][whichAlignment].status = NotFound;
                            location = InvalidGenomeLocation;
                            singleResults[whichRead][whichAlignment].score = -1;
                            singleResults[whichRead][whichAlignment].direction = FORWARD;
                        }
                        else {
                            if (addFrontClipping > 0) {
                                cumulativePositiveAddFrontClipping += addFrontClipping;
                                reads[whichRead]->setAdditionalFrontClipping(cumulativePositiveAddFrontClipping);
                            }
                            location += addFrontClipping;
                        }
                    }
                }

                finalLocations[whichRead][nResults + whichAlignment] = location;
                used += usedBuffer[whichRead][nResults + whichAlignment];
            } // For each single alignment of a read
        } // For each read

        //
        // They all fit into the buffer.
        //

        //
        // Commit the updates for the pairs.
        //
        for (int whichReadPair = 0; whichReadPair < nResults; whichReadPair++) {
            for (int firstOrSecond = 0; firstOrSecond < NUM_READS_PER_PAIR; firstOrSecond++) {
                // adjust for write order
                int writeFirstOrSecond = (!!firstOrSecond) ^ (finalLocations[0][whichReadPair] > finalLocations[1][whichReadPair]); // goofy looking !! converts int to bool
                writer->advance((unsigned)usedBuffer[firstOrSecond][whichReadPair],
                    finalLocations[writeFirstOrSecond][whichReadPair] == InvalidGenomeLocation ? finalLocations[1 - writeFirstOrSecond][whichReadPair] : finalLocations[writeFirstOrSecond][whichReadPair]);
            }
        }

        //
        // Now commit the updates for the single reads.
        //
        for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            for (int whichAlignment = 0; whichAlignment < nSingleResults[whichRead]; whichAlignment++) {
                writer->advance((unsigned)usedBuffer[whichRead][nResults + whichAlignment], finalLocations[whichRead][nResults + whichAlignment]);
            }
        }

        retVal = true;
        break;

blownBuffer:
        if (pass > 0) {
            WriteErrorMessage("Unable to fit all alignments for one read pair into a single write buffer.  Increase the size of the write buffer with -wbs, or reduce the number of alignments with -om or -omax\n");
            WriteErrorMessage("Read id: '%.*s'\n", reads[0]->getIdLength(), reads[0]->getId());
            soft_exit(1);
        }

        if (!writer->nextBatch()) {
            goto done;
        }
            
    } // For each buffer full pass



done:
    if (usedBuffer[0] != staticUsedBuffer[0]) {
        delete[] usedBuffer[0];
        usedBuffer[0] = usedBuffer[1] = NULL;

        delete[] finalLocations[0];
        finalLocations[0] = finalLocations[1] = NULL;
    }

    reads[0]->setAdditionalFrontClipping(0);
    reads[1]->setAdditionalFrontClipping(0);

    return retVal;
}

    void
SimpleReadWriter::close()
{
    writer->close();
}

class SimpleReadWriterSupplier : public ReadWriterSupplier
{
public:
    SimpleReadWriterSupplier(const FileFormat* i_format, DataWriterSupplier* i_dataSupplier, const Genome* i_genome, bool i_killIfTooSlow, 
                             bool i_emitInternalScore, char *i_internalScoreTag, bool i_ignoreAlignmentAdjustmentsForOm, int i_matchReward, 
                             int i_subPenalty, int i_gapOpenPenalty, int i_gapExtendPenalty, bool i_attachAlignmentTimes)
        :
        format(i_format),
        dataSupplier(i_dataSupplier),
        genome(i_genome),
        killIfTooSlow(i_killIfTooSlow),
        emitInternalScore(i_emitInternalScore),
        ignoreAlignmentAdjustmentsForOm(i_ignoreAlignmentAdjustmentsForOm),
        matchReward(i_matchReward),
        subPenalty(i_subPenalty),
        gapOpenPenalty(i_gapOpenPenalty),
        gapExtendPenalty(i_gapExtendPenalty),
        attachAlignmentTimes(i_attachAlignmentTimes)
    {
        if (emitInternalScore) {
            if (strlen(i_internalScoreTag) != 2) {
                WriteErrorMessage("SimpleReadWriterSupplier: bad internal score tag\n");
                soft_exit(1);
            }
            strcpy(internalScoreTag, i_internalScoreTag);
        } else {
            internalScoreTag[0] = '\0';
        }
    }

    ~SimpleReadWriterSupplier()
    {
        delete dataSupplier;
    }

    virtual ReadWriter* getWriter()
    {
        return new SimpleReadWriter(format, dataSupplier->getWriter(), genome, killIfTooSlow, emitInternalScore, internalScoreTag, 
                                    ignoreAlignmentAdjustmentsForOm, matchReward, subPenalty, gapOpenPenalty, gapExtendPenalty, attachAlignmentTimes);
    }

    virtual void close()
    {
        dataSupplier->close();
    }

private:
    const FileFormat* format;
    DataWriterSupplier* dataSupplier;
    const Genome* genome;
    bool killIfTooSlow;
    bool emitInternalScore;
    char internalScoreTag[3];
    bool ignoreAlignmentAdjustmentsForOm;
    int matchReward;
    int subPenalty;
    int gapOpenPenalty;
    int gapExtendPenalty;
    bool attachAlignmentTimes;
};

    ReadWriterSupplier*
ReadWriterSupplier::create(
    const FileFormat* format,
    DataWriterSupplier* dataSupplier,
    const Genome* genome,
    bool killIfTooSlow,
    bool emitInternalScore, 
    char *internalScoreTag,
    bool ignoreAlignmentAdjustmentsForOm,
    int matchReward,
    int subPenalty,
    int gapOpenPenalty,
    int gapExtendPenalty,
    bool attachAlignmentTimes)
{
    return new SimpleReadWriterSupplier(format, dataSupplier, genome, killIfTooSlow, emitInternalScore, internalScoreTag, 
                                        ignoreAlignmentAdjustmentsForOm, matchReward, subPenalty, gapOpenPenalty, 
                                        gapExtendPenalty, attachAlignmentTimes);
}

