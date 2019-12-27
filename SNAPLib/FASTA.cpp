/*++

Module Name:

    FASTA.cpp

Abstract:

    FASTA reader

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#include "stdafx.h"
#include "Compat.h"
#include "FASTA.h"
#include "Error.h"
#include "exit.h"
#include "Util.h"

using namespace std;

void
MarkALTContigIfAppropriate(
	const char		*contigName,
	GenomeDistance	 contigSize,
	const char		*opt_in_alt_names,
	int				 opt_in_alt_names_count,
	const char		*opt_out_alt_names,
	int				 opt_out_alt_names_count,
	GenomeDistance	 maxSizeForAutomaticALT,
	Genome			*genome)
{
	for (int i = 0; i < opt_out_alt_names_count; i++) {
		if (!stricmp(opt_out_alt_names+i, contigName)) {
			return;
		}
	} // opt out

	if (contigSize <= maxSizeForAutomaticALT) {
		genome->markContigALT(contigName);
		return;
	}

	for (int i = 0; i < opt_in_alt_names_count; i++) {
		if (!stricmp(opt_in_alt_names + i, contigName)) {
			genome->markContigALT(contigName);
			return;
		} // match
	} // opt in
} // MarkALTContigIfAppropriate


//
// There are several ways of specifying ALT contigs.  There is an opt-in list of ALTs, an opt-out list of regular chromosomes (these must be mutually
// exclusive), and a size cutoff below which is contig is an ALT.  The opt-in and opt-out lists supersede the size cutoff.
//

    const Genome *
ReadFASTAGenome(
	const char		*fileName,
	const char		*pieceNameTerminatorCharacters,
	bool			 spaceIsAPieceNameTerminator,
	unsigned		 chromosomePaddingSize,
	const char		*opt_in_alt_names,
	int				 opt_in_alt_names_count,
	const char		*opt_out_alt_names,
	int				 opt_out_alt_names_count,
	GenomeDistance	 maxSizeForAutomaticALT)
{
    //
    // We need to know a bound on the size of the genome before we create the Genome object.
    // A bound is the number of bytes in the FASTA file, because we store at most one base per
    // byte.  Get the file size to use for this bound.
    //
    _int64 fileSize = QueryFileSize(fileName);
    bool isValidGenomeCharacter[256];

    for (int i = 0; i < 256; i++) {
        isValidGenomeCharacter[i] = false;
    }

    isValidGenomeCharacter['A'] = isValidGenomeCharacter['T'] = isValidGenomeCharacter['C'] = isValidGenomeCharacter['G'] = isValidGenomeCharacter['N'] = true;
    isValidGenomeCharacter['a'] = isValidGenomeCharacter['t'] = isValidGenomeCharacter['c'] = isValidGenomeCharacter['g'] = isValidGenomeCharacter['n'] = true;

    FILE *fastaFile = fopen(fileName, "r");
    if (fastaFile == NULL) {
        WriteErrorMessage("Unable to open FASTA file '%s' (even though we already got its size)\n",fileName);
        return NULL;
    }

    int lineBufferSize = 0;
    char *lineBuffer;
 
    //
    // Count the chromosomes
    //
    unsigned nChromosomes = 0;

    while (NULL != reallocatingFgets(&lineBuffer,&lineBufferSize,fastaFile)) {
        if (lineBuffer[0] == '>') {
            nChromosomes++;
        }
    }
    rewind(fastaFile);

    Genome *genome = new Genome(fileSize + (nChromosomes+1) * (size_t)chromosomePaddingSize, fileSize + (nChromosomes+1) * (size_t)chromosomePaddingSize, chromosomePaddingSize, nChromosomes + 1);

    char *paddingBuffer = new char[chromosomePaddingSize+1];
    for (unsigned i = 0; i < chromosomePaddingSize; i++) {
        paddingBuffer[i] = 'n';
    }
    paddingBuffer[chromosomePaddingSize] = '\0';

    bool warningIssued = false;
    bool inAContig = false;
	GenomeDistance contigLength = 0;
	char *lastContigName = new char[lineBufferSize];

    while (NULL != reallocatingFgets(&lineBuffer, &lineBufferSize, fastaFile)) {
        if (lineBuffer[0] == '>') {
			if (inAContig) {
				MarkALTContigIfAppropriate(lastContigName, contigLength, opt_in_alt_names, opt_in_alt_names_count, opt_out_alt_names,
					opt_out_alt_names_count, maxSizeForAutomaticALT, genome);
			}

            inAContig = true;
			contigLength = 0;

            //
            // A new contig.  Add in the padding first.
            //
            genome->addData(paddingBuffer);

            //
            // Now supply the chromosome name.
            //
            if (NULL != pieceNameTerminatorCharacters) {
                for (int i = 0; i < strlen(pieceNameTerminatorCharacters); i++) {
                    char *terminator = strchr(lineBuffer+1, pieceNameTerminatorCharacters[i]);
                    if (NULL != terminator) {
                        *terminator = '\0';
                    }
                }
            }

            if (spaceIsAPieceNameTerminator) {
                char *terminator = strchr(lineBuffer, ' ');
                if (NULL != terminator) {
                    *terminator = '\0';
                }
                terminator = strchr(lineBuffer, '\t');
                if (NULL != terminator) {
                    *terminator = '\0';
                }
            }

            char *terminator = strchr(lineBuffer, '\n');
            if (NULL != terminator) {
                *terminator = '\0';
            }

            terminator = strchr(lineBuffer, '\r');
            if (NULL != terminator) {
                *terminator = '\0';
            }

            genome->startContig(lineBuffer+1);
			strcpy(lastContigName, lineBuffer + 1);
        } else {
            if (!inAContig) {
                WriteErrorMessage("\nFASTA file doesn't beging with a contig name (i.e., the first line doesn't start with '>').\n");
                soft_exit(1);
            }

            //
            // Convert it to upper case and truncate the newline before adding it to the genome.
            //

            char *newline = strchr(lineBuffer, '\n');
            if (NULL != newline) {
                *newline = 0;
            }

            size_t lineLen = strlen(lineBuffer);

			for (unsigned i = 0; i < lineLen; i++) {
              lineBuffer[i] = toupper(lineBuffer[i]);
            }

			contigLength += lineLen;

			for (unsigned i = 0; i < lineLen; i++) {
                if (!isValidGenomeCharacter[(unsigned char)lineBuffer[i]]) {
                    if (!warningIssued) {
                        WriteErrorMessage("\nFASTA file contained a character that's not a valid base (or N): '%c', full line '%s'; \nconverting to 'N'.  This may happen again, but there will be no more warnings.\n", lineBuffer[i], lineBuffer);
                        warningIssued = true;
                    }
                    lineBuffer[i] = 'N';
                }
            }
            genome->addData(lineBuffer);
        }
    }

	if (!inAContig) {
		WriteErrorMessage("The FASTA file was empty.");
		return NULL;
	}

	MarkALTContigIfAppropriate(lastContigName, contigLength, opt_in_alt_names, opt_in_alt_names_count, opt_out_alt_names,
		opt_out_alt_names_count, maxSizeForAutomaticALT, genome);

    //
    // And finally add padding at the end of the genome.
    //
    genome->addData(paddingBuffer);
    genome->fillInContigLengths();
    genome->sortContigsByName();

    fclose(fastaFile);
    delete [] paddingBuffer;
    delete [] lineBuffer;
	delete [] lastContigName;
    return genome;
}

//
// TODO: Reduce code duplication with the mutator.
//
bool AppendFASTAGenome(const Genome *genome, FILE *fasta, const char *prefix="")
{
    int nContigs = genome->getNumContigs();
    const Genome::Contig *contigs = genome->getContigs();
    for (int i = 0; i < nContigs; ++i) {
        const Genome::Contig &contig = contigs[i];
        GenomeLocation start = contig.beginningLocation;
        GenomeLocation end = i + 1 < nContigs ? contigs[i + 1].beginningLocation : genome->getCountOfBases();
        GenomeDistance size = end - start;
        const char *bases = genome->getSubstring(start, size);

        fprintf(fasta, ">%s%s\n", prefix, contig.name);
        fwrite(bases, 1, size, fasta);
        fputc('\n', fasta);
    }
    return !ferror(fasta);
}
