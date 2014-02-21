// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the Reads and ReadsLoader classes.
// ==========================================================================
// NOTE(esiragusa): the FragmentStore should provide these functionalities.

#ifndef APP_CUDAMAPPER_STORE_READS_H_
#define APP_CUDAMAPPER_STORE_READS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadsConfig
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct ReadsConfig
{
    typedef String<Dna5Q>           TReadSeq;
    typedef Owner<ConcatDirect<> >  TReadSpec;
    typedef Owner<ConcatDirect<> >  TReadNameSpec;
};

// ----------------------------------------------------------------------------
// Class Reads
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = ReadsConfig<TSpec> >
struct Reads
{
    typedef typename TConfig::TReadSeq                  TReadSeq;
    typedef typename TConfig::TReadSpec                 TReadSpec;
    typedef typename TConfig::TReadNameSpec             TReadNameSpec;

    typedef StringSet<TReadSeq, TReadSpec>              TReadSeqs;
	typedef StringSet<CharString, TReadNameSpec>        TReadNames;
	typedef NameStoreCache<TReadNames, CharString>      TReadNamesCache;

    TReadSeqs           seqs;
    TReadNames          names;
    TReadNamesCache     namesCache;

    Reads() :
        seqs(),
        names(),
		namesCache(names)
    {}
};

// ----------------------------------------------------------------------------
// Class ReadsLoader
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = ReadsConfig<> >
struct ReadsLoader
{
    typedef std::fstream                            TStream;
    typedef RecordReader<TStream, SinglePass<> >    TRecordReader;

    TStream                         _file;
    AutoSeqStreamFormat             _fileFormat;
    std::auto_ptr<TRecordReader>    _reader;
};

// ----------------------------------------------------------------------------
// Class ReadsLoader; PairedEnd
// ----------------------------------------------------------------------------

template <typename TConfig>
struct ReadsLoader<PairedEnd, TConfig>
{
    typedef std::fstream                            TStream;
    typedef RecordReader<TStream, SinglePass<> >    TRecordReader;

    TStream                             _file1;
    TStream                             _file2;
    Pair<AutoSeqStreamFormat>           _fileFormat;
    Pair<std::auto_ptr<TRecordReader> > _reader;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void clear(Reads<TSpec, TConfig> & me)
{
    clear(me.seqs);
    clear(me.names);
//    clear(me.namesCache);
}

// ----------------------------------------------------------------------------
// Function load()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize>
void load(Reads<TSpec, TConfig> & me, ReadsLoader<TSpec, TConfig> & loader, TSize count)
{
    _load(me, count, *(loader._reader), loader._fileFormat);
    appendReverseComplement(me);
}

template <typename TConfig, typename TSize>
void load(Reads<PairedEnd, TConfig> & me, ReadsLoader<PairedEnd, TConfig> & loader, TSize count)
{
    _load(me, count, *(loader._reader.i1), loader._fileFormat.i1);
    _load(me, count, *(loader._reader.i2), loader._fileFormat.i2);
    appendReverseComplement(me);
}

template <typename TSpec, typename TConfig, typename TSize, typename TReader, typename TFormat>
void _load(Reads<TSpec, TConfig> & me, TSize count, TReader & reader, TFormat & format)
{
    typedef Reads<TSpec, TConfig>           TReads;
    typedef typename TReads::TReadSeq       TReadSeq;

    CharString  seqName;
    TReadSeq    seq;

    // Read records.
    for (; count > 0 && !atEnd(reader); count--)
    {
        if (readRecord(seqName, seq, reader, format) != 0)
            throw RuntimeError("Error while reading read record.");

        appendValue(me.seqs, seq, Generous());
        appendValue(me.names, prefix(seqName, lastOf(seqName, IsSpace())), Generous());
    }
}

// ----------------------------------------------------------------------------
// Function appendReverseComplement()
// ----------------------------------------------------------------------------
// Append reverse complemented reads.

template <typename TSpec, typename TConfig>
void appendReverseComplement(Reads<TSpec, TConfig> & me)
{
    typedef Reads<TSpec, TConfig>           TReads;
    typedef typename TReads::TReadSeqs      TReadSeqs;
    typedef typename Value<TReadSeqs>::Type TReadSeq;
    typedef typename Size<TReadSeqs>::Type  TReadSeqId;

    TReadSeqId readSeqsCount = length(me.seqs);

    reserve(me.seqs, 2 * readSeqsCount, Exact());

    for (TReadSeqId readSeqId = 0; readSeqId < readSeqsCount; ++readSeqId)
    {
        TReadSeq const & read = me.seqs[readSeqId];
        appendValue(me.seqs, read, Exact());
        reverseComplement(back(me.seqs));
    }
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TString>
void open(ReadsLoader<TSpec, TConfig> & loader, TString const & readsFile)
{
    typedef ReadsLoader<TSpec, TConfig>             TReadsLoader;
    typedef typename TReadsLoader::TRecordReader    TRecordReader;

    // Open file.
    loader._file.open(toCString(readsFile), std::ios::binary | std::ios::in);

    if (!loader._file.is_open())
        throw RuntimeError("Error while opening reads file.");

    // Initialize record reader.
    loader._reader.reset(new TRecordReader(loader._file));

    // Autodetect file format.
    if (!guessStreamFormat(*(loader._reader), loader._fileFormat))
        throw RuntimeError("Error while guessing reads file format.");
}

template <typename TConfig, typename TString>
void open(ReadsLoader<PairedEnd, TConfig> & loader, Pair<TString> const & readsFile)
{
    typedef ReadsLoader<PairedEnd, TConfig>         TReadsLoader;
    typedef typename TReadsLoader::TRecordReader    TRecordReader;

    // Open files.
    loader._file1.open(toCString(readsFile.i1), std::ios::binary | std::ios::in);
    loader._file2.open(toCString(readsFile.i2), std::ios::binary | std::ios::in);

    if (!loader._file1.is_open() || !loader._file2.is_open())
        throw RuntimeError("Error while opening reads file.");

    // Initialize record reader.
    loader._reader.i1.reset(new TRecordReader(loader._file1));
    loader._reader.i2.reset(new TRecordReader(loader._file2));

    // Autodetect file format.
    if (!guessStreamFormat(*(loader._reader.i1), loader._fileFormat.i1) ||
        !guessStreamFormat(*(loader._reader.i2), loader._fileFormat.i2))
        throw RuntimeError("Error while guessing reads file format.");
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void close(ReadsLoader<TSpec, TConfig> & loader)
{
    loader._file.close();
}

template <typename TConfig>
void close(ReadsLoader<PairedEnd, TConfig> & loader)
{
    loader._file1.close();
    loader._file2.close();
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline bool atEnd(ReadsLoader<TSpec, TConfig> & reads)
{
    return atEnd(*(reads._reader));
}

template <typename TConfig>
inline bool atEnd(ReadsLoader<PairedEnd, TConfig> & reads)
{
    return atEnd(*(reads._reader.i1)) && atEnd(*(reads._reader.i2));
}

// ----------------------------------------------------------------------------
// Functions on ReadSeqsStore
// ----------------------------------------------------------------------------

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getReadSeqsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs);
}

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getReadsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs) / 2;
}

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getPairsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs) / 4;
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isFwdReadSeq(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return readSeqId < getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isRevReadSeq(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return !isFwdReadSeq(readSeqs, readSeqId);
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isFirstMate(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return isFwdReadSeq(readSeqs, readSeqId) ? readSeqId < getPairsCount(readSeqs) : readSeqId < getPairsCount(readSeqs) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isSecondMate(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return !isFirstMate(readSeqs, readSeqId);
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getFirstMateFwdSeqId(TReadSeqs const & /* readSeqs */, TPairId pairId)
{
    return pairId;
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getSecondMateFwdSeqId(TReadSeqs const & readSeqs, TPairId pairId)
{
    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    return pairId + getPairsCount(readSeqs);
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getFirstMateRevSeqId(TReadSeqs const & readSeqs, TPairId pairId)
{
//    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    return getFirstMateFwdSeqId(readSeqs, pairId) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getSecondMateRevSeqId(TReadSeqs const & readSeqs, TPairId pairId)
{
    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    return getSecondMateFwdSeqId(readSeqs, pairId) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline typename Size<TReadSeqs>::Type
getReadId(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    return isFwdReadSeq(readSeqs, readSeqId) ? readSeqId : readSeqId - getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline typename Size<TReadSeqs>::Type
getPairId(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));

    typename Size<TReadSeqs>::Type pairId = readSeqId;

    if (isRevReadSeq(readSeqs, readSeqId))
        pairId -= getReadsCount(readSeqs);

    if (isSecondMate(readSeqs, readSeqId))
        pairId -= getPairsCount(readSeqs);

    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    
    return pairId;
}

template <typename TReadSeqs, typename TReadSeqId>
inline typename Size<TReadSeqs>::Type
getMateSeqId(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));

    typename Size<TReadSeqs>::Type pairId = getPairId(readSeqs, readSeqId);

    if (isFirstMate(readSeqs, readSeqId))
    {
        if (isFwdReadSeq(readSeqs, readSeqId))
            return getSecondMateRevSeqId(readSeqs, pairId);
        else
            return getSecondMateFwdSeqId(readSeqs, pairId);
    }
    else
    {
        if (isFwdReadSeq(readSeqs, readSeqId))
            return getFirstMateRevSeqId(readSeqs, pairId);
        else
            return getFirstMateFwdSeqId(readSeqs, pairId);
    }
}

#endif  // #ifndef APP_CUDAMAPPER_STORE_READS_H_
