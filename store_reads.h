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

#ifndef APP_CUDAMAPPER_STORE_READS_H_
#define APP_CUDAMAPPER_STORE_READS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadsConfig
// ----------------------------------------------------------------------------

template <typename TUseReadStore_       = True,
          typename TUseReadNameStore_   = True,
          typename TForward_            = True,
          typename TReverse_            = True,
          typename TFragStoreConfig_    = void>
struct ReadsConfig
{
    typedef TUseReadStore_      TUseReadStore;
    typedef TUseReadNameStore_  TUseReadNameStore;
    typedef TForward_           TForward;
    typedef TReverse_           TReverse;
    typedef TFragStoreConfig_   TFragStoreConfig;
};

// ----------------------------------------------------------------------------
// Class Reads
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = ReadsConfig<> >
struct Reads
{
    typedef typename TConfig::TFragStoreConfig      TFragStoreConfig_;
    typedef FragmentStore<void, TFragStoreConfig_>  TFragmentStore_;
    typedef typename TFragmentStore_::TReadSeqStore TReadSeqStore;

    Holder<TFragmentStore_>         _store;
    unsigned                        readsCount;

    Reads() :
        _store(),
        readsCount(0)
    {}

    template <typename TFragmentStore>
    Reads(TFragmentStore & store) :
        _store(store),
        readsCount(0)
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
    typedef Reads<TSpec, TConfig>                   TReads;

    TStream                         _file;
    AutoSeqStreamFormat             _fileFormat;
    std::auto_ptr<TRecordReader>    _reader;
    Holder<TReads>                  reads;

    ReadsLoader() {}

    ReadsLoader(TReads & reads) :
        reads(reads)
    {}
};

// ----------------------------------------------------------------------------
// Class ReadsLoader; PaiedEnd
// ----------------------------------------------------------------------------

template <typename TConfig>
struct ReadsLoader<PairedEnd, TConfig>
{
    typedef std::fstream                            TStream;
    typedef RecordReader<TStream, SinglePass<> >    TRecordReader;
    typedef Reads<PairedEnd, TConfig>               TReads;

    TStream                             _file1;
    TStream                             _file2;
    Pair<AutoSeqStreamFormat>           _fileFormat;
    Pair<std::auto_ptr<TRecordReader> > _reader;
    Holder<TReads>                      reads;

    ReadsLoader(TReads & reads) :
        reads(reads)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ReadsHost<T>::Type                                    [TObject]
// ----------------------------------------------------------------------------

template <typename TObject>
struct ReadsHost {};

template <typename TObject>
struct ReadsHost<TObject const> :
    ReadsHost<TObject> {};

// ----------------------------------------------------------------------------
// Metafunction ReadsHost<T>::Type                                [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
struct ReadsHost<ReadsLoader<TSpec, TConfig> >
{
    typedef Reads<TSpec, TConfig>   Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setReads()                                                [TObject]
// ----------------------------------------------------------------------------

template <typename TObject, typename TReads>
inline void
setReads(TObject & object, TReads /* const */ & reads)
{
    setValue(object.reads, reads);
}

// ----------------------------------------------------------------------------
// Function getReads()                                                [TObject]
// ----------------------------------------------------------------------------

template <typename TObject>
inline typename ReadsHost<TObject>::Type &
getReads(TObject & object)
{
    return value(object.reads);
}

// ----------------------------------------------------------------------------
// Function clear()                                                     [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void clear(Reads<TSpec, TConfig> & reads)
{
    clearReads(value(reads._store));
    reads.readsCount = 0;
}

// ----------------------------------------------------------------------------
// Function appendSeq()                                                 [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeq>
inline void appendSeq(Reads<TSpec, TConfig> & reads, TReadSeq const & seq)
{
    appendValue(getSeqs(reads), seq, Generous());
}

// ----------------------------------------------------------------------------
// Function appendName()                                                [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadName>
inline void appendName(Reads<TSpec, TConfig> & /* reads */, TReadName const & /* seqName */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TReadName>
inline void appendName(Reads<TSpec, TConfig> & reads, TReadName const & seqName, True const & /* tag */)
{
    appendValue(getNames(reads), seqName, Generous());
}

// ----------------------------------------------------------------------------
// Function appendId()                                                  [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadId>
inline void appendId(Reads<TSpec, TConfig> & /* reads */, TReadId const & /* matePairId */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TReadId>
inline void appendId(Reads<TSpec, TConfig> & reads, TReadId const & matePairId, True const & /* tag */)
{
    typedef FragmentStore<TSpec, typename TConfig::TFragStoreConfig>    TFragmentStore;
    typedef typename Value<typename TFragmentStore::TReadStore>::Type   TReadStoreElement;

	TReadStoreElement r;
	r.matePairId = matePairId;

	appendValue(getIds(reads), r, Generous());
}

// ----------------------------------------------------------------------------
// Function getSeqs()                                                   [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline typename FragmentStore<TSpec, typename TConfig::TFragStoreConfig>::TReadSeqStore &
getSeqs(Reads<TSpec, TConfig> const & reads)
{
    return value(reads._store).readSeqStore;
}

// ----------------------------------------------------------------------------
// Function getNames()                                                  [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline typename FragmentStore<TSpec, typename TConfig::TFragStoreConfig>::TReadNameStore &
getNames(Reads<TSpec, TConfig> const & reads)
{
    return value(reads._store).readNameStore;
}

// ----------------------------------------------------------------------------
// Function getIds()                                                    [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline typename FragmentStore<TSpec, typename TConfig::TFragStoreConfig>::TReadStore &
getIds(Reads<TSpec, TConfig> const & reads)
{
    return value(reads._store).readStore;
}

// ----------------------------------------------------------------------------
// Function getReadId()                                                 [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSeqId>
inline TSeqId getReadId(Reads<TSpec, TConfig> const & reads, TSeqId seqId)
{
    // Deal with reverse complemented reads.
    return isForward(reads, seqId) ? seqId : seqId - reads.readsCount;
}

// ----------------------------------------------------------------------------
// Function isForward()                                                 [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSeqId>
inline bool isForward(Reads<TSpec, TConfig> const & reads, TSeqId seqId)
{
    return seqId < reads.readsCount;
}

// ----------------------------------------------------------------------------
// Function isReverse()                                                 [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSeqId>
inline bool isReverse(Reads<TSpec, TConfig> const & reads, TSeqId seqId)
{
    return seqId >= reads.readsCount;
}

// ----------------------------------------------------------------------------
// Function open()                                                [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TString>
void open(ReadsLoader<TSpec, TConfig> & loader, TString const & readsFile)
{
    typedef ReadsLoader<TSpec, TConfig>             TReadsLoader;
    typedef typename TReadsLoader::TRecordReader    TRecordReader;

    // Open file.
    loader._file.open(toCString(readsFile), std::ios::binary | std::ios::in);

    if (!loader._file.is_open())
        throw std::runtime_error("Error while opening reads file.");

    // Initialize record reader.
    loader._reader.reset(new TRecordReader(loader._file));

    // Autodetect file format.
    if (!guessStreamFormat(*(loader._reader), loader._fileFormat))
        throw std::runtime_error("Error while guessing reads file format.");
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
        throw std::runtime_error("Error while opening reads file.");

    // Initialize record reader.
    loader._reader.i1.reset(new TRecordReader(loader._file1));
    loader._reader.i2.reset(new TRecordReader(loader._file2));

    // Autodetect file format.
    if (!guessStreamFormat(*(loader._reader.i1), loader._fileFormat.i1) ||
        !guessStreamFormat(*(loader._reader.i2), loader._fileFormat.i2))
        throw std::runtime_error("Error while guessing reads file format.");
}

// ----------------------------------------------------------------------------
// Function close()                                               [ReadsLoader]
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
// Function load()                                                [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize, typename TReader, typename TFormat>
void _load(ReadsLoader<TSpec, TConfig> & loader, TSize count, TReader & reader, TFormat & format)
{
    typedef FragmentStore<TSpec, typename TConfig::TFragStoreConfig>    TFragmentStore;
    typedef typename Value<typename TFragmentStore::TReadStore>::Type   TReadStoreElement;
    typedef typename TConfig::TFragStoreConfig::TReadSeq                TReadSeq;

    CharString  seqName;
    TReadSeq    seq;

    // Read records.
    for (; count > 0 && !atEnd(reader); count--)
    {
        if (readRecord(seqName, seq, reader, format) != 0)
            throw std::runtime_error("Error while reading read record.");

        appendSeq(getReads(loader), seq);
        appendName(getReads(loader), seqName, typename TConfig::TUseReadNameStore());
        appendId(getReads(loader), TReadStoreElement::INVALID_ID, typename TConfig::TUseReadStore());
    }
}

template <typename TSpec, typename TConfig, typename TSize>
void load(ReadsLoader<TSpec, TConfig> & loader, TSize count)
{
    _load(loader, count, *(loader._reader), loader._fileFormat);

    getReads(loader).readsCount = length(getSeqs(getReads(loader)));

    // Append reverse complemented reads.
    _loadReverseComplement(loader);
}

template <typename TConfig, typename TSize>
void load(ReadsLoader<PairedEnd, TConfig> & loader, TSize count)
{
    _load(loader, count, *(loader._reader.i1), loader._fileFormat.i1);
    _load(loader, count, *(loader._reader.i2), loader._fileFormat.i2);

    getReads(loader).readsCount = length(getSeqs(getReads(loader)));

    // Append reverse complemented reads.
    _loadReverseComplement(loader);
}

// ----------------------------------------------------------------------------
// Function _loadReverseComplement()                              [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void _loadReverseComplement(ReadsLoader<TSpec, TConfig> & loader)
{
    typedef FragmentStore<TSpec, typename TConfig::TFragStoreConfig>    TFragmentStore;
    typedef typename Size<typename TFragmentStore::TReadSeqStore>::Type TReadSeqStoreSize;
    typedef typename TConfig::TFragStoreConfig::TReadSeq                TReadSeq;

    for (TReadSeqStoreSize readId = 0; readId < getReads(loader).readsCount; ++readId)
    {
        TReadSeq const & read = getSeqs(getReads(loader))[readId];
        appendSeq(getReads(loader), read);
        reverseComplement(back(getSeqs(getReads(loader))));
    }
}

// ----------------------------------------------------------------------------
// Function atEnd()                                               [ReadsLoader]
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
    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
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
