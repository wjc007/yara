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

#ifndef SEQAN_EXTRAS_MASAI_READS_H_
#define SEQAN_EXTRAS_MASAI_READS_H_

#include <stdexcept>

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
    unsigned                        _avgSeqLengthEstimate;
    unsigned                        _avgNameLengthEstimate;
    unsigned                        _countEstimate;
    unsigned                        readsCount;

    Reads() :
        _store(),
        _avgSeqLengthEstimate(0),
        _avgNameLengthEstimate(0),
        _countEstimate(0),
        readsCount(0)
    {}

    template <typename TFragmentStore>
    Reads(TFragmentStore & store) :
        _store(store),
        _avgSeqLengthEstimate(0),
        _avgNameLengthEstimate(0),
        _countEstimate(0),
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
    unsigned long                   _fileSize;
    AutoSeqStreamFormat             _fileFormat;
    std::auto_ptr<TRecordReader>    _reader;
    Holder<TReads>                  reads;

    ReadsLoader() :
        _fileSize(0)
    {}

    ReadsLoader(TReads & reads) :
        _fileSize(0),
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
// Function reserve()                                                   [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void reserve(Reads<TSpec, TConfig> & reads)
{
    reserve(reads, reads._countEstimate);
}

template <typename TConfig>
inline void reserve(Reads<PairedEnd, TConfig> & /* reads */)
{}

template <typename TSpec, typename TConfig, typename TSize>
inline void reserve(Reads<TSpec, TConfig> & reads, TSize count)
{
    // Reserve space in the readSeqStore, also considering reverse complemented reads.
    reserve(value(reads._store).readSeqStore.concat, 2 * count * reads._avgSeqLengthEstimate, Exact());
    reserve(value(reads._store).readSeqStore, 2 * count, Exact());

    // Reserve space for ids.
    reserveIds(reads, count, typename TConfig::TUseReadStore());

    // Reserve space for names.
    reserveNames(reads, count, reads._avgNameLengthEstimate, typename TConfig::TUseReadNameStore());
}

// ----------------------------------------------------------------------------
// Function reserveIds()                                                [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize>
inline void reserveIds(Reads<TSpec, TConfig> & /* reads */, TSize /* space */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TSize>
inline void reserveIds(Reads<TSpec, TConfig> & reads, TSize space, True const & /* tag */)
{
    reserve(getSeqs(reads), space, Exact());
}

// ----------------------------------------------------------------------------
// Function reserveNames()                                              [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize, typename TLength>
inline void reserveNames(Reads<TSpec, TConfig> & /* reads */, TSize /* count */, TLength /* length */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TSize, typename TLength>
inline void reserveNames(Reads<TSpec, TConfig> & reads, TSize count, TLength length, True const & /* tag */)
{
    reserve(getNames(reads).concat, count * length, Exact());
    reserve(getNames(reads), count, Exact());
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
// Function avgSeqLength()                                              [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline typename Size<typename FragmentStore<TSpec, typename TConfig::TFragStoreConfig>::TReadSeqStore>::Type
avgSeqLength(Reads<TSpec, TConfig> const & reads)
{
    return reads._avgSeqLengthEstimate;
}

// ----------------------------------------------------------------------------
// Function _estimateRecordSize()                                 [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
unsigned long _estimateRecordSize(ReadsLoader<TSpec, TConfig> const & loader, Fastq const & /* tag */)
{
    // 6 stands for: @, +, and four \n.
    return getReads(loader)._avgNameLengthEstimate + 2 * getReads(loader)._avgSeqLengthEstimate + 6;
}

template <typename TSpec, typename TConfig>
unsigned long _estimateRecordSize(ReadsLoader<TSpec, TConfig> const & loader, Fasta const & /* tag */)
{
    // 3 stands for: >, and two \n.
    return getReads(loader)._avgNameLengthEstimate + getReads(loader)._avgSeqLengthEstimate + 3;
}

// ----------------------------------------------------------------------------
// Function _estimateReadsStatistics()                            [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void _estimateReadsStatistics(ReadsLoader<TSpec, TConfig> & loader)
{
    typedef typename TConfig::TFragStoreConfig::TReadSeq    TReadSeq;

    CharString  seqName;
    TReadSeq    seq;

    // Read first record.
    if (readRecord(seqName, seq, *(loader._reader), loader._fileFormat) != 0)
        return;

    // Estimate read seqs and names length.
    getReads(loader)._avgSeqLengthEstimate = length(seq);
    getReads(loader)._avgNameLengthEstimate = length(seqName);

    // Estimate record size.
    unsigned long recordSize;
    switch (loader._fileFormat.tagId)
    {
    case Find<AutoSeqStreamFormat, Fasta>::VALUE:
        recordSize = _estimateRecordSize(loader, Fasta());
        break;
    case Find<AutoSeqStreamFormat, Fastq>::VALUE:
        recordSize = _estimateRecordSize(loader, Fastq());
        break;
    default:
        recordSize = 0;
        break;
    }

    // Estimate number of reads in file.
    if (recordSize > 0)
        getReads(loader)._countEstimate = loader._fileSize / recordSize;
    else
        getReads(loader)._countEstimate = 0;
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

    // Compute file size.
    loader._file.seekg(0, std::ios::end);
    loader._fileSize = loader._file.tellg();
    loader._file.seekg(0, std::ios::beg);

    // Initialize record reader.
    loader._reader.reset(new TRecordReader(loader._file));

    // Autodetect file format.
    if (!guessStreamFormat(*(loader._reader), loader._fileFormat))
        throw std::runtime_error("Error while guessing reads file format.");

    // Estimate statistics for reads in file.
    _estimateReadsStatistics(loader);

    // Reopen file and reinitialize record reader, the ugly way...
    loader._file.close();
    loader._file.open(toCString(readsFile), std::ios::binary | std::ios::in);
    loader._reader.reset(new TRecordReader(loader._file));
}

template <typename TConfig, typename TString>
void open(ReadsLoader<PairedEnd, TConfig> & loader, TString const & readsLeftFile, TString const & readsRightFile)
{
    if (!loadReads(value(getReads(loader)._store), readsLeftFile, readsRightFile))
        throw std::runtime_error("Error while loading reads file.");

    getReads(loader).readsCount = length(getSeqs(getReads(loader)));

    _loadReverseComplement(loader);
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
void close(ReadsLoader<PairedEnd, TConfig> & /* loader */)
{
}

// ----------------------------------------------------------------------------
// Function load()                                                [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void load(ReadsLoader<TSpec, TConfig> & loader)
{
    load(loader, MaxValue<unsigned long>::VALUE);
}

template <typename TConfig>
void load(ReadsLoader<PairedEnd, TConfig> & /* loader */)
{
}

template <typename TSpec, typename TConfig, typename TSize>
void load(ReadsLoader<TSpec, TConfig> & loader, TSize count)
{
    switch (loader._fileFormat.tagId)
    {
    case Find<AutoSeqStreamFormat, Fasta>::VALUE:
        load(loader, count, Fasta());
        break;
    case Find<AutoSeqStreamFormat, Fastq>::VALUE:
        load(loader, count, Fastq());
        break;
    default:
        throw std::runtime_error("Unsupported reads file format.");
        break;
    }
}

template <typename TSpec, typename TConfig, typename TSize, typename TFormat>
void load(ReadsLoader<TSpec, TConfig> & loader, TSize count, TFormat const & /* tag */)
{
    typedef FragmentStore<TSpec, typename TConfig::TFragStoreConfig>    TFragmentStore;
    typedef typename Value<typename TFragmentStore::TReadStore>::Type   TReadStoreElement;
    typedef typename TConfig::TFragStoreConfig::TReadSeq                TReadSeq;

    CharString  seqName;
    TReadSeq    seq;

    // Read records.
    for (; count > 0 && !atEnd(loader); count--)
    {
        // NOTE(esiragusa): AutoFormat seems to thrash memory allocation.
//        if (readRecord(seqName, seq, *(loader._reader), loader._fileFormat) != 0)

        if (readRecord(seqName, seq, *(loader._reader), TFormat()) != 0)
            throw std::runtime_error("Error while reading read record.");

        appendSeq(getReads(loader), seq);
        appendName(getReads(loader), seqName, typename TConfig::TUseReadNameStore());
        appendId(getReads(loader), TReadStoreElement::INVALID_ID, typename TConfig::TUseReadStore());
    }

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

#endif  // #ifndef SEQAN_EXTRAS_MASAI_READS_H_
