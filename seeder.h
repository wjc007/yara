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

#ifndef APP_CUDAMAPPER_SEEDER_H_
#define APP_CUDAMAPPER_SEEDER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SeedsCounter
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec = void>
struct SeedsCounter
{
    String<TSize> seedsPerRead;

    template <typename TCount>
    SeedsCounter(TCount readsCount)
    {
        resize(seedsPerRead, readsCount, 0, Exact());
    }

    template <typename TPos, typename TLength>
    void operator() (TPos pos, TLength /* length */)
    {
        seedsPerRead[getSeqNo(pos)]++;
    }
};

// ----------------------------------------------------------------------------
// Class SeedsManager
// ----------------------------------------------------------------------------

template <typename TSeeds, typename TSeedsCount, typename TSpec = void>
struct SeedsManager
{
    TSeeds & seeds;
    TSeedsCount & seedsPerRead;

    SeedsManager(TSeeds & seeds, TSeedsCount & seedsPerRead) :
        seeds(seeds),
        seedsPerRead(seedsPerRead)
    {
        std::partial_sum(begin(seedsPerRead, Standard()), end(seedsPerRead, Standard()), begin(seedsPerRead, Standard()));

        // Resize space for seeds.
        clear(seeds);
        resize(seeds, back(seedsPerRead), Exact());
    }

    template <typename TPos, typename TLength>
    void operator() (TPos pos, TLength len)
    {
        --seedsPerRead[getSeqNo(pos)];
        assignInfixWithLength(seeds, seedsPerRead[getSeqNo(pos)], pos, len);
    }

    ~SeedsManager()
    {
        _refreshStringSetLimits(seeds);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function selectSeeds()
// ----------------------------------------------------------------------------

template <typename TReadSeqs, typename TReadSeqId, typename TDelegate>
inline void selectSeeds(TReadSeqs const & readSeqs, TReadSeqId readSeqId, TDelegate & delegate)
{
    typedef typename StringSetPosition<TReadSeqs>::Type     TPos;
    typedef typename Value<TReadSeqs>::Type                 TReadSeq;
    typedef typename Size<TReadSeq>::Type                   TSize;

    TSize readLength = length(readSeqs[readSeqId]);
//    TSize errorsPerRead = std::ceil(readsLength * (options.errorRate / 100.0));
    TSize seedsPerRead = 6;//errorsPerRead + 1;
    TSize seedsLength = readLength / seedsPerRead;

    for (TSize seedId = 0; seedId < seedsPerRead; ++seedId)
        delegate(TPos(readSeqId, seedId * seedsLength), seedsLength);
}

template <typename TReadSeqs, typename TSpec>
inline void selectSeeds(StringSet<TReadSeqs, Segment<TSpec> > & seeds, TReadSeqs & readSeqs)
{
    typedef StringSet<TReadSeqs, Segment<TSpec> >       TSeeds;
    typedef typename Value<TReadSeqs>::Type             TReadSeq;
    typedef typename Id<TReadSeqs>::Type                TId;
    typedef typename Size<TReadSeq>::Type               TSize;
    typedef SeedsCounter<TSize>                         TCounter;
    typedef SeedsManager<TSeeds, String<TSize> >        TManager;

    TId readsCount = length(readSeqs);

    TCounter counter(readsCount);
    for (TId readSeqId = 0; readSeqId < readsCount; ++readSeqId)
        selectSeeds(readSeqs, readSeqId, counter);

    TManager manager(seeds, counter.seedsPerRead);
    for (TId readSeqId = 0; readSeqId < readsCount; ++readSeqId)
        selectSeeds(readSeqs, readSeqId, manager);
}

// --------------------------------------------------------------------------
// Function getHostPos()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TPos>
inline typename StringSetPosition<THost>::Type
getHostPos(StringSet<THost, Segment<TSpec> > const & me, TPos pos)
{
    return me.positions[pos];
}

// ----------------------------------------------------------------------------
// Function getSeedIds()
// ----------------------------------------------------------------------------
// getSeqNosByHostSeqNo(seqNo)

template <typename THost, typename TSpec, typename TReadId>
inline Pair<typename Id<StringSet<THost, Segment<TSpec> > const>::Type>
getSeedIds(StringSet<THost, Segment<TSpec> > const & seeds, TReadId readId)
{
    typedef StringSet<THost, Segment<TSpec> > const         TStringSet;
    typedef typename Id<TStringSet>::Type                   TId;
    typedef typename StringSetPositions<THost>::Type const  TPositions;
    typedef typename Iterator<TPositions, Standard>::Type   TPositionsIter;
    typedef typename Value<TPositions>::Type                TPos;

    TPositionsIter posBegin = begin(seeds.positions, Standard());
    TPositionsIter posEnd = end(seeds.positions, Standard());
    TPos pos(readId, 0);

    TPositionsIter seedsBegin = std::lower_bound(posBegin, posEnd, TPos(readId, 0));
    TPositionsIter seedsEnd = std::lower_bound(posBegin, posEnd, TPos(readId + 1, 0));

    return Pair<TId>(position(seedsBegin, seeds.positions), position(seedsEnd, seeds.positions));
}

// ----------------------------------------------------------------------------
// Function getReadSeqId()
// ----------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TSeedId>
inline typename Id<THost>::Type
getReadSeqId(StringSet<THost, Segment<TSpec> > const & seeds, TSeedId seedId)
{
    return getSeqNo(getHostPos(seeds, seedId));
}

// ----------------------------------------------------------------------------
// Function getPosInRead()
// ----------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TSeedId>
inline Pair<typename Position<typename Value<THost>::Type>::Type>
getPosInRead(StringSet<THost, Segment<TSpec> > const & seeds, TSeedId seedId)
{
    typedef typename Position<typename Value<THost>::Type>::Type TPos;

    TPos seedBegin = getSeqOffset(getHostPos(seeds, seedId));
    TPos seedEnd = seedBegin + length(value(seeds, seedId));

    return Pair<TPos>(seedBegin, seedEnd);
}

#endif  // #ifndef APP_CUDAMAPPER_SEEDER_H_
