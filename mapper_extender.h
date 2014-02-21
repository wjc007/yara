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

#ifndef APP_CUDAMAPPER_MAPPER_EXTENDER_H_
#define APP_CUDAMAPPER_MAPPER_EXTENDER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class HitsExtender
// ----------------------------------------------------------------------------
// One instance per thread.

template <typename TSpec, typename Traits>
struct HitsExtender
{
    typedef typename Traits::TContigSeqs       TContigSeqs;
    typedef typename Traits::TContigsPos       TContigsPos;
    typedef typename Traits::TReadSeqs         TReadSeqs;
    typedef typename Traits::TReadSeq          TReadSeq;
    typedef typename Traits::TReadsContext     TReadsContext;
    typedef typename Traits::TMatches          TMatches;
    typedef typename Traits::TMatch            TMatch;
    typedef typename Traits::TSeeds            TSeeds;
    typedef typename Traits::THits             THits;
    typedef typename Traits::TSA               TSA;

    typedef AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_> TMyersSpec;
    typedef Myers<TMyersSpec, True, void>               TAlgorithm;
    typedef Extender<TContigSeqs, TReadSeq, TAlgorithm> TExtender;

    // Thread-private data.
    TExtender           extender;
    TMatch              prototype;

    // Shared-memory read-write data.
    TReadsContext &     ctx;
    TMatches &          matches;

    // Shared-memory read-only data.
    TContigSeqs const & contigSeqs;
    TReadSeqs &         readSeqs;
    TSeeds const &      seeds;
    THits const &       hits;
    TSA const &         sa;
    Options const &     options;

    HitsExtender(TReadsContext & ctx,
                 TMatches & matches,
                 TContigSeqs const & contigSeqs,
                 TSeeds const & seeds,
                 THits const & hits,
                 TSA const & sa,
                 Options const & options) :
        extender(contigSeqs),
        prototype(),
        ctx(ctx),
        matches(matches),
        contigSeqs(contigSeqs),
        readSeqs(host(seeds)),
        seeds(seeds),
        hits(hits),
        sa(sa),
        options(options)
    {
        // Iterate over all hits.
        iterate(hits, *this, Rooted(), typename Traits::TThreading());
    }

    template <typename THitsIterator>
    void operator() (THitsIterator const & hitsIt)
    {
        _extendHitImpl(*this, hitsIt);
    }

    template <typename TMatchPos, typename TMatchErrors>
    void operator() (TMatchPos matchBegin, TMatchPos matchEnd, TMatchErrors matchErrors)
    {
        _addMatchImpl(*this, matchBegin, matchEnd, matchErrors);
    }

//    void operator() (TExtender const & /* extender */)
//    {
//        _addMatchImpl(*this);
//    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _extendHitImpl()
// ----------------------------------------------------------------------------
// Extends one hit.

template <typename TSpec, typename Traits, typename THitsIterator>
inline void _extendHitImpl(HitsExtender<TSpec, Traits> & me, THitsIterator const & hitsIt)
{
    typedef HitsExtender<TSpec, Traits>                 THitsExtender;

    typedef typename THitsExtender::TContigSeqs         TContigSeqs;
    typedef typename Size<TContigSeqs>::Type            TContigId;
    typedef typename THitsExtender::TContigsPos         TContigsPos;

    typedef typename THitsExtender::TReadSeqs           TReadSeqs;
    typedef typename THitsExtender::TReadSeq            TReadSeq;
    typedef typename Id<TReadSeqs>::Type                TReadId;
    typedef Pair<typename Position<TReadSeq>::Type>     TReadPos;
    typedef typename Size<TReadSeq>::Type               TReadSeqSize;
    typedef typename Size<TReadSeq>::Type               TErrors;

    typedef typename THitsExtender::TSeeds              TSeeds;
    typedef typename Id<TSeeds>::Type                   TSeedId;

    typedef typename THitsExtender::THits               THits;
    typedef typename Value<THits>::Type                 THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef typename Position<THit>::Type               THitRange;
    typedef unsigned char                               THitErrors;

    typedef typename THitsExtender::TMatches            TMatches;
    typedef typename Value<TMatches>::Type              TMatch;

    typedef typename THitsExtender::TSA                 TSA;
    typedef typename Size<TSA>::Type                    TSAPos;
    typedef typename Value<TSA>::Type                   TSAValue;

    typedef typename THitsExtender::TReadsContext       TReadsContext;

    // Get hit id.
    THitId hitId = position(hitsIt);

    // Extract hit info.
    TSeedId seedId = getSeedId(me.hits, hitId);
    THitRange hitRange = getRange(me.hits, hitId);
    THitErrors hitErrors = getErrors(me.hits, hitId);

    // Get read.
    TReadId readSeqId = getReadSeqId(me.seeds, seedId);
    TReadSeq readSeq = me.readSeqs[readSeqId];

    // Skip unseeded and mapped reads.
    if (getStatus(me.ctx, readSeqId) != STATUS_SEEDED) return;

    // Fill readSeqId.
    setReadId(me.prototype, me.readSeqs, readSeqId);

    // Get position in read.
    TReadPos readPos = getPosInRead(me.seeds, seedId);
    TReadSeqSize seedLength = getValueI2(readPos) - getValueI1(readPos);

    for (TSAPos saPos = getValueI1(hitRange); saPos < getValueI2(hitRange); ++saPos)
    {
        // Invert SA value.
        TSAValue saValue = me.sa[saPos];
        SEQAN_ASSERT_GEQ(suffixLength(saValue, me.contigSeqs), seedLength);
        if (suffixLength(saValue, me.contigSeqs) < seedLength) continue;
        setSeqOffset(saValue, suffixLength(saValue, me.contigSeqs) - seedLength);

        // Compute position in contig.
        TContigsPos contigBegin = saValue;
        TContigsPos contigEnd = posAdd(contigBegin, seedLength);

        // Get absolute number of errors.
        TErrors maxErrors = getReadErrors(me.options, length(readSeq));

        extend(me.extender,
               readSeq,
               contigBegin, contigEnd,
               getValueI1(readPos), getValueI2(readPos),
               hitErrors, maxErrors,
               me);

        // Stop when the read has been mapped.
//        if (getStatus(me.ctx, readSeqId) == STATUS_MAPPED) break;
    }

    // Full stratum analyzed.
    // TODO(esiragusa): generalize stratum for approximate seeds, where one seed can have multiple hits.
    // TODO(esiragusa): extend any fwd + any rev seed consecutively, then increase the stratum of both fwd & rev read.
//    incStratum(me.ctx, readSeqId);
}

// ----------------------------------------------------------------------------
// Function _addMatchImpl()
// ----------------------------------------------------------------------------
// Adds one match.

template <typename TSpec, typename Traits, typename TMatchPos, typename TMatchErrors>
inline void _addMatchImpl(HitsExtender<TSpec, Traits> & me,
                          TMatchPos matchBegin,
                          TMatchPos matchEnd,
                          TMatchErrors matchErrors)
{
    setContigPosition(me.prototype, matchBegin, matchEnd);
    me.prototype.errors = matchErrors;
    appendValue(me.matches, me.prototype, Insist(), typename Traits::TThreading());
}

/*
template <typename TSpec, typename Traits, typename TMatchPos, typename TMatchErrors>
inline void _addMatchImpl(HitsExtender<TSpec, Traits> & me,
                          TMatchPos matchBegin,
                          TMatchPos matchEnd,
                          TMatchErrors matchErrors)
{
    typedef typename Size<TReadSeqs>::Type   TReadSeqId;

    // TODO(esiragusa): rename prototype.readId member to prototype.readSeqId
    TReadSeqId readId = getReadId(me.readSeqs, me.prototype.readId);

//    minErrors[readId] = std::min(minErrors[readId](me.ctx, readId), matchErrors);
    setMinErrors(me.ctx, readId, matchErrors);

    // One optimal match has been reported.
    if (getMinErrors(me.ctx, readId) <= getStratum(ctx, prototype.readId))
    {
        // Mark both forward and reverse sequence as mapped.
        setStatus(me.ctx, getFirstMateFwdSeqId(me.readSeqs, readId), STATUS_MAPPED);
        setStatus(me.ctx, getFirstMateRevSeqId(me.readSeqs, readId), STATUS_MAPPED);
    }
}
*/

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_EXTENDER_H_
