// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
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

#ifndef APP_YARA_MAPPER_ALIGNER_H_
#define APP_YARA_MAPPER_ALIGNER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class MatchesAligner
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
struct MatchesAligner
{
    typedef typename Traits::TContigSeqs       TContigSeqs;
    typedef typename Traits::TReadSeqs         TReadSeqs;
    typedef typename Traits::TMatchesSet       TMatchesSet;
    typedef typename Traits::TMatches          TMatches;

    typedef String<GapAnchor<int> >            TGapAnchors;
    typedef AnchorGaps<TGapAnchors>            TAnchorGaps;

    // Thread-private data.
    TGapAnchors contigAnchors;
    TGapAnchors readAnchors;
    CharString cigar;
    CharString md;

    // Shared-memory read-write data.
//    TGaps

    // Shared-memory read-only data.
    TMatchesSet const &     matchesSet;
    TMatches const &        pairs;
    TContigSeqs const &     contigSeqs;
    TReadSeqs const &       readSeqs;
    Options const &         options;

    MatchesAligner(TMatchesSet const & matchesSet,
                   TMatches const & pairs,
                   TContigSeqs const & contigSeqs,
                   TReadSeqs const & readSeqs,
                   Options const & options) :
        matchesSet(matchesSet),
        pairs(pairs),
        contigSeqs(contigSeqs),
        readSeqs(readSeqs),
        options(options)
    {
        iterate(matchesSet, *this, Standard(), typename Traits::TThreading());
    }

    template <typename TIterator>
    void operator() (TIterator const & it)
    {
        _alignMatchImpl(*this, it);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _alignMatchImpl()
// ----------------------------------------------------------------------------
// Aligns one match.

template <typename TSpec, typename Traits, typename TMatchesIt>
inline void _alignMatchImpl(MatchesAligner<TSpec, Traits> & me, TMatchesIt const & it)
{
    typedef typename Value<TMatchesIt const>::Type  TMatches;
    typedef typename Value<TMatches>::Type          TMatch;

    typedef typename Traits::TContigSeqs            TContigSeqs;
    typedef typename Value<TContigSeqs const>::Type TContigSeq;
    typedef typename Infix<TContigSeq>::Type        TContigInfix;

    typedef typename Traits::TReadSeqs              TReadSeqs;
    typedef typename Value<TReadSeqs const>::Type   TReadSeq;

    typedef MatchesAligner<TSpec, Traits>           TMatchesAligner;
    typedef typename TMatchesAligner::TAnchorGaps   TAnchorGaps;
    typedef Gaps<TContigInfix, TAnchorGaps>         TContigGaps;
    typedef Gaps<TReadSeq, TAnchorGaps>             TReadGaps;

    TMatches const & matches = value(it);

    if (empty(matches)) return;

    // The first match is the primary one.
    TMatch const & match = front(matches);

    unsigned errors = getErrors(match);

    TReadSeq const & readSeq = me.readSeqs[getReadSeqId(match, me.readSeqs)];

    TContigInfix const & contigInfix = infix(me.contigSeqs[getContigId(match)],
                                             getContigBegin(match),
                                             getContigEnd(match));

    clear(me.contigAnchors);
    clear(me.readAnchors);
    TContigGaps contigGaps(contigInfix, me.contigAnchors);
    TReadGaps readGaps(readSeq, me.readAnchors);

    // Do not align if the match contains no gaps.
    // TODO(esiragusa): reuse DP matrix.
    if (!(errors == 0 || (errors == 1 && length(contigInfix) == length(readSeq))))
        globalAlignment(contigGaps, readGaps, Score<short, EditDistance>(), -(int)errors, (int)errors);

    clear(me.cigar);
    clear(me.md);
    getCigarString(me.cigar, contigGaps, readGaps);
//    getMDString(me.md, contigGaps, readGaps);
}

#endif  // #ifndef APP_YARA_MAPPER_ALIGNER_H_
