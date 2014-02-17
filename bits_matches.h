// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
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
// This file contains classes for storing and manipulating matches.
// ==========================================================================

#ifndef APP_CUDAMAPPER_BITS_MATCHES_H_
#define APP_CUDAMAPPER_BITS_MATCHES_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Match
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Match
{
    unsigned        readId;
    unsigned        contigBegin;
    unsigned        contigEnd;
    unsigned char   contigId;
    unsigned char   errors;
};

template <typename TSpec>
inline unsigned getReadId(Match<TSpec> const & me)
{
    return me.readId;
}

template <typename TSpec>
inline unsigned getContigId(Match<TSpec> const & me)
{
    return me.contigId;
}

template <typename TSpec>
inline unsigned getContigBegin(Match<TSpec> const & me)
{
    return me.contigBegin;
}

template <typename TSpec>
inline unsigned getContigEnd(Match<TSpec> const & me)
{
    return me.contigEnd;
}

template <typename TSpec>
inline bool onForwardStrand(Match<TSpec> const & me)
{
    return true;
}

template <typename TSpec>
inline bool onReverseStrand(Match<TSpec> const & me)
{
    return !onForwardStrand(me);
}

template <typename TSpec>
inline unsigned char getScore(Match<TSpec> const & /* me */)
{
    return 0;
}
template <typename TSpec>
inline unsigned char getErrors(Match<TSpec> const & me)
{
    return me.errors;
}

// ----------------------------------------------------------------------------
// Class MatchSorterByXXX
// ----------------------------------------------------------------------------

template <typename TMatch>
struct MatchSorterByReadId
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return a.readId <= b.readId;
    }
};

template <typename TMatch>
struct MatchSorterByBeginPos
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return (a.contigId < b.contigId) || (a.contigId == b.contigId && a.contigBegin < b.contigBegin);
    }
};

template <typename TMatch>
struct MatchSorterByEndPos
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return (a.contigId < b.contigId) || (a.contigId == b.contigId && a.contigEnd < b.contigEnd);
    }
};

template <typename TMatch>
struct MatchSorterByErrors
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return a.errors <= b.errors;
    }
};

// ----------------------------------------------------------------------------
// Class MatchesCounter
// ----------------------------------------------------------------------------

template <typename TReadSeqs, typename TSequencing = SingleEnd>
struct MatchesCounter
{
    TReadSeqs const &   readSeqs;
    String<bool>        matched;

    MatchesCounter(TReadSeqs const & readSeqs) :
        readSeqs(readSeqs)
    {
        resize(matched, getReadsCount(readSeqs), false, Exact());
    }

    template <typename TMatch>
    void operator() (TMatch const & match)
    {
        matched[getReadId(readSeqs, match.readId)] = true;
    }
};

template <typename TReadSeqs>
struct MatchesCounter<TReadSeqs, PairedEnd>
{
    TReadSeqs const &   readSeqs;
    String<bool>        matched;

    MatchesCounter(TReadSeqs const & readSeqs) :
        readSeqs(readSeqs)
    {
        resize(matched, getPairsCount(readSeqs), false, Exact());
    }

    template <typename TMatch>
    void operator() (TMatch const & match)
    {
        matched[getPairId(readSeqs, match.readId)] = true;
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function isDuplicateBegin()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicateBegin(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return a.readId == b.readId && a.contigId == b.contigId && a.contigBegin == b.contigBegin;
}

// ----------------------------------------------------------------------------
// Function isDuplicateEnd()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicateEnd(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return a.readId == b.readId && a.contigId == b.contigId && a.contigEnd == b.contigEnd;
}

// ----------------------------------------------------------------------------
// Function removeDuplicateMatches()
// ----------------------------------------------------------------------------

template <typename TMatches, typename TThreading>
inline void removeDuplicateMatches(TMatches & matches, TThreading const & threading)
{
    typedef typename Iterator<TMatches, Standard>::Type         TMatchesIterator;
    typedef typename Value<TMatches>::Type                      TMatch;

    TMatchesIterator matchesBegin;
    TMatchesIterator matchesEnd;
    TMatchesIterator newIt;
    TMatchesIterator oldIt;

    matchesBegin = begin(matches, Standard());
    matchesEnd   = end(matches, Standard());
    newIt = matchesBegin;
    oldIt = matchesBegin;

//    stableSort(matches, MatchSorterByReadId<TMatch>(), threading);

    stableSort(matches, MatchSorterByEndPos<TMatch>(), threading);

    // Remove duplicates by end position.
    while (oldIt != matchesEnd)
    {
        *newIt = *oldIt;

        ++oldIt;

        while (oldIt != matchesEnd && isDuplicateEnd(*newIt, *oldIt)) ++oldIt;

        ++newIt;
    }

    resize(matches, newIt - matchesBegin, Exact());

    matchesBegin = begin(matches, Standard());
    matchesEnd   = end(matches, Standard());
    newIt = matchesBegin;
    oldIt = matchesBegin;

    stableSort(matches, MatchSorterByBeginPos<TMatch>(), threading);

    // Remove duplicates by begin position.
    while (oldIt != matchesEnd)
    {
        *newIt = *oldIt;

        ++oldIt;

        while (oldIt != matchesEnd && isDuplicateBegin(*newIt, *oldIt)) ++oldIt;

        ++newIt;
    }

    resize(matches, newIt - matchesBegin, Exact());
}

// ----------------------------------------------------------------------------
// Function sortByErrors()
// ----------------------------------------------------------------------------

template <typename TMatches, typename TThreading>
inline void sortByErrors(TMatches & matches, TThreading const & threading)
{
    typedef typename Value<TMatches>::Type  TMatch;

    sort(matches, MatchSorterByErrors<TMatch>(), threading);
}

// ----------------------------------------------------------------------------
// Function getCount()
// ----------------------------------------------------------------------------

template <typename TReadSeqs, typename TSequencing, typename TThreading>
inline typename Size<TReadSeqs>::Type
getCount(MatchesCounter<TReadSeqs, TSequencing> const & counter, TThreading const & threading)
{
    return count(counter.matched, true, threading);
}

// ----------------------------------------------------------------------------
// Function countMatches()
// ----------------------------------------------------------------------------

template <typename TReadSeqs, typename TMatches, typename TSequencing, typename TThreading>
inline typename Size<TReadSeqs>::Type
countMatches(TReadSeqs const & readSeqs, TMatches const & matches, TSequencing const & /* tag */, TThreading const & threading)
{
    return getCount(forEach(matches, MatchesCounter<TReadSeqs, TSequencing>(readSeqs), threading), threading);
}

#endif  // #ifndef APP_CUDAMAPPER_BITS_MATCHES_H_
