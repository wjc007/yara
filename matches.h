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

#ifndef APP_CUDAMAPPER_MATCHES_H_
#define APP_CUDAMAPPER_MATCHES_H_

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
    unsigned        beginPos;
    short           endPosDelta;
    unsigned char   contigId;
    unsigned char   errors;
};

// ----------------------------------------------------------------------------
// Class MatchSorterByXXX
// ----------------------------------------------------------------------------

template <typename TMatch>
struct MatchSorterByReadId :
    std::binary_function<TMatch, TMatch, int>
{
    inline int operator()(TMatch const & a, TMatch const & b) const
    {
        return a.readId - b.readId;
    }
};

template <typename TMatch>
struct MatchSorterByBeginPos
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return (a.contigId < b.contigId) ||
               (a.contigId == b.contigId && (isForward(a) && isReverse(b))) ||
               (a.contigId == b.contigId && !(isReverse(a) && isForward(b)) && a.beginPos < b.beginPos);
    }
};

template <typename TMatch>
struct MatchSorterByEndPos
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return (a.contigId < b.contigId) ||
               (a.contigId == b.contigId && (isForward(a) && isReverse(b))) ||
               (a.contigId == b.contigId && !(isReverse(a) && isForward(b)) && endPos(a) < endPos(b));
    }
};

template <typename TMatch>
struct MatchSorterByErrors
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return a.errors < b.errors;
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function fill()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void fill(Match<TSpec> & match,
                 TContigId contigId,
                 Pair<TContigPos> contigPos,
                 TReadId readId,
                 TErrors errors,
                 bool isReverse)
{
    match.readId = readId;
    match.beginPos = getValueI1(contigPos);
    match.endPosDelta = getValueI2(contigPos) - getValueI1(contigPos);
    if (isReverse)
        match.endPosDelta = -match.endPosDelta;
    match.contigId = contigId;
    match.errors = errors;
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void assign(Match<TSpec> & dest, Match<TSpec> const & source)
{
    dest.readId = source.readId;
    dest.beginPos = source.beginPos;
    dest.endPosDelta = source.endPosDelta;
    dest.contigId = source.contigId;
    dest.errors = source.errors;
}

// ----------------------------------------------------------------------------
// Function isForward()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isForward(Match<TSpec> const & match)
{
    return match.endPosDelta > 0;
}

// ----------------------------------------------------------------------------
// Function isReverse()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isReverse(Match<TSpec> const & match)
{
    return match.endPosDelta < 0;
}

// ----------------------------------------------------------------------------
// Function isConcordant()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isConcordant(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return (isForward(a) && isForward(b)) || (isReverse(a) && isReverse(b));
}

// ----------------------------------------------------------------------------
// Function isDuplicateBegin()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicateBegin(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return a.readId == b.readId && a.contigId == b.contigId && isConcordant(a, b) && a.beginPos == b.beginPos;
}

// ----------------------------------------------------------------------------
// Function isDuplicateEnd()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicateEnd(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return a.readId == b.readId && a.contigId == b.contigId && isConcordant(a, b) && endPos(a) == endPos(b);
}

// ----------------------------------------------------------------------------
// Function endPos()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned endPos(Match<TSpec> const & match)
{
    return match.beginPos + abs(match.endPosDelta);
}

// ----------------------------------------------------------------------------
// Function getContigPos()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline Pair<unsigned> getContigPos(Match<TSpec> const & match)
{
    return Pair<unsigned>(match.beginPos, endPos(match));
}

// ----------------------------------------------------------------------------
// Function removeDuplicateMatches()
// ----------------------------------------------------------------------------

template <typename TMatches>
inline void removeDuplicateMatches(TMatches & matches)
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

//    std::sort(matchesBegin, matchesEnd, MatchSorterByReadId<TMatch>());

    std::stable_sort(matchesBegin, matchesEnd, MatchSorterByEndPos<TMatch>());

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

    std::stable_sort(matchesBegin, matchesEnd, MatchSorterByBeginPos<TMatch>());

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

template <typename TMatches>
inline void sortByErrors(TMatches & matches)
{
    typedef typename Iterator<TMatches, Standard>::Type         TMatchesIterator;
    typedef typename Value<TMatches>::Type                      TMatch;

    TMatchesIterator matchesBegin = begin(matches, Standard());
    TMatchesIterator matchesEnd = end(matches, Standard());

    std::sort(matchesBegin, matchesEnd, MatchSorterByErrors<TMatch>());
}

#endif  // #ifndef APP_CUDAMAPPER_MATCHES_H_
