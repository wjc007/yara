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
// Forwards
// ============================================================================

//template <typename THaystack, typename TNeedle, typename TSpec>
//struct Extender;
//
//template <typename THaystack, typename TNeedle, typename TSpec>
//struct Verifier;

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

template <typename TReadSeqs, typename TSpec = void>
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

// ----------------------------------------------------------------------------
// Class MatchesManager
// ----------------------------------------------------------------------------

template <typename TMatches, typename TConfig = void>
struct MatchesManager
{
    typedef typename Value<TMatches>::Type  TMatch;

    TMatches & matches;
    TMatch prototype;

    MatchesManager(TMatches & matches) :
        matches(matches),
        prototype()
    {}

    template <typename THaystackPos, typename TErrors>
    void operator() (THaystackPos matchBegin, THaystackPos matchEnd, TErrors errors)
    {
        SEQAN_ASSERT_EQ(getValueI1(matchBegin), getValueI1(matchEnd));

        prototype.contigId = getValueI1(matchBegin);
        prototype.contigBegin = getValueI2(matchBegin);
        prototype.contigEnd = getValueI2(matchEnd);
        prototype.errors = errors;
        appendValue(matches, prototype);
    }

//    template <typename THaystack, typename TNeedle, typename TSpec>
//    void operator() (Extender<THaystack, TNeedle, TSpec> const & extender)
//    {
//        SEQAN_ASSERT_EQ(getValueI1(extender.matchBegin), getValueI1(extender.matchEnd));
//
//        prototype.contigBegin = getValueI2(extender.matchBegin);
//        prototype.contigEnd = getValueI2(extender.matchEnd);
//        prototype.errors = extender.errors;
//        appendValue(matches, prototype);
//    }

//    template <typename THaystack, typename TNeedle, typename TSpec>
//    void operator() (Verifier<THaystack, TNeedle, TSpec> const & verifier)
//    {
//        appendValue(matches, prototype);
//    }
};

// ----------------------------------------------------------------------------
// Class MatchesManager
// ----------------------------------------------------------------------------

template <typename TMatches>
struct MatchesManager<TMatches, AnyBest>
{
    typedef typename Value<TMatches>::Type  TMatch;

    String<unsigned char> minErrors;
    TMatch prototype;

    MatchesManager(TMatches & /* matches */, TReadSeqs & readSeqs) :
        prototype()
    {
        resize(minErrors, getReadsCount(readSeqs), MaxValue<unsigned char>::VALUE, Exact());
    }

    template <typename THaystackPos, typename TErrors>
    void operator() (THaystackPos /* matchBegin */, THaystackPos /* matchEnd */, TErrors errors)
    {
        minErrors[prototype.readId] = _min(minErrors[prototype.readId], errors);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function fill()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TReadId, typename TContigId, typename TContigPos, typename TErrors>
inline void fill(Match<TSpec> & match,
                 TReadId readId,
                 TContigId contigId,
                 TContigPos contigBegin,
                 TContigPos contigEnd,
                 TErrors errors)
{
    match.readId = readId;
    match.contigBegin = contigBegin;
    match.contigEnd = contigEnd;
    match.contigId = contigId;
    match.errors = errors;
}

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

//    std::stable_sort(matchesBegin, matchesEnd, MatchSorterByReadId<TMatch>());

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

// ----------------------------------------------------------------------------
// Function getCount()
// ----------------------------------------------------------------------------

template <typename TReadSeqs, typename TSpec>
inline typename Size<TReadSeqs>::Type
getCount(MatchesCounter<TReadSeqs, TSpec> const & counter)
{
    return std::count(begin(counter.matched, Standard()), end(counter.matched, Standard()), true);
}

// ----------------------------------------------------------------------------
// Function countMatches()
// ----------------------------------------------------------------------------

template <typename TReadSeqs, typename TMatches, typename TSpec>
inline typename Size<TReadSeqs>::Type
countMatches(TReadSeqs const & readSeqs, TMatches const & matches, TSpec const & /* tag */)
{
    return getCount(std::for_each(begin(matches, Standard()),
                                  end(matches, Standard()),
                                  MatchesCounter<TReadSeqs, TSpec>(readSeqs)));
}

#endif  // #ifndef APP_CUDAMAPPER_MATCHES_H_
