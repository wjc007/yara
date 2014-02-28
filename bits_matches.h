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
// Tags
// ============================================================================

struct SortErrors_;
typedef Tag<SortErrors_> const SortErrors;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Match
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Match
{
    unsigned        readId       : 22;
    unsigned char   contigId     : 8;
    unsigned        contigBegin  : 30;
    unsigned short  contigEnd    : 14;
    bool            isFwd        : 1;
    unsigned char   errors       : 5;
}
__attribute__((packed));

// ----------------------------------------------------------------------------
// Class MatchReadId
// ----------------------------------------------------------------------------
// NOTE(esiragusa): MatchReadId<TMatch> functor could be MemberGetter<TMatch, TMember>
// NOTE(esiragusa): getReadId<TMatch>() could be getMember(TMatch, TMember())

template <typename TMatch>
struct MatchReadId
{
    inline unsigned operator()(TMatch const & me) const
    {
        return getReadId(me);
    }
};

// ----------------------------------------------------------------------------
// Class MatchSorter
// ----------------------------------------------------------------------------

template <typename TMatch, typename TSpec>
struct MatchSorter;

template <typename TMatch>
struct MatchSorter<TMatch, SortReadId>
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return getReadId(a) < getReadId(b);
    }
};

template <typename TMatch>
struct MatchSorter<TMatch, SortBeginPos>
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return contigLess(a, b) || (contigEqual(a, b) && getContigBegin(a) < getContigBegin(b));
    }
};

template <typename TMatch>
struct MatchSorter<TMatch, SortEndPos>
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return contigLess(a, b) || (contigEqual(a, b) && getContigEnd(a) < getContigEnd(b));
    }
};

template <typename TMatch>
struct MatchSorter<TMatch, SortErrors>
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return getErrors(a) < getErrors(b);
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
        matched[getReadId(match)] = true;
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
        matched[getPairId(readSeqs, getReadId(match))] = true;
    }
};

// ----------------------------------------------------------------------------
// Class DuplicateRemover
// ----------------------------------------------------------------------------

template <typename TCounts, typename TPosition>
struct DuplicateRemover
{
    TCounts &    unique;

    DuplicateRemover(TCounts & unique) :
        unique(unique)
    {}

    template <typename TIterator>
    void operator() (TIterator & it)
    {
        typedef typename Value<TIterator>::Type     TMatches;
        typedef typename Value<TMatches>::Type      TMatch;

        TMatches const & matches = value(it);

        sort(matches, MatchSorter<TMatch, TPosition>());
        unique[position(it) + 1] = compactUniqueMatches(matches, TPosition());
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Match Setters
// ----------------------------------------------------------------------------

template <typename TSpec, typename TReadSeqs, typename TReadSeqId>
inline void setReadId(Match<TSpec> & me, TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    me.readId = getReadId(readSeqs, readSeqId);
    me.isFwd = isFwdReadSeq(readSeqs, readSeqId);
}

template <typename TSpec, typename TContigPos>
inline void setContigPosition(Match<TSpec> & me, TContigPos contigBegin, TContigPos contigEnd)
{
    SEQAN_ASSERT_EQ(getValueI1(contigBegin), getValueI1(contigEnd));
    SEQAN_ASSERT_LT(getValueI2(contigBegin), getValueI2(contigEnd));

    me.contigId = getValueI1(contigBegin);
    me.contigBegin = getValueI2(contigBegin);
    me.contigEnd = getValueI2(contigEnd) - getValueI2(contigBegin);
}

// ----------------------------------------------------------------------------
// Match Getters
// ----------------------------------------------------------------------------

template <typename TReadSeqs, typename TSpec>
inline typename Size<TReadSeqs>::Type
getReadSeqId(Match<TSpec> const & me, TReadSeqs const & readSeqs)
{
    return onForwardStrand(me) ? getFirstMateFwdSeqId(readSeqs, me.readId) : getFirstMateRevSeqId(readSeqs, me.readId);
}

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
    return me.contigBegin + me.contigEnd;
}

template <typename TSpec>
inline bool onForwardStrand(Match<TSpec> const & me)
{
    return me.isFwd;
}

template <typename TSpec>
inline bool onReverseStrand(Match<TSpec> const & me)
{
    return !onForwardStrand(me);
}

template <typename TSpec>
inline unsigned char getScore(Match<TSpec> const & /* me */)
{
    return 254;
}

template <typename TSpec>
inline unsigned char getErrors(Match<TSpec> const & me)
{
    return me.errors;
}

// ----------------------------------------------------------------------------
// Function strandEqual()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool strandEqual(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return !(onForwardStrand(a) ^ onForwardStrand(b));
}

// ----------------------------------------------------------------------------
// Function strandLess()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool strandLess(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return onForwardStrand(a) && onReverseStrand(b);
}

// ----------------------------------------------------------------------------
// Function contigEqual()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool contigEqual(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return getContigId(a) == getContigId(b) && strandEqual(a, b);
}

// ----------------------------------------------------------------------------
// Function contigLess()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool contigLess(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return getContigId(a) < getContigId(b) || (getContigId(a) == getContigId(b) && strandLess(a, b));
}

// ----------------------------------------------------------------------------
// Function isDuplicate(SortBeginPos)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicate(Match<TSpec> const & a, Match<TSpec> const & b, SortBeginPos)
{
    return contigEqual(a, b) && getContigBegin(a) == getContigBegin(b);
}

// ----------------------------------------------------------------------------
// Function isDuplicate(SortEndPos)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicate(Match<TSpec> const & a, Match<TSpec> const & b, SortEndPos)
{
    return contigEqual(a, b) && getContigEnd(a) == getContigEnd(b);
}

// ----------------------------------------------------------------------------
// Function compactUniqueMatches()
// ----------------------------------------------------------------------------
// Compact unique matches at the beginning of their container.

template <typename TMatches, typename TPosition>
inline typename Size<TMatches>::Type
compactUniqueMatches(TMatches & matches, Tag<TPosition> const & posTag)
{
    typedef typename Iterator<TMatches, Standard>::Type         TMatchesIterator;

    TMatchesIterator matchesBegin = begin(matches, Standard());
    TMatchesIterator matchesEnd = end(matches, Standard());
    TMatchesIterator newIt = matchesBegin;
    TMatchesIterator oldIt = matchesBegin;

    while (oldIt != matchesEnd)
    {
        *newIt = *oldIt;

        ++oldIt;

        while (oldIt != matchesEnd && isDuplicate(*newIt, *oldIt, posTag)) ++oldIt;

        ++newIt;
    }

    return newIt - matchesBegin;
}

// ----------------------------------------------------------------------------
// Function removeDuplicates()
// ----------------------------------------------------------------------------
// Remove duplicate matches from a set of matches.

template <typename TMatchesSet, typename TThreading>
inline void removeDuplicates(TMatchesSet & matchesSet, TThreading const & threading)
{
    typedef typename StringSetLimits<TMatchesSet>::Type         TLimits;

    TLimits newLimits;
    resize(newLimits, length(stringSetLimits(matchesSet)));
    front(newLimits) = 0;

    // Sort matches by end position and move unique matches at the beginning.
    iterate(matchesSet, DuplicateRemover<TLimits, SortEndPos>(newLimits), Rooted(), threading);

    // Exclude duplicate matches at the end.
    assign(stringSetLimits(matchesSet), newLimits);
    _refreshStringSetLimits(matchesSet, threading);

    // Sort matches by begin position and move unique matches at the beginning.
    iterate(matchesSet, DuplicateRemover<TLimits, SortBeginPos>(newLimits), Rooted(), threading);

    // Exclude duplicate matches at the end.
    assign(stringSetLimits(matchesSet), newLimits);
    _refreshStringSetLimits(matchesSet, threading);
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

// ----------------------------------------------------------------------------
// Function countBestMatches()
// ----------------------------------------------------------------------------
// Count the number of cooptimal matches - ordering by errors is required.

template <typename TMatches>
inline typename Size<TMatches>::Type
countBestMatches(TMatches const & matches)
{
    typedef typename Iterator<TMatches const, Standard>::Type   TIter;
    typedef typename Size<TMatches>::Type                       TCount;

    TIter itBegin = begin(matches, Standard());
    TIter itEnd = end(matches, Standard());

    TCount count = 0;

    for (TIter it = itBegin; it != itEnd && getErrors(*it) <= getErrors(*itBegin); it++, count++) ;

    return count;
}

// ----------------------------------------------------------------------------
// Function sortMatches()
// ----------------------------------------------------------------------------

//template <typename TMatches, typename TKey>
//inline void sortMatches(TMatches & matches)
//{
//    typedef typename Value<TMatches>::Type  TMatch;
//
//    sort(matches, MatchSorter<TMatch, TKey>());
//}

template <typename TIterator, typename TKey>
inline void sortMatches(TIterator & it)
{
    typedef typename Value<TIterator>::Type TMatches;
    typedef typename Value<TMatches>::Type  TMatch;

    TMatches matches = value(it);
    sort(matches, MatchSorter<TMatch, TKey>());
}

// ----------------------------------------------------------------------------
// Function findSameContig()
// ----------------------------------------------------------------------------

template <typename TMatchesIterator, typename TMatches>
inline void findSameContig(TMatchesIterator & leftIt, TMatchesIterator & rightIt,
                           TMatches const & left, TMatches const & right)
{
    while (!atEnd(leftIt, left) && !atEnd(rightIt, right))
    {
        if (getContigId(*leftIt) < getContigId(*rightIt))
            findNextContig(leftIt, left, getContigId(*leftIt));
        else if (getContigId(*leftIt) > getContigId(*rightIt))
            findNextContig(rightIt, right, getContigId(*rightIt));
        else
            break;
    }
}

// ----------------------------------------------------------------------------
// Function findNextContig()
// ----------------------------------------------------------------------------

template <typename TMatchesIterator, typename TMatches, typename TContigId>
inline void findNextContig(TMatchesIterator & it, TMatches const & matches, TContigId contigId)
{
    while (!atEnd(it, matches) && getContigId(*it) <= contigId) ++it;
}

// ----------------------------------------------------------------------------
// Function findReverseStrand()
// ----------------------------------------------------------------------------

template <typename TMatchesIterator, typename TMatches, typename TContigId>
inline void goDifferentStrand(TMatchesIterator & it, TMatches const & matches, TContigId contigId)
{
    while (!atEnd(it, matches) && (getContigId(*it) <= contigId) && onForwardStrand(*it)) ++it;
}

// ----------------------------------------------------------------------------
// Function pairMatches()
// ----------------------------------------------------------------------------

template <typename TMatches, typename TDelegate>
inline void pairMatches(TMatches const & left, TMatches const & right, TDelegate & delegate)
{
    typedef typename Iterator<TMatches const, Standard>::Type   TIterator;
    typedef typename Infix<TMatches const>::Type                TInfix;

    TIterator leftIt = begin(left, Standard());
    TIterator rightIt = begin(right, Standard());

    do
    {
        TIterator leftBegin;
        TIterator rightBegin;

        // Find matches on the same contig.
        findSameContig(leftIt, rightIt, left, right);

        // Find matches on forward strand.
        leftBegin = leftIt;
        rightBegin = rightIt;
        findReverseStrand(leftIt, left, getContigId(*leftBegin));
        findReverseStrand(rightIt, right, getContigId(*rightBegin));
        TInfix leftFwd = infix(left, leftBegin, leftIt);
        TInfix rightFwd = infix(right, rightBegin, rightIt);

        // Find matches on reverse strand.
        leftBegin = leftIt;
        rightBegin = rightIt;
        findNextContig(leftIt, left, getContigId(*leftBegin));
        findNextContig(rightIt, right, getContigId(*rightBegin));
        TInfix leftRev = infix(left, leftBegin, leftIt);
        TInfix rightRev = infix(right, rightBegin, rightIt);

        enumeratePairs(leftFwd, rightRev, LeftFile(), RightFile(), delegate);
        enumeratePairs(rightFwd, leftRev, RightFile(), LeftFile(), delegate);
    }
    while (!atEnd(leftIt, left) && !atEnd(rightIt, right));
}

// ----------------------------------------------------------------------------
// Function enumeratePairs()
// ----------------------------------------------------------------------------

template <typename TMatches, typename TMateFwd, typename TMateRev, typename TDelegate>
inline void enumeratePairs(TMatches const & matchesFwd, TMatches const & matchesRev,
                           TMateFwd const & /* tag */, TMateRev const & /* tag */,
                           TDelegate & delegate)
{
    typedef typename Iterator<TMatches const, Standard>::Type TIterator;

    unsigned libraryLength = 200;
    unsigned libraryError = 100;

    TIterator matchesFwdBegin = begin(matchesFwd, Standard());
    TIterator matchesFwdEnd = end(matchesFwd, Standard());
    TIterator matchesRevBegin = begin(matchesRev, Standard());
    TIterator matchesRevEnd = end(matchesRev, Standard());

    // NOTE matchesFwd queue C= matchesRev queue, i.e. matchesFwdTail >= matchesRevTail && matchesFwdHead <= matchesRevHead

    TIterator matchesFwdIt = matchesFwdBegin;

    if (matchesFwdIt == matchesFwdEnd) return;

    unsigned matchesFwdLength = getContigEnd(*matchesFwdBegin) - getContigBegin(*matchesFwdBegin);

    // Get next matchesRev match.
    for (TIterator matchesRevIt = matchesRevBegin; matchesRevIt != matchesRevEnd; ++matchesRevIt)
    {
        // Compute matchesRev match queue.
        unsigned matchesRevHead = getContigEnd(*matchesRevIt) - libraryLength + libraryError + matchesFwdLength;
        unsigned matchesRevTail = getContigEnd(*matchesRevIt) - libraryLength - libraryError;

        // Seek next matchesFwd match after the tail of matchesRev match queue.
        while (matchesFwdIt != matchesFwdEnd && getContigBegin(*matchesFwdIt) < matchesRevTail)
            ++matchesFwdIt;

        TIterator matchesFwdTailIt = matchesFwdIt;
        if (matchesFwdTailIt == matchesFwdEnd)
            break;

        // Continue if matchesFwd tail is beyond matchesRev head.
        unsigned matchesFwdTail = getContigBegin(*matchesFwdTailIt);
        if (matchesFwdTail >= matchesRevHead)
            continue;

        // Seek next matchesFwd match after the head of matchesRev match queue.
        while (matchesFwdIt != matchesFwdEnd && getContigEnd(*matchesFwdIt) <= matchesRevHead)
            ++matchesFwdIt;
        TIterator matchesFwdHeadIt = matchesFwdIt;

        // Couple all fwds mates vs rev mate.
        for (TIterator matchesFwdQueueIt = matchesFwdTailIt; matchesFwdQueueIt != matchesFwdHeadIt; ++matchesFwdQueueIt)
            delegate(*matchesFwdQueueIt, *matchesRevIt);

        // Empty matchesFwd queue.
//        if (matchesFwdTailIt == matchesFwdHeadIt) continue;
    }
}

#endif  // #ifndef APP_CUDAMAPPER_BITS_MATCHES_H_
