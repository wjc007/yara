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

#ifndef APP_CUDAMAPPER_EXTENDER_H_
#define APP_CUDAMAPPER_EXTENDER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ExtenderConfig
// ----------------------------------------------------------------------------

template <typename TOptions_, typename TIndex_, typename TContigs_, typename TReadSeqs_>
struct ExtenderConfig
{
    typedef TOptions_       TOptions;
    typedef TIndex_         TIndex;
    typedef TContigs_       TContigs;
    typedef TReadSeqs_      TReadSeqs;
};

// ----------------------------------------------------------------------------
// Class Extender
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig>
struct Extender
{
    typedef typename TConfig::TOptions                                  TOptions;
    typedef typename TConfig::TIndex                                    TIndex;
    typedef typename TConfig::TContigs                                  TContigs;
    typedef typename TConfig::TReadSeqs                                 TReadSeqs;

    typedef typename Value<TContigs>::Type                              TContig;
    typedef typename Value<TReadSeqs>::Type                             TReadSeq;
    typedef typename Infix<TReadSeq>::Type                              TReadInfix;
    typedef ModifiedString<TReadInfix, ModReverse>                      TReadInfixRev;

    typedef AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_>   TMyersSpec;
    typedef Myers<TMyersSpec, True, void>                               TAlgorithm;
    typedef PatternState_<TReadInfix, TAlgorithm>                       TPatternState;
    typedef PatternState_<TReadInfixRev, TAlgorithm>                    TPatternStateRev;

    typedef Verifier<TExecSpace, TConfig>                               TVerifier;

    TOptions const &    options;
    TIndex &            index;
    TContigs &          contigs;

    TPatternState       patternState;
    TPatternStateRev    patternStateRev;

    TVerifier verifier;

    unsigned readErrors;
    unsigned seedErrors;
    unsigned seedLength;
    unsigned hitsThreshold;
    unsigned long verificationsCount;
    unsigned long matchesCount;

    Extender(TOptions const & options, TIndex & index, TContigs & contigs) :
        options(options),
        index(index),
        contigs(contigs),
        verifier(options, contigs),
        readErrors(5),
        seedErrors(0),
        seedLength(16),
        hitsThreshold(300),
        verificationsCount(0),
        matchesCount(0)
    {}
};

// ----------------------------------------------------------------------------
// Function _extendLeft()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TPatternState,
          typename TContigInfix, typename TReadInfix, typename TContigPos, typename TErrors>
inline bool _extendLeft(Extender<TExecSpace, TConfig> & extender,
                        TPatternState & patternState,
                        TContigInfix & contigInfix,
                        TReadInfix & readInfix,
                        TContigPos & matchBegin,
                        TErrors & errors)
{
    typedef ModifiedString<TReadInfix, ModReverse>          TReadInfixRev;
    typedef ModifiedString<TContigInfix, ModReverse>        TContigInfixRev;
    typedef Finder<TContigInfixRev>                         TFinder;

    // Lcp trick.
    TContigPos lcp = 0;
    {  // TODO(holtgrew): Workaround to storing and returning copies in host() for nested infixes/modified strings. This is ugly and should be fixed later.
        TReadInfixRev readInfixRev(readInfix);
        TContigInfixRev contigInfixRev(contigInfix);
        lcp = lcpLength(contigInfixRev, readInfixRev);
    }
    if (lcp == length(readInfix))
    {
        matchBegin -= lcp;
        return true;
    }
    setEndPosition(contigInfix, endPosition(contigInfix) - lcp);
    setEndPosition(readInfix, endPosition(readInfix) - lcp);

    TErrors remainingErrors = extender.readErrors - errors;
    TErrors minErrors = remainingErrors + 1;
    TContigPos endPos = 0;

    // Stop seed extension.
    if (!remainingErrors)
        return false;

    // Align.
    TReadInfixRev readInfixRev(readInfix);
    TContigInfixRev contigInfixRev(contigInfix);
    TFinder finder(contigInfixRev);
    patternState.leftClip = remainingErrors;

    // TODO(esiragusa): Use a generic type for errors.
    while (find(finder, readInfixRev, patternState, -static_cast<int>(remainingErrors)))
    {
        TErrors currentErrors = -getScore(patternState);

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = position(finder) + 1;
        }
    }

    errors += minErrors;
    matchBegin -= endPos + lcp;

    return errors <= extender.readErrors;
}

// ----------------------------------------------------------------------------
// Function _extendRight()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TPatternState,
          typename TContigInfix, typename TReadInfix, typename TContigPos, typename TErrors>
inline bool _extendRight(Extender<TExecSpace, TConfig> & extender,
                         TPatternState & patternState,
                         TContigInfix & contigInfix,
                         TReadInfix & readInfix,
                         TContigPos & matchEnd,
                         TErrors & errors)
{
    typedef Finder<TContigInfix>    TFinder;

    // Lcp trick.
    TContigPos lcp = lcpLength(contigInfix, readInfix);
    if (lcp == length(readInfix))
    {
        matchEnd += lcp;
        return true;
    }
    else if (lcp == length(contigInfix))
    {
        errors += length(readInfix) - length(contigInfix);
        matchEnd += length(readInfix);
        return errors <= extender.readErrors;
    }
    setBeginPosition(contigInfix, beginPosition(contigInfix) + lcp);
    setBeginPosition(readInfix, beginPosition(readInfix) + lcp);

    // NOTE Uncomment this to disable lcp trick.
//    TContigPos lcp = 0;

    TErrors remainingErrors = extender.readErrors - errors;
    TErrors minErrors = remainingErrors + 1;
    TContigPos endPos = 0;

    // NOTE Comment this to disable lcp trick.
    // Stop seed extension.
    if (!remainingErrors)
        return false;

    // Remove last base.
    TContigInfix contigPrefix(contigInfix);
    TReadInfix readPrefix(readInfix);
    setEndPosition(contigPrefix, endPosition(contigPrefix) - 1);
    setEndPosition(readPrefix, endPosition(readPrefix) - 1);

    // Align.
    TFinder finder(contigPrefix);
    patternState.leftClip = remainingErrors;

    while (find(finder, readPrefix, patternState, -static_cast<int>(remainingErrors)))
    {
        TContigPos currentEnd = position(finder) + 1;
        TErrors currentErrors = -getScore(patternState);

        // Compare last base.
        if (contigInfix[currentEnd] != back(readInfix))
            if (++currentErrors > remainingErrors)
                continue;

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = currentEnd;
        }
    }

    errors += minErrors;
    matchEnd += endPos + lcp + 1;

    return errors <= extender.readErrors;
}

// ----------------------------------------------------------------------------
// Function verifyHit()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs,
          typename TReadId, typename TReadPos, typename TContigId, typename TContigPos, typename TErrors>
inline bool verifyHit(Extender<TExecSpace, TConfig> & extender,
                      TReadSeqs & readSeqs,
                      TReadId readId,
                      TReadPos readBegin,
                      TReadPos readEnd,
                      TContigId contigId,
                      TContigPos contigBegin,
                      TContigPos contigEnd,
                      TErrors hitErrors)
{
    typedef Extender<TExecSpace, TConfig>                               TExtender;
    typedef typename TExtender::TContig                                 TContig;
    typedef typename Size<TContig>::Type                                TContigSize;
    typedef typename Infix<TContig>::Type                               TContigInfix;
    typedef typename TExtender::TReadSeq                                TReadSeq;
    typedef typename Infix<TReadSeq>::Type                              TReadInfix;

    TContig contig = extender.contigs[contigId];
    TReadSeq readSeq = readSeqs[readId];
    TContigSize contigLength = length(contig);
    TErrors readLength = length(readSeq);
    TErrors readErrors = hitErrors;

    // Extend left.
    TContigPos matchBegin = contigBegin;

    if (readBegin > 0)
    {
        TContigPos contigLeftBegin = 0;
        if (contigBegin > readBegin + extender.readErrors - readErrors)
            contigLeftBegin = contigBegin - (readBegin + extender.readErrors - readErrors);

        TContigInfix contigLeft = infix(contig, contigLeftBegin, contigBegin);
        TReadInfix readLeft = infix(readSeq, 0, readBegin);

        if (!_extendLeft(extender, extender.patternStateRev, contigLeft, readLeft, matchBegin, readErrors))
            return false;
    }

    // Extend right.
    TContigPos matchEnd = contigEnd;

    if (readEnd < readLength)
    {
        TContigPos contigRightEnd = contigLength;
        if (contigRightEnd > contigBegin + readLength - readBegin + extender.readErrors - readErrors)
            contigRightEnd = contigBegin + readLength - readBegin + extender.readErrors - readErrors;

        if (contigBegin + extender.seedLength >= contigRightEnd)
            return false;

        TContigInfix contigRight = infix(contig, contigEnd, contigRightEnd);
        TReadInfix readRight = infix(readSeq, readEnd, readLength);

        if (!_extendRight(extender, extender.patternState, contigRight, readRight, matchEnd, readErrors))
            return false;
    }

    return true;
}

// ----------------------------------------------------------------------------
// Function anchorRead()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename THits, typename TSA, typename TReadId>
inline void anchorRead(Extender<TExecSpace, TConfig> & extender, TReadSeqs & readSeqs, THits const & hits, TSA const & sa, TReadId anchorId, TReadId mateId)
{
    typedef Extender<TExecSpace, TConfig>                               TExtender;
    typedef typename TExtender::TContig                                 TContig;
    typedef typename TExtender::TReadSeq                                TReadSeq;
    typedef typename Value<TSA>::Type                                   THit;
    typedef typename Size<THits>::Type                                  THitId;
    typedef typename Size<TSA>::Type                                    THitPos;
    typedef typename Value<THit, 1>::Type                               TContigId;
    typedef typename Size<TContig>::Type                                TContigPos;
    typedef typename Size<TReadSeq>::Type                               TReadPos;

    // Consider the hits of all seeds of the anchor.
    THitId hitsBegin = anchorId * (extender.readErrors + 1);
    THitId hitsEnd = (anchorId + 1) * (extender.readErrors + 1);
    for (THitId hitId = hitsBegin; hitId < hitsEnd; ++hitId)
    {
        // Verify all hits of a seed of the anchor.
        for (THitPos hitPos = getValueI1(hits.ranges[hitId]); hitPos < getValueI2(hits.ranges[hitId]); ++hitPos)
        {
            THit hit = sa[hitPos];
            setSeqOffset(hit, suffixLength(hit, extender.contigs) - extender.seedLength);
//            THit hit = toSuffixPosition(extender.index, sa[hitPos], extender.seedLength);

            TReadPos readBegin = (hitId - hitsBegin) * extender.seedLength;
            TReadPos readEnd = (hitId - hitsBegin + 1) * extender.seedLength;
            TContigId contigId = getValueI1(hit);
            TContigPos contigBegin = getValueI2(hit);
            TContigPos contigEnd = getValueI2(hit) + extender.seedLength;

            if (verifyHit(extender, readSeqs,
                          anchorId, readBegin, readEnd,
                          contigId, contigBegin, contigEnd,
                          extender.seedErrors))
            {
                extender.matchesCount++;
                TContigPos matchBegin = contigBegin;
                TContigPos matchEnd = contigEnd;

                if (findMate(extender.verifier, readSeqs, contigId, matchBegin, matchEnd, mateId))
                    extender.verifier.matchesCount++;
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function anchorPairs()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename THits, typename TSA>
inline void anchorPairs(Extender<TExecSpace, TConfig> & extender, TReadSeqs & readSeqs, THits const & hits, TSA const & sa)
{
    typedef typename Size<TReadSeqs>::Type                              TReadId;

    TReadId pairsCount = length(readSeqs) / 4;

    for (TReadId pairId = 0; pairId < pairsCount; ++pairId)
    {
        // Get mates ids.
        TReadId fwdOneId = pairId;
        TReadId fwdTwoId = pairId + pairsCount;
        TReadId revOneId = pairId + 2 * pairsCount;
        TReadId revTwoId = pairId + 3 * pairsCount;

        // Choose the anchor.
        unsigned long fwdOneHits = countHits(hits, fwdOneId);
        unsigned long fwdTwoHits = countHits(hits, fwdTwoId);
        unsigned long revOneHits = countHits(hits, revOneId);
        unsigned long revTwoHits = countHits(hits, revTwoId);

        unsigned long pairOneTwoHits = std::min(fwdOneHits, revTwoHits);
        unsigned long pairTwoOneHits = std::min(fwdTwoHits, revOneHits);

        // Skip the pair if the anchor is hard.
        if (pairOneTwoHits + pairTwoOneHits > extender.hitsThreshold) continue;

        TReadId anchorOneTwoId = (pairOneTwoHits == fwdOneHits) ? fwdOneId : revTwoId;
        TReadId anchorTwoOneId = (pairTwoOneHits == fwdTwoHits) ? fwdTwoId : revOneId;
        TReadId mateOneTwoId = (pairOneTwoHits == fwdOneHits) ? revTwoId : fwdOneId;
        TReadId mateTwoOneId = (pairTwoOneHits == fwdTwoHits) ? revOneId : fwdTwoId;

        extender.verificationsCount += pairOneTwoHits + pairTwoOneHits;

        anchorRead(extender, readSeqs, hits, sa, anchorOneTwoId, mateOneTwoId);
        anchorRead(extender, readSeqs, hits, sa, anchorTwoOneId, mateTwoOneId);
    }
}

#endif  // #ifndef APP_CUDAMAPPER_EXTENDER_H_
