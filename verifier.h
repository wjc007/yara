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

#ifndef APP_CUDAMAPPER_VERIFIER_H_
#define APP_CUDAMAPPER_VERIFIER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class VerifierConfig
// ----------------------------------------------------------------------------

template <typename TOptions_, typename TIndex_, typename TContigs_, typename TReadSeqs_>
struct VerifierConfig
{
    typedef TOptions_       TOptions;
    typedef TIndex_         TIndex;
    typedef TContigs_       TContigs;
    typedef TReadSeqs_      TReadSeqs;
};

// ----------------------------------------------------------------------------
// Class Verifier
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig>
struct Verifier
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

    TOptions const &    options;
    TIndex &            index;
    TContigs &          contigs;

    TPatternState       patternState;
    TPatternStateRev    patternStateRev;

    unsigned readErrors;
    unsigned seedErrors;
    unsigned seedLength;
    unsigned hitsThreshold;
    unsigned long verificationsCount;
    unsigned long matchesCount;

    Verifier(TOptions const & options, TIndex & index, TContigs & contigs) :
        options(options),
        index(index),
        contigs(contigs),
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
inline bool _extendLeft(Verifier<TExecSpace, TConfig> & verifier,
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

    TErrors remainingErrors = verifier.readErrors - errors;
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

    return errors <= verifier.readErrors;
}

// ----------------------------------------------------------------------------
// Function _extendRight()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TPatternState,
          typename TContigInfix, typename TReadInfix, typename TContigPos, typename TErrors>
inline bool _extendRight(Verifier<TExecSpace, TConfig> & verifier,
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
        return errors <= verifier.readErrors;
    }
    setBeginPosition(contigInfix, beginPosition(contigInfix) + lcp);
    setBeginPosition(readInfix, beginPosition(readInfix) + lcp);

    // NOTE Uncomment this to disable lcp trick.
//    TContigPos lcp = 0;

    TErrors remainingErrors = verifier.readErrors - errors;
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

    return errors <= verifier.readErrors;
}

// ----------------------------------------------------------------------------
// Function verifyHit()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs,
          typename TReadId, typename TReadPos, typename TContigId, typename TContigPos, typename TErrors>
inline bool verifyHit(Verifier<TExecSpace, TConfig> & verifier,
                      TReadSeqs & readSeqs,
                      TReadId readId,
                      TReadPos readBegin,
                      TReadPos readEnd,
                      TContigId contigId,
                      TContigPos contigBegin,
                      TContigPos contigEnd,
                      TErrors hitErrors)
{
    typedef Verifier<TExecSpace, TConfig>                               TVerifier;
    typedef typename TVerifier::TContig                                 TContig;
    typedef typename Size<TContig>::Type                                TContigSize;
    typedef typename Infix<TContig>::Type                               TContigInfix;
    typedef typename TVerifier::TReadSeq                                TReadSeq;
    typedef typename Infix<TReadSeq>::Type                              TReadInfix;

    TContig contig = verifier.contigs[contigId];
    TReadSeq readSeq = readSeqs[readId];
    TContigSize contigLength = length(contig);
    TErrors readLength = length(readSeq);
    TErrors readErrors = hitErrors;

    // Extend left.
    TContigPos matchBegin = contigBegin;

    if (readBegin > 0)
    {
        TContigPos contigLeftBegin = 0;
        if (contigBegin > readBegin + verifier.readErrors - readErrors)
            contigLeftBegin = contigBegin - (readBegin + verifier.readErrors - readErrors);

        TContigInfix contigLeft = infix(contig, contigLeftBegin, contigBegin);
        TReadInfix readLeft = infix(readSeq, 0, readBegin);

        if (!_extendLeft(verifier, verifier.patternStateRev, contigLeft, readLeft, matchBegin, readErrors))
            return false;
    }

    // Extend right.
    TContigPos matchEnd = contigEnd;

    if (readEnd < readLength)
    {
        TContigPos contigRightEnd = contigLength;
        if (contigRightEnd > contigBegin + readLength - readBegin + verifier.readErrors - readErrors)
            contigRightEnd = contigBegin + readLength - readBegin + verifier.readErrors - readErrors;

        if (contigBegin + verifier.seedLength >= contigRightEnd)
            return false;

        TContigInfix contigRight = infix(contig, contigEnd, contigRightEnd);
        TReadInfix readRight = infix(readSeq, readEnd, readLength);

        if (!_extendRight(verifier, verifier.patternState, contigRight, readRight, matchEnd, readErrors))
            return false;
    }

    return true;
}

// ----------------------------------------------------------------------------
// Function anchorRead()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename THits, typename TSA, typename TReadId>
inline void anchorRead(Verifier<TExecSpace, TConfig> & verifier, TReadSeqs & readSeqs, THits const & hits, TSA const & sa, TReadId anchorId)
{
    typedef Verifier<TExecSpace, TConfig>                               TVerifier;
    typedef typename TVerifier::TContig                                 TContig;
    typedef typename Size<TContig>::Type                                THitPos;
    typedef typename Size<THits>::Type                                  THitId;
    typedef typename Value<TSA>::Type                                   THit;

    // Consider the hits of all seeds of the anchor.
    THitId hitsBegin = anchorId * (verifier.readErrors + 1);
    THitId hitsEnd = (anchorId + 1) * (verifier.readErrors + 1);
    for (THitId hitId = hitsBegin; hitId < hitsEnd; ++hitId)
    {
        // Verify all hits of a seed of the anchor.
        for (THitPos hitPos = getValueI1(hits.ranges[hitId]); hitPos < getValueI2(hits.ranges[hitId]); ++hitPos)
        {
            THit hit = toSuffixPosition(verifier.index, sa[hitPos], verifier.seedLength);

            if (verifyHit(verifier,
                          readSeqs,
                          anchorId, (hitId - hitsBegin) * verifier.seedLength, (hitId - hitsBegin + 1) * verifier.seedLength,
                          getValueI1(hit), getValueI2(hit), getValueI2(hit) + verifier.seedLength,
                          verifier.seedErrors))
            {
                verifier.matchesCount++;
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function anchorPairs()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename THits, typename TSA>
inline void anchorPairs(Verifier<TExecSpace, TConfig> & verifier, TReadSeqs & readSeqs, THits const & hits, TSA const & sa)
{
    typedef typename Size<TReadSeqs>::Type                              TReadId;

#ifdef ENABLE_GENOME_LOADING
    setValue(verifier.index.text, verifier.contigs);
#endif

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

        TReadId anchorOneTwoId = (pairOneTwoHits == fwdOneHits) ? fwdOneId : revTwoId;
        TReadId anchorTwoOneId = (pairTwoOneHits == fwdTwoHits) ? fwdTwoId : revOneId;

        // Skip the pair if the anchor is hard.
        if (pairOneTwoHits + pairTwoOneHits > verifier.hitsThreshold) continue;

        verifier.verificationsCount += pairOneTwoHits + pairTwoOneHits;

#ifdef ENABLE_GENOME_LOADING
        anchorRead(verifier, readSeqs, hits, sa, anchorOneTwoId);
        anchorRead(verifier, readSeqs, hits, sa, anchorTwoOneId);
#endif
    }
}

#endif  // #ifndef APP_CUDAMAPPER_VERIFIER_H_
