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
// Class Mater
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig>
struct Mater
{
    typedef typename TConfig::TOptions                                  TOptions;
    typedef typename TConfig::TContigs                                  TContigs;
    typedef typename TConfig::TReadSeqs                                 TReadSeqs;

    typedef typename Value<TContigs>::Type                              TContig;
    typedef typename Value<TReadSeqs>::Type                             TReadSeq;

    typedef Myers<>                                                     TAlgorithm;
    typedef Pattern<TReadSeq, TAlgorithm>                               TPattern;

    TOptions const &    options;
    TContigs &          contigs;

    TPattern            pattern;

    unsigned readErrors;
    unsigned long matchesCount;

    Mater(TOptions const & options, TContigs & contigs) :
        options(options),
        contigs(contigs),
        readErrors(5),
        matchesCount(0)
    {
        _patternMatchNOfPattern(pattern, false);
        _patternMatchNOfFinder(pattern, false);
    }
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

    typedef Mater<TExecSpace, TConfig>                                  TMater;

    TOptions const &    options;
    TIndex &            index;
    TContigs &          contigs;

    TPatternState       patternState;
    TPatternStateRev    patternStateRev;

    TMater mater;

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
        mater(options, contigs),
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
inline void anchorRead(Verifier<TExecSpace, TConfig> & verifier, TReadSeqs & readSeqs, THits const & hits, TSA const & sa, TReadId anchorId, TReadId mateId)
{
    typedef Verifier<TExecSpace, TConfig>                               TVerifier;
    typedef typename TVerifier::TContig                                 TContig;
    typedef typename TVerifier::TReadSeq                                TReadSeq;
    typedef typename Value<TSA>::Type                                   THit;
    typedef typename Size<THits>::Type                                  THitId;
    typedef typename Size<TSA>::Type                                    THitPos;
    typedef typename Value<THit, 1>::Type                               TContigId;
    typedef typename Size<TContig>::Type                                TContigPos;
    typedef typename Size<TReadSeq>::Type                               TReadPos;

    // Consider the hits of all seeds of the anchor.
    THitId hitsBegin = anchorId * (verifier.readErrors + 1);
    THitId hitsEnd = (anchorId + 1) * (verifier.readErrors + 1);
    for (THitId hitId = hitsBegin; hitId < hitsEnd; ++hitId)
    {
        // Verify all hits of a seed of the anchor.
        for (THitPos hitPos = getValueI1(hits.ranges[hitId]); hitPos < getValueI2(hits.ranges[hitId]); ++hitPos)
        {
            THit hit = sa[hitPos];
            setSeqOffset(hit, suffixLength(hit, verifier.contigs) - verifier.seedLength);
//            THit hit = toSuffixPosition(verifier.index, sa[hitPos], verifier.seedLength);

            TReadPos readBegin = (hitId - hitsBegin) * verifier.seedLength;
            TReadPos readEnd = (hitId - hitsBegin + 1) * verifier.seedLength;
            TContigId contigId = getValueI1(hit);
            TContigPos contigBegin = getValueI2(hit);
            TContigPos contigEnd = getValueI2(hit) + verifier.seedLength;

            if (verifyHit(verifier, readSeqs,
                          anchorId, readBegin, readEnd,
                          contigId, contigBegin, contigEnd,
                          verifier.seedErrors))
            {
                verifier.matchesCount++;
                TContigPos matchBegin = contigBegin;
                TContigPos matchEnd = contigEnd;

                if (findMate(verifier.mater, readSeqs, contigId, matchBegin, matchEnd, mateId))
                    verifier.mater.matchesCount++;
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
        if (pairOneTwoHits + pairTwoOneHits > verifier.hitsThreshold) continue;

        TReadId anchorOneTwoId = (pairOneTwoHits == fwdOneHits) ? fwdOneId : revTwoId;
        TReadId anchorTwoOneId = (pairTwoOneHits == fwdTwoHits) ? fwdTwoId : revOneId;
        TReadId mateOneTwoId = (pairOneTwoHits == fwdOneHits) ? revTwoId : fwdOneId;
        TReadId mateTwoOneId = (pairTwoOneHits == fwdTwoHits) ? revOneId : fwdTwoId;

        verifier.verificationsCount += pairOneTwoHits + pairTwoOneHits;

        anchorRead(verifier, readSeqs, hits, sa, anchorOneTwoId, mateOneTwoId);
        anchorRead(verifier, readSeqs, hits, sa, anchorTwoOneId, mateTwoOneId);
    }
}

// ----------------------------------------------------------------------------
// Function _getContigInfix()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TContigId, typename TContigPos>
inline void _getContigInfix(Mater<TExecSpace, TConfig> & verifier,
                            TContigId contigId,
                            TContigPos matchBegin,
                            TContigPos /* matchEnd */,
                            TContigPos & infixBegin,
                            TContigPos & infixEnd,
                            RightMate)
{
    typedef Mater<TExecSpace, TConfig>              TMater;
    typedef typename TMater::TContig                TContig;
    typedef typename Size<TContig>::Type            TContigSize;

    TContigSize contigLength = length(verifier.contigs[contigId]);

    infixBegin = contigLength;
    if (infixBegin < matchBegin + verifier.options.libraryLength - verifier.options.libraryError)
        infixBegin = matchBegin + verifier.options.libraryLength - verifier.options.libraryError;

    infixEnd = contigLength;
    if (infixEnd < matchBegin + verifier.options.libraryLength + verifier.options.libraryError)
        infixEnd = matchBegin + verifier.options.libraryLength + verifier.options.libraryError;

    SEQAN_ASSERT_LEQ(infixBegin, infixEnd);
    SEQAN_ASSERT_LEQ(infixEnd - infixBegin, 2 * verifier.options.libraryError);
}

template <typename TExecSpace, typename TConfig, typename TContigId, typename TContigPos>
inline void _getContigInfix(Mater<TExecSpace, TConfig> & verifier,
                            TContigId /* contigId */,
                            TContigPos /* matchBegin */,
                            TContigPos matchEnd,
                            TContigPos & infixBegin,
                            TContigPos & infixEnd,
                            LeftMate)
{
    infixBegin = 0;
    if (matchEnd > verifier.options.libraryLength + verifier.options.libraryError)
        infixBegin = matchEnd - verifier.options.libraryLength - verifier.options.libraryError;

    infixEnd = 0;
    if (matchEnd + verifier.options.libraryError > verifier.options.libraryLength)
        infixEnd = matchEnd - verifier.options.libraryLength + verifier.options.libraryError;

    SEQAN_ASSERT_LEQ(infixBegin, infixEnd);
    SEQAN_ASSERT_LEQ(infixEnd - infixBegin, 2 * verifier.options.libraryError);
}

// ----------------------------------------------------------------------------
// Function findMate()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename TContigId, typename TContigPos, typename TReadId>
inline bool findMate(Mater<TExecSpace, TConfig> & verifier,
                     TReadSeqs & readSeqs,
                     TContigId contigId,
                     TContigPos matchBegin,
                     TContigPos matchEnd,
                     TReadId mateId)
{
    typedef Mater<TExecSpace, TConfig>              TMater;
    typedef typename TMater::TReadSeq               TReadSeq;
    typedef typename TMater::TContig                TContig;
    typedef typename Infix<TContig>::Type           TContigInfix;
    typedef Finder<TContigInfix>                    TFinder;

    TReadId readsCount = length(readSeqs) / 2;

    TContig contig = verifier.contigs[contigId];
    TReadSeq mateSeq = readSeqs[mateId];

    bool reverseComplemented = (mateId >= readsCount);

    TContigPos contigBegin;
    TContigPos contigEnd;

    if (reverseComplemented)
        _getContigInfix(verifier, contigId, matchBegin, matchEnd, contigBegin, contigEnd, LeftMate());
    else
        _getContigInfix(verifier, contigId, matchBegin, matchEnd, contigBegin, contigEnd, RightMate());

    TContigInfix contigInfix = infix(contig, contigBegin, contigEnd);

    TFinder finder(contigInfix);
    setHost(verifier.pattern, mateSeq);

    bool paired = false;
    while (find(finder, verifier.pattern, -static_cast<int>(verifier.readErrors)))
        paired = true;

    return paired;
}

#endif  // #ifndef APP_CUDAMAPPER_VERIFIER_H_
