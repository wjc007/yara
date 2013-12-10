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

template <typename TOptions_, typename TContigs_, typename TReadSeqs_, typename TSeeder_>
struct ExtenderConfig
{
    typedef TOptions_       TOptions;
    typedef TContigs_       TContigs;
    typedef TReadSeqs_      TReadSeqs;
    typedef TSeeder_        TSeeder;
};

// ----------------------------------------------------------------------------
// Class Extender
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig>
struct Extender
{
    typedef typename TConfig::TOptions                                  TOptions;
    typedef typename TConfig::TContigs                                  TContigs;
    typedef typename TConfig::TReadSeqs                                 TReadSeqs;
    typedef typename TConfig::TSeeder                                   TSeeder;

    typedef typename Value<TContigs>::Type                              TContig;
    typedef typename Value<TReadSeqs>::Type                             TReadSeq;
    typedef typename Infix<TReadSeq>::Type                              TReadInfix;
    typedef ModifiedString<TReadInfix, ModReverse>                      TReadInfixRev;

    typedef AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_>   TMyersSpec;
    typedef Myers<TMyersSpec, True, void>                               TAlgorithm;
    typedef PatternState_<TReadInfix, TAlgorithm>                       TPatternState;
    typedef PatternState_<TReadInfixRev, TAlgorithm>                    TPatternStateRev;

    TOptions const &    options;
    TContigs &          contigs;
    TSeeder &           seeder;

    TPatternState       patternState;
    TPatternStateRev    patternStateRev;

    unsigned readErrors;
    unsigned matchesCount;

    Extender(TOptions const & options, TContigs & contigs, TSeeder & seeder) :
        options(options),
        contigs(contigs),
        seeder(seeder),
        readErrors(5),
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
// Function extendHit()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs,
          typename TReadId, typename TReadPos, typename TContigId, typename TContigPos, typename TErrors>
inline bool extendHit(Extender<TExecSpace, TConfig> & extender,
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

        TContigInfix contigRight = infix(contig, contigEnd, contigRightEnd);
        TReadInfix readRight = infix(readSeq, readEnd, readLength);

        if (!_extendRight(extender, extender.patternState, contigRight, readRight, matchEnd, readErrors))
            return false;
    }

    return true;
}

// ----------------------------------------------------------------------------
// Function extendHits()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename THits, typename TSA, typename TMatches>
inline void extendHits(Extender<TExecSpace, TConfig> & extender, TReadSeqs & readSeqs, THits const & hits, TSA const & sa, TMatches & matches)
{
    typedef Extender<TExecSpace, TConfig>               TExtender;
    typedef typename TExtender::TSeeder                 TSeeder;
    typedef typename TExtender::TContigs                TContigs;
    typedef typename TExtender::TContig                 TContig;
    typedef typename Size<TContigs>::Type               TContigId;
    typedef typename Size<TContig>::Type                TContigPos;
    typedef typename Size<TReadSeqs>::Type              TReadId;
    typedef typename TSeeder::TReadPos                  TReadPos;
    typedef typename TSeeder::TReadSeqSize              TReadSeqSize;
    typedef typename TSeeder::TSeedIds                  TSeedIds;
    typedef typename TSeeder::TSeedId                   TSeedId;
    typedef typename THits::THitId                      THitId;
    typedef typename THits::THitRange                   THitRange;
    typedef typename THits::THitErrors                  THitErrors;
    typedef typename Size<TSA>::Type                    TSAPos;
    typedef typename Value<TSA>::Type                   TSAValue;

    TReadId readsCount = length(readSeqs);
    for (TReadId readId = 0; readId < readsCount; ++readId)
    {
        TSeedIds seedIds = getSeedIds(extender.seeder, readId);

        for (TSeedId seedId = getValueI1(seedIds); seedId < getValueI2(seedIds); ++seedId)
        {
            // Get position in read.
            TReadPos readPos = getPosInRead(extender.seeder, seedId);
            TReadSeqSize seedLength = getValueI2(readPos) - getValueI1(readPos);

            // TODO(esiragusa): iterate over all hits of the seed.
            // THitIds hitIds = getHitIds(extender.seeder, seedId);
            {
                THitId hitId = seedId;

                THitRange hitRange = getHitRange(hits, hitId);
                THitErrors hitErrors = getHitErrors(hits, hitId);

                for (TSAPos saPos = getValueI1(hitRange); saPos < getValueI2(hitRange); ++saPos)
                {
                    // Invert SA value.
                    TSAValue saValue = sa[saPos];
                    setSeqOffset(saValue, suffixLength(saValue, extender.contigs) - seedLength);

                    // Compute position in contig.
                    TContigId contigId = getValueI1(saValue);
                    TContigPos contigBegin = getValueI2(saValue);
                    TContigPos contigEnd = getValueI2(saValue) + seedLength;

                    if (extendHit(extender, readSeqs,
                                  readId, getValueI1(readPos), getValueI2(readPos),
                                  contigId, contigBegin, contigEnd,
                                  hitErrors))
                    {
                        extender.matchesCount++;
//                        appendValue(matches, match);
                    }
                }
            }
        }
    }
}

#endif  // #ifndef APP_CUDAMAPPER_EXTENDER_H_
