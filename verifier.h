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

template <typename TOptions_, typename TContigs_, typename TReadSeqs_>
struct VerifierConfig
{
    typedef TOptions_       TOptions;
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
    TContigs const &    contigs;

    TPatternState       patternState;
    TPatternStateRev    patternStateRev;

    Verifier(TOptions const & options, TContigs const & contigs) :
        options(options),
        contigs(contigs)
    {}
};

template <typename TExecSpace, typename TConfig, typename TReadSeq, typename TPos, typename TErrors>
inline void verifyHit(Verifier<TExecSpace, TConfig> & verifier, TReadSeq & readSeq,
                      TPos hitBegin, TPos hitEnd, TErrors hitErrors)
{
//    TContig const & contig = verifier.contigs[getValueI1(hitBegin)];

//    TReadSeqSize readLength = length(readSeq);

}

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename THits, typename TSA>
inline void verifyHits(Verifier<TExecSpace, TConfig> & verifier, TReadSeqs & readSeqs, THits const & hits, TSA const & sa)
{
    typedef Verifier<TExecSpace, TConfig>                               TVerifier;
    typedef typename TVerifier::TContig                                 TContig;
    typedef typename TVerifier::TReadSeq                                TReadSeq;
    typedef typename Size<TReadSeqs>::Type                              TReadId;
    typedef typename Size<THits>::Type                                  THitId;
    typedef typename Size<TContig>::Type                                THitPos;
    typedef typename Value<TSA>::Type                                   THit;

    TReadId pairsCount = length(readSeqs) / 4;

    for (TReadId pairId = 0; pairId < pairsCount; ++pairId)
    {
        // Get mates ids.
        TReadId fwdId = pairId;
        TReadId revId = pairId + 3 * pairsCount;

        // Choose the anchor.
        unsigned long fwdHits = countHits(hits, fwdId);
        unsigned long revHits = countHits(hits, revId);
        unsigned long anchorHits = std::min(fwdHits, revHits);
        TReadId anchorId = (anchorHits == fwdHits) ? fwdId : revId;
        TReadSeq anchor = readSeqs[anchorId];

        // Skip the pair if the anchor is hard.
        if (anchorHits > 300) continue;

        // Consider the hits of all seeds of the anchor.
        THitId hitsBegin = anchorId * (5u + 1);
        THitId hitsEnd = (anchorId + 1) * (5u + 1);
        for (THitId hitId = hitsBegin; hitId < hitsEnd; ++hitId)
        {
            // Verify all hits of a seed of the anchor.
            for (THitPos hitPos = getValueI1(hits.ranges[hitId]); hitPos < getValueI2(hits.ranges[hitId]); ++hitPos)
            {
                THit hit = sa[hitPos];

                verifyHit(verifier, anchor, hit, posAdd(hit, 16u), 0u);
            }
        }
    }
}

#endif  // #ifndef APP_CUDAMAPPER_VERIFIER_H_
