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

    typedef Myers<>                                                     TAlgorithm;
    typedef Pattern<TReadSeq, TAlgorithm>                               TPattern;

    TOptions const &    options;
    TContigs &          contigs;

    TPattern            pattern;

    unsigned readErrors;
    unsigned long matchesCount;

    Verifier(TOptions const & options, TContigs & contigs) :
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
// Function _getContigInfix()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TContigId, typename TContigPos>
inline void _getContigInfix(Verifier<TExecSpace, TConfig> & verifier,
                            TContigId contigId,
                            TContigPos matchBegin,
                            TContigPos /* matchEnd */,
                            TContigPos & infixBegin,
                            TContigPos & infixEnd,
                            RightMate)
{
    typedef Verifier<TExecSpace, TConfig>           TVerifier;
    typedef typename TVerifier::TContig             TContig;
    typedef typename Size<TContig>::Type            TContigSize;

    TContigSize contigLength = length(verifier.contigs[contigId]);

    infixBegin = 0;
    if (matchBegin + verifier.options.libraryLength > verifier.options.libraryError)
        infixBegin = matchBegin + verifier.options.libraryLength - verifier.options.libraryError;
    infixBegin = _min(infixBegin, contigLength);

    infixEnd = _min(matchBegin + verifier.options.libraryLength + verifier.options.libraryError, contigLength);

    SEQAN_ASSERT_LEQ(infixBegin, infixEnd);
    SEQAN_ASSERT_LEQ(infixEnd - infixBegin, 2 * verifier.options.libraryError);
}

template <typename TExecSpace, typename TConfig, typename TContigId, typename TContigPos>
inline void _getContigInfix(Verifier<TExecSpace, TConfig> & verifier,
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
inline bool findMate(Verifier<TExecSpace, TConfig> & verifier,
                     TReadSeqs & readSeqs,
                     TContigId contigId,
                     TContigPos matchBegin,
                     TContigPos matchEnd,
                     TReadId mateId)
{
    typedef Verifier<TExecSpace, TConfig>           TVerifier;
    typedef typename TVerifier::TReadSeq            TReadSeq;
    typedef typename TVerifier::TContig             TContig;
    typedef typename Infix<TContig>::Type           TContigInfix;
    typedef Finder<TContigInfix>                    TFinder;

    TContig contig = verifier.contigs[contigId];
    TReadSeq mateSeq = readSeqs[mateId];

    TContigPos contigBegin;
    TContigPos contigEnd;

    if (isRevReadSeq(readSeqs, mateId))
        _getContigInfix(verifier, contigId, matchBegin, matchEnd, contigBegin, contigEnd, RightMate());
    else
        _getContigInfix(verifier, contigId, matchBegin, matchEnd, contigBegin, contigEnd, LeftMate());

    TContigInfix contigInfix = infix(contig, contigBegin, contigEnd);

    TFinder finder(contigInfix);
    setHost(verifier.pattern, mateSeq);

    bool paired = false;
    while (find(finder, verifier.pattern, -static_cast<int>(verifier.readErrors)))
        paired = true;

    return paired;
}

// ----------------------------------------------------------------------------
// Function findMate()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename TMatches>
inline void verifyMatches(Verifier<TExecSpace, TConfig> & verifier,
                          TReadSeqs & readSeqs,
                          TMatches const & anchors,
                          TMatches & /* mates */)
{
    typedef Verifier<TExecSpace, TConfig>               TVerifier;
    typedef typename Size<TMatches>::Type               TMatchId;
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    TMatchId matchesCount = length(anchors);

    for (TMatchId matchId = 0; matchId < matchesCount; ++matchId)
    {
        TMatch const & match = anchors[matchId];
        TReadId mateId = getMateSeqId(readSeqs, match.readId);

        if (findMate(verifier, readSeqs, match.contigId, match.contigBegin, match.contigEnd, mateId))
            verifier.matchesCount++;
    }
}

#endif  // #ifndef APP_CUDAMAPPER_VERIFIER_H_
