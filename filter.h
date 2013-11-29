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

#ifndef APP_CUDAMAPPER_FILTER_H_
#define APP_CUDAMAPPER_FILTER_H_

using namespace seqan;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Seed
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec = void>
struct Seed
{
    typedef typename Value<TNeedles>::Type  TNeedle_;
    typedef typename View<TNeedle_>::Type   Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TPattern>
inline void
init(Hits<TIndex, TSpec> & hits, TPattern const & pattern)
{
    reserve(hits.ranges, length(needle(pattern)), Exact());
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
typename View<Hits<TIndex, TSpec> >::Type
view(Hits<TIndex, TSpec> & hits)
{
    typename View<Hits<TIndex, TSpec> >::Type hitsView;

    hitsView.ranges = view(hits.ranges);

    return hitsView;
}

// ----------------------------------------------------------------------------
// Function appendRange()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TFinder>
inline void
appendRange(Hits<TIndex, TSpec> & hits, TFinder const & finder)
{
    SEQAN_OMP_PRAGMA(critical(Hits_appendRange))
    appendValue(hits.ranges, range(textIterator(finder)));
}

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec, typename TFinder>
inline SEQAN_HOST_DEVICE void
appendRange(Hits<TIndex, View<Device<TSpec> > > & /* hits */, TFinder const & /* finder */)
{
    // TODO(esiragusa): Global lock.
//    appendValue(hits.ranges, range(textIterator(finder)));
}
#endif

// ----------------------------------------------------------------------------
// Kernel _fillSeedsKernel()
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TSeeds, typename TReadSeqs, typename TSize>
SEQAN_GLOBAL void
_fillSeedsKernel(TSeeds seeds, TReadSeqs readSeqs, TSize seedLength, TSize readSeqsCount, TSize seedsPerReadSeq)
{
    typedef typename Value<TReadSeqs>::Type                 TReadSeq;

    TSize readSeqId = getThreadId();

    if (readSeqId >= readSeqsCount) return;

    TReadSeq const & readSeq = readSeqs[readSeqId];

    for (TSize seedId = 0; seedId < seedsPerReadSeq; ++seedId)
        seeds[readSeqId * seedsPerReadSeq + seedId] = infix(readSeq, seedId * seedLength, (seedId + 1) * seedLength);
}
#endif

// ----------------------------------------------------------------------------
// Function _fillSeeds(); ExecDevice
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TSeeds, typename TReadSeqs, typename TSize>
inline void
_fillSeeds(TSeeds & seeds, TReadSeqs /* const */ & readSeqs,
           TSize seedLength, TSize readSeqsCount, TSize seedsPerReadSeq,
           ExecDevice const & /* tag */)
{
    // Compute grid size.
    unsigned ctaSize = 256;
    unsigned activeBlocks = (readSeqsCount + ctaSize - 1) / ctaSize;

    _fillSeedsKernel<<<activeBlocks, ctaSize>>>(view(seeds), view(readSeqs), seedLength, readSeqsCount, seedsPerReadSeq);
}
#endif

// ----------------------------------------------------------------------------
// Function _fillSeeds(); ExecHost
// ----------------------------------------------------------------------------

template <typename TSeeds, typename TReadSeqs, typename TSize>
inline void
_fillSeeds(TSeeds & seeds, TReadSeqs /* const */ & readSeqs,
           TSize seedLength, TSize readSeqsCount, TSize seedsPerReadSeq,
           ExecHost const & /* tag */)
{
    typedef typename Value<TReadSeqs>::Type                 TReadSeq;
    typedef typename Infix<TReadSeqs>::Type                 TReadSeqInfix;

    for (TSize readSeqId = 0; readSeqId != readSeqsCount; ++readSeqId)
    {
        TReadSeq const & readSeq = readSeqs[readSeqId];

        for (TSize seedId = 0; seedId < seedsPerReadSeq; ++seedId)
        {
            TReadSeqInfix seedInfix = infix(readSeq, seedId * seedLength, (seedId + 1) * seedLength);
            seeds[readSeqId * seedsPerReadSeq + seedId] = view(seedInfix);
        }
    }
}

// ----------------------------------------------------------------------------
// Function fillSeeds()
// ----------------------------------------------------------------------------

template <typename TSeeds, typename TReadSeqs, typename TSeedLength, typename TExecSpace>
inline void
fillSeeds(TSeeds & seeds, TReadSeqs /* const */ & readSeqs, TSeedLength seedLength, TExecSpace const & tag)
{
    typedef typename Size<TReadSeqs>::Type  TSize;

    TSize readSeqsCount = length(readSeqs);
    TSize readSeqLength = length(back(readSeqs));
    TSize seedsPerReadSeq = readSeqLength / seedLength;

    resize(seeds, readSeqsCount * seedsPerReadSeq, Exact());

    _fillSeeds(seeds, readSeqs, static_cast<TSize>(seedLength), readSeqsCount, seedsPerReadSeq, tag);
}

#endif  // #ifndef APP_CUDAMAPPER_FILTER_H_
