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

#ifndef APP_CUDAMAPPER_SEEDER_H_
#define APP_CUDAMAPPER_SEEDER_H_

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

// ----------------------------------------------------------------------------
// Metafunction SeedingAlgorithm
// ----------------------------------------------------------------------------

template <typename TDistance, typename TSpec = void>
struct SeedingAlgorithm
{
    typedef typename If<IsSameType<TDistance, Exact>,
                        FinderSTree, Backtracking<HammingDistance>
                       >::Type  Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SeedSet
// ----------------------------------------------------------------------------

template <typename TReadSeqs, typename TSpec = void>
struct SeedSet
{
    typedef String<typename Seed<TReadSeqs>::Type>          THostSeeds;
    typedef typename Space<THostSeeds, TExecSpace>::Type    TSeeds;

    typedef typename Size<TSeeds>::Type                     TSeedId;
    typedef Pair<TSeedId>                                   TSeedIds;

    typedef typename Size<TReadSeqs>::Type                  TReadSeqId;
    typedef typename Value<TReadSeqs>::Type                 TReadSeq;
    typedef typename Size<TReadSeq>::Type                   TReadSeqSize;
    typedef Pair<TReadSeqSize>                              TReadPos;

    TSeeds              seeds;

    TReadSeqSize        errorsPerSeed;
    TSeedId             readsCount;
    TReadSeqSize        readsLength;
    TReadSeqSize        seedsLength;
    TReadSeqSize        seedsPerRead;
    TReadSeqSize        errorsPerRead;

    SeedSet() :
        errorsPerSeed(0),
        readsCount(0),
        readsLength(0),
        seedsLength(0),
        seedsPerRead(0),
        errorsPerRead(0)
    {}
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TReadSeqs, typename TSpec>
struct Value<SeedSet<TReadSeqs, TSpec> >
{
    typedef Type;
};
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Kernel _fillSeedsKernel()
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TSeeds, typename TReadSeqs, typename TSize>
SEQAN_GLOBAL void
_fillSeedsKernel(TSeeds seeds, TReadSeqs readSeqs, TSize readSeqsCount, TSize seedsPerReadSeq, TSize seedLength)
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
template <typename TConfig, typename TReadSeqs, typename TSize>
inline void _fillSeeds(SeedSet<TReadSeqs, TSpec> & seedSet, TReadSeqs & readSeqs,
                       TSize readSeqsCount, TSize seedsPerReadSeq, TSize seedLength)
{
    // Compute grid size.
    unsigned ctaSize = 256;
    unsigned activeBlocks = (readSeqsCount + ctaSize - 1) / ctaSize;

    _fillSeedsKernel<<<activeBlocks, ctaSize>>>(view(seedSet.seeds), view(readSeqs),
                                                readSeqsCount, seedsPerReadSeq, seedLength);
}
#endif

// ----------------------------------------------------------------------------
// Function _fillSeeds(); ExecHost
// ----------------------------------------------------------------------------

template <typename TConfig, typename TReadSeqs, typename TSize>
inline void _fillSeeds(SeedSet<ExecHost, TConfig> & seedSet, TReadSeqs & readSeqs,
                       TSize readSeqsCount, TSize seedsPerReadSeq, TSize seedLength)
{
    typedef typename Value<TReadSeqs>::Type                 TReadSeq;
    typedef typename Infix<TReadSeqs>::Type                 TReadSeqInfix;

    for (TSize readSeqId = 0; readSeqId != readSeqsCount; ++readSeqId)
    {
        TReadSeq const & readSeq = readSeqs[readSeqId];

        for (TSize seedId = 0; seedId < seedsPerReadSeq; ++seedId)
        {
            TReadSeqInfix seedInfix = infix(readSeq, seedId * seedLength, (seedId + 1) * seedLength);

            seedSet.seeds[readSeqId * seedsPerReadSeq + seedId] = view(seedInfix);
        }
    }
}

// ----------------------------------------------------------------------------
// Function fill()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs>
inline void fill(SeedSet<TExecSpace, TConfig> & seedSet, TReadSeqs & readSeqs)
{
    seedSet.readsCount = length(readSeqs);
    seedSet.readsLength = length(front(readSeqs));
    seedSet.errorsPerRead = std::ceil(seedSet.readsLength * (seedSet.options.errorRate / 100.0));
    seedSet.seedsPerRead = seedSet.errorsPerRead + 1;
    seedSet.seedsLength = seedSet.readsLength / seedSet.seedsPerRead;

    // Resize space for seeds.
    clear(seedSet.seeds);
    resize(seedSet.seeds, seedSet.readsCount * seedSet.seedsPerRead, Exact());

    _fillSeeds(seedSet, readSeqs, seedSet.readsCount, seedSet.seedsPerRead, seedSet.seedsLength);
}

// ----------------------------------------------------------------------------
// Function getSeedIds()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadId>
inline typename SeedSet<TExecSpace, TConfig>::TSeedIds
getSeedIds(SeedSet<TExecSpace, TConfig> const & seedSet, TReadId readId)
{
    typedef SeedSet<TExecSpace, TConfig>     TSeedSet;
    typedef typename TSeedSet::TSeedIds      TSeedIds;

    return TSeedIds(readId * seedSet.seedsPerRead, (readId + 1) * seedSet.seedsPerRead);
}

// ----------------------------------------------------------------------------
// Function getReadSeqId()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TSeedId>
inline typename SeedSet<TExecSpace, TConfig>::TReadSeqId
getReadSeqId(SeedSet<TExecSpace, TConfig> const & seedSet, TSeedId seedId)
{
    return seedId / seedSet.seedsPerRead;
}

// ----------------------------------------------------------------------------
// Function getLocalSeedId()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TSeedId>
inline typename SeedSet<TExecSpace, TConfig>::TSeedId
getLocalSeedId(SeedSet<TExecSpace, TConfig> const & seedSet, TSeedId seedId)
{
    return seedId % seedSet.seedsPerRead;
}

// ----------------------------------------------------------------------------
// Function getPosInRead()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TSeedId>
inline typename SeedSet<TExecSpace, TConfig>::TReadPos
getPosInRead(SeedSet<TExecSpace, TConfig> const & seedSet, TSeedId seedId)
{
    typedef SeedSet<TExecSpace, TConfig>     TSeedSet;
    typedef typename TSeedSet::TReadPos      TReadPos;

    TSeedId localSeedId = getLocalSeedId(seedSet, seedId);

    return TReadPos(localSeedId * seedSet.seedsLength, (localSeedId + 1) * seedSet.seedsLength);
}

#endif  // #ifndef APP_CUDAMAPPER_SEEDER_H_
