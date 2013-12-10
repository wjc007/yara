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

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SeederConfig
// ----------------------------------------------------------------------------

template <typename TOptions_, typename TIndex_, typename TReadSeqs_, typename TDistance_>
struct SeederConfig
{
    typedef TOptions_    TOptions;
    typedef TIndex_      TIndex;
    typedef TReadSeqs_   TReadSeqs;
    typedef TDistance_   TDistance;
};

// ----------------------------------------------------------------------------
// Class Seeder
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig>
struct Seeder
{
    typedef typename TConfig::TOptions                      TOptions;
    typedef typename TConfig::TIndex                        TIndex;
    typedef typename TConfig::TReadSeqs                     TReadSeqs;
    typedef typename TConfig::TDistance                     TDistance;

    typedef String<typename Seed<TReadSeqs>::Type>          THostSeeds;
    typedef typename Space<THostSeeds, TExecSpace>::Type    TSeeds;
    typedef typename Size<TSeeds>::Type                     TSeedId;
    typedef Pair<TSeedId>                                   TSeedIds;
    typedef typename Size<TReadSeqs>::Type                  TReadSeqSize;
    typedef Pair<TReadSeqSize>                              TReadPos;

    typedef typename If<IsSameType<TDistance, Exact>,
                        FinderSTree,
                        Backtracking<HammingDistance>
                       >::Type                              TAlgorithmSpec;
    typedef Multiple<TAlgorithmSpec>                        TAlgorithm;
    typedef Pattern<TSeeds, TAlgorithm>                     TPattern;
    typedef Finder2<TIndex, TPattern, TAlgorithm>           TFinder;

    TOptions const &    options;
    TFinder             finder;
    TSeeds              seeds;
    TReadSeqSize        errorsPerSeed;

    TSeedId             readsCount;
    TReadSeqSize        readsLength;
    TReadSeqSize        seedsLength;
    TReadSeqSize        seedsPerRead;
    TReadSeqSize        errorsPerRead;

    template <typename TErrors>
    Seeder(TOptions const & options, TIndex & index, TErrors errorsPerSeed) :
        options(options),
        finder(index),
        errorsPerSeed(errorsPerSeed),
        readsCount(0),
        readsLength(0),
        seedsLength(0),
        seedsPerRead(0),
        errorsPerRead(0)
    {
        // Set the error threshold.
    //    setScoreThreshold(finder, errorsPerSeed);
    }
};

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
inline void _fillSeeds(Seeder<ExecDevice, TConfig> & seeder, TReadSeqs & readSeqs,
                       TSize readSeqsCount, TSize seedsPerReadSeq, TSize seedLength)
{
    // Compute grid size.
    unsigned ctaSize = 256;
    unsigned activeBlocks = (readSeqsCount + ctaSize - 1) / ctaSize;

    _fillSeedsKernel<<<activeBlocks, ctaSize>>>(view(seeder.seeds), view(readSeqs),
                                                readSeqsCount, seedsPerReadSeq, seedLength);
}
#endif

// ----------------------------------------------------------------------------
// Function _fillSeeds(); ExecHost
// ----------------------------------------------------------------------------

template <typename TConfig, typename TReadSeqs, typename TSize>
inline void _fillSeeds(Seeder<ExecHost, TConfig> & seeder, TReadSeqs & readSeqs,
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

            seeder.seeds[readSeqId * seedsPerReadSeq + seedId] = view(seedInfix);
        }
    }
}

// ----------------------------------------------------------------------------
// Function fillSeeds()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs>
inline void fillSeeds(Seeder<TExecSpace, TConfig> & seeder, TReadSeqs & readSeqs)
{
    seeder.readsCount = length(readSeqs);
    seeder.readsLength = length(front(readSeqs));
    seeder.errorsPerRead = std::ceil(seeder.readsLength * (seeder.options.errorRate / 100.0));
    seeder.seedsPerRead = seeder.errorsPerRead + 1;
    seeder.seedsLength = seeder.readsLength / seeder.seedsPerRead;

    // Resize space for seeds.
    clear(seeder.seeds);
    resize(seeder.seeds, seeder.readsCount * seeder.seedsPerRead, Exact());

    _fillSeeds(seeder, readSeqs, seeder.readsCount, seeder.seedsPerRead, seeder.seedsLength);
}

// ----------------------------------------------------------------------------
// Function findSeeds()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TDelegate>
inline void findSeeds(Seeder<TExecSpace, TConfig> & seeder, TDelegate & delegate)
{
    typedef Seeder<TExecSpace, TConfig>     TSeeder;
    typedef typename TSeeder::TPattern      TPattern;

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif

    // Instantiate a pattern object.
    TPattern pattern(seeder.seeds);

    // Initialize the delegate.
    init(delegate, pattern);

    // Find hits.
    find(seeder.finder, pattern, delegate);

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif
}

// ----------------------------------------------------------------------------
// Function getSeedIds()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadId>
inline typename Seeder<TExecSpace, TConfig>::TSeedIds
getSeedIds(Seeder<TExecSpace, TConfig> & seeder, TReadId readId)
{
    typedef Seeder<TExecSpace, TConfig>     TSeeder;
    typedef typename TSeeder::TSeedIds      TSeedIds;

    return TSeedIds(readId * seeder.seedsPerRead, (readId + 1) * seeder.seedsPerRead);
}

// ----------------------------------------------------------------------------
// Function getLocalSeedId()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TSeedId>
inline typename Seeder<TExecSpace, TConfig>::TSeedId
getLocalSeedId(Seeder<TExecSpace, TConfig> & seeder, TSeedId seedId)
{
    return seedId % seeder.seedsPerRead;
}

// ----------------------------------------------------------------------------
// Function getReadId()
// ----------------------------------------------------------------------------

//template <typename TExecSpace, typename TConfig, typename TSeedId>
//inline typename Seeder<TExecSpace, TConfig>::TReadId
//getReadId(Seeder<TExecSpace, TConfig> & seeder, TSeedId seedId)
//{
//    return ...;
//}

// ----------------------------------------------------------------------------
// Function getPosInRead()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TSeedId>
inline typename Seeder<TExecSpace, TConfig>::TReadPos
getPosInRead(Seeder<TExecSpace, TConfig> & seeder, TSeedId seedId)
{
    typedef Seeder<TExecSpace, TConfig>     TSeeder;
    typedef typename TSeeder::TReadPos      TReadPos;

    TSeedId localSeedId = getLocalSeedId(seeder, seedId);

    return TReadPos(localSeedId * seeder.seedsLength, (localSeedId + 1) * seeder.seedsLength);
}

#endif  // #ifndef APP_CUDAMAPPER_SEEDER_H_
