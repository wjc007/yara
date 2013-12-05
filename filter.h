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
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SeederConfig
// ----------------------------------------------------------------------------

template <typename TOptions_, typename TIndex_, typename TReadSeqs_>
struct SeederConfig
{
    typedef TOptions_    TOptions;
    typedef TIndex_      TIndex;
    typedef TReadSeqs_   TReadSeqs;
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
//    typedef typename TConfig::TAlgorithm                    TAlgorithm;

    typedef String<typename Seed<TReadSeqs>::Type>          TSeedsString;
    typedef typename Space<TSeedsString, TExecSpace>::Type  TSeeds;

//    typedef Multiple<Backtracking<HammingDistance> >        TAlgorithm;
    typedef Multiple<FinderSTree>                           TAlgorithm;
    typedef Pattern<TSeeds, TAlgorithm>                     TPattern;
    typedef Finder2<TIndex, TPattern, TAlgorithm>           TFinder;

    TSeeds              seeds;
    TFinder             finder;
    TOptions const &    options;

    Seeder(TIndex & index, TOptions const & options) :
        finder(index),
        options(options)
    {
        // Set the error threshold.
    //    setScoreThreshold(finder, options.errorsPerSeed);
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
// Function _fillSeeds()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs>
inline void _fillSeeds(Seeder<TExecSpace, TConfig> & seeder, TReadSeqs & readSeqs)
{
    typedef typename Size<TReadSeqs>::Type  TSize;

    TSize readSeqsCount = length(readSeqs);
    TSize readSeqLength = length(back(readSeqs));
    TSize readSeqErrors = std::ceil(readSeqLength * (seeder.options.errorRate / 100.0));
    TSize seedLength    = readSeqLength / (readSeqErrors + 1);
    TSize seedsPerReadSeq = readSeqLength / seedLength;

    resize(seeder.seeds, readSeqsCount * seedsPerReadSeq, Exact());

    _fillSeeds(seeder, readSeqs, readSeqsCount, seedsPerReadSeq, seedLength);
}

// ----------------------------------------------------------------------------
// Function runSeeder()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TDelegate>
void runSeeder(Seeder<TExecSpace, TConfig> & seeder, typename TConfig::TReadSeqs & readSeqs, TDelegate & delegate)
{
    typedef Seeder<TExecSpace, TConfig>     TSeeder;
    typedef typename TSeeder::TPattern      TPattern;

    clear(seeder.seeds);

    // Collect seeds from read seqs.
    _fillSeeds(seeder, readSeqs);
    std::cout << "Seeds count:\t\t\t" << length(seeder.seeds) << std::endl;

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif

    // Instantiate a pattern object.
    TPattern pattern(seeder.seeds);

    // Resize space for hits.
    init(delegate, pattern);

    // Find hits.
    find(seeder.finder, pattern, delegate);

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif
}

#endif  // #ifndef APP_CUDAMAPPER_FILTER_H_
