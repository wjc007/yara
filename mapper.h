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

#ifndef APP_CUDAMAPPER_MAPPER_H_
#define APP_CUDAMAPPER_MAPPER_H_

using namespace seqan;


// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Ranges_
// ----------------------------------------------------------------------------

struct Ranges_;

template <typename TSpec = void>
struct Count;

// ----------------------------------------------------------------------------
// Class Hits
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec = void>
struct Hits
{
    typename Member<Hits, Ranges_>::Type    ranges;

    template <typename TFinder>
    inline SEQAN_HOST_DEVICE void
    operator() (TFinder const & finder)
    {
        appendRange(*this, finder);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Member Ranges_
// ----------------------------------------------------------------------------

struct Ranges_;

namespace seqan {
template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, TSpec>, Ranges_>
{
    typedef Pair<typename Size<TIndex>::Type>   TRange_;
    typedef String<TRange_>                     Type;
};

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, Device<TSpec> >, Ranges_>
{
    typedef Pair<typename Size<TIndex>::Type>   TRange_;
    typedef thrust::device_vector<TRange_>      Type;
};
#endif

template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, View<TSpec> >, Ranges_>
{
    typedef typename Member<Hits<TIndex, TSpec>, Ranges_>::Type TRanges_;
    typedef ContainerView<TRanges_, Resizable<TSpec> >          Type;
};
}

// ----------------------------------------------------------------------------
// Member Ranges_; Count
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, Count<TSpec> >, Ranges_>
{
    typedef typename Size<TIndex>::Type         TCount_;
    typedef String<TCount_>                     Type;
};

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, Device<Count<TSpec> > >, Ranges_>
{
    typedef typename Size<TIndex>::Type         TCount_;
    typedef thrust::device_vector<TCount_>      Type;
};
#endif

template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, View<Device<Count<TSpec> > > >, Ranges_>
{
    typedef typename View<typename Member<Hits<TIndex, Device<Count<TSpec> > >, Ranges_>::Type>::Type    Type;
};
}

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TIndex, typename TSpec>
struct View<Hits<TIndex, TSpec> >
{
    typedef Hits<TIndex, View<TSpec> >  Type;
};
}

// ----------------------------------------------------------------------------
// Metafunction Device
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TIndex, typename TSpec>
struct Device<Hits<TIndex, TSpec> >
{
    typedef Hits<TIndex, Device<TSpec> >  Type;
};
}

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
// Metafunction Seeds
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec = void>
struct Seeds
{
    typedef typename Seed<TNeedles>::Type   TSeed_;
    typedef String<TSeed_>                  Type;
};

#ifdef PLATFORM_CUDA
template <typename TNeedles, typename TSpec>
struct Seeds<TNeedles, Device<TSpec> >
{
    typedef typename Seed<TNeedles>::Type   TSeed_;
    typedef thrust::device_vector<TSeed_>   Type;
};
#endif

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
// Function init()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TPattern>
inline void
init(Hits<TIndex, Count<TSpec> > & hits, TPattern const & /* pattern */)
{
    resize(hits.ranges, omp_get_max_threads(), 0, Exact());
}

template <typename TIndex, typename TSpec, typename TPattern>
inline void
init(Hits<TIndex, Device<Count<TSpec> > > & hits, TPattern const & pattern)
{
    resize(hits.ranges, length(needle(pattern)), Exact());
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
// Function appendRange()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TFinder>
inline SEQAN_HOST_DEVICE void
appendRange(Hits<TIndex, Count<TSpec> > & hits, TFinder const & finder)
{
    hits.ranges[getThreadId()] += countOccurrences(textIterator(finder));
}

template <typename TIndex, typename TSpec, typename TFinder>
inline SEQAN_HOST_DEVICE void
appendRange(Hits<TIndex, View<Device<Count<TSpec> > > > & hits, TFinder const & finder)
{
    hits.ranges[getThreadId()] += countOccurrences(textIterator(finder));
}

// ----------------------------------------------------------------------------
// Function getCount()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
inline typename Size<TIndex>::Type
getCount(Hits<TIndex, Count<TSpec> > & hits)
{
    typedef typename Size<TIndex>::Type TSize;

    // TODO(esiragusa): Add function reduce().
    TSize count = 0;
    for (TSize i = 0; i < length(hits.ranges); ++i)
        count += hits.ranges[i];

    return count;
}

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec>
inline typename Size<TIndex>::Type
getCount(Hits<TIndex, Device<Count<TSpec> > > & hits)
{
    return thrust::reduce(hits.ranges.begin(), hits.ranges.end());
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



// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString  genomeFile;
    CharString  genomeIndexFile;
    CharString  readsFile;

    bool        noCuda;
    unsigned    threadsCount;
    int         mappingBlock;
    unsigned    seedLength;
    unsigned    errorsPerSeed;

    Options() :
        noCuda(false),
        threadsCount(1),
        mappingBlock(200000),
        seedLength(33),
        errorsPerSeed(0)
    {}
};

// ----------------------------------------------------------------------------
// Class Mapper
// ----------------------------------------------------------------------------

template <typename TExecSpace>
struct Mapper
{
    typedef Genome<void, CUDAStoreConfig>                           TGenome;
    typedef GenomeLoader<void, CUDAStoreConfig>                     TGenomeLoader;
    typedef Index<typename Contigs<TGenome>::Type, TGenomeIndexSpec> TIndex;
    typedef FragmentStore<void, CUDAStoreConfig>                    TStore;
    typedef ReadsConfig<False, False, True, True, CUDAStoreConfig>  TReadsConfig;
    typedef Reads<void, TReadsConfig>                               TReads;
    typedef ReadsLoader<void, TReadsConfig>                         TReadsLoader;

    TGenome             genome;
#ifdef ENABLE_GENOME_LOADING
    TGenomeLoader       genomeLoader;
#endif
    TIndex              index;
    TStore              store;
    TReads              reads;
    TReadsLoader        readsLoader;

    Timer<double>       timer;

    Mapper() :
        genome(),
#ifdef ENABLE_GENOME_LOADING
        genomeLoader(genome),
#endif
        index(contigs(genome)),
        store(),
        reads(store),
        readsLoader(reads)
    {};
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function configureThreads()
// ----------------------------------------------------------------------------
// Sets the number of threads that OpenMP can spawn.

template <typename TExecSpace>
void configureThreads(Mapper<TExecSpace> & /* mapper */, Options const & options)
{
    omp_set_num_threads(options.threadsCount);
    std::cout << "Threads count:\t\t\t" << omp_get_max_threads() << std::endl;
}

// ----------------------------------------------------------------------------
// Function loadGenome()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void loadGenome(Mapper<TExecSpace> & mapper, Options const & options)
{
    open(mapper.genomeLoader, options.genomeFile);

    std::cout << "Loading genome:\t\t\t" << std::flush;
    start(mapper.timer);
    load(mapper.genomeLoader);
    stop(mapper.timer);
    std::cout << mapper.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function loadGenomeIndex()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void loadGenomeIndex(Mapper<TExecSpace> & mapper, Options const & options)
{
    std::cout << "Loading genome index:\t\t" << std::flush;
    start(mapper.timer);

    if (!open(mapper.index, toCString(options.genomeIndexFile)))
        throw std::runtime_error("Error while opening genome index file.");

    stop(mapper.timer);
    std::cout << mapper.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function loadReads()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void loadReads(Mapper<TExecSpace> & mapper, Options const & options)
{
    std::cout << "Loading reads:\t\t\t" << std::flush;
    start(mapper.timer);
    load(mapper.readsLoader, options.mappingBlock);
    stop(mapper.timer);
    std::cout << mapper.timer << std::endl;

    std::cout << "Reads count:\t\t\t" << mapper.reads.readsCount << std::endl;
}

// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void mapReads(Mapper<TExecSpace> & mapper, Options const & options)
{
    _mapReads(mapper, options, mapper.index, getSeqs(mapper.reads));
}

// ----------------------------------------------------------------------------
// Function _mapReads()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TIndex, typename TReadSeqs>
void _mapReads(Mapper<TExecSpace> & mapper, Options const & options, TIndex & index, TReadSeqs & readSeqs)
{
    typedef typename ExecSpec<TReadSeqs, void>::Type    TSeedsSpec;
    typedef typename Seeds<TReadSeqs, TSeedsSpec>::Type TSeeds;
//    typedef Multiple<Backtracking<HammingDistance> >    TAlgoSpec;
    typedef Multiple<FinderSTree>                       TAlgoSpec;
    typedef Pattern<TSeeds, TAlgoSpec>                  TPattern;
    typedef Finder2<TIndex, TPattern, TAlgoSpec>        TFinder;
    typedef typename ExecSpec<TIndex, Count<> >::Type   THitsSpec;
    typedef Hits<TIndex, THitsSpec>                     THits;

#ifdef PLATFORM_CUDA
    cudaDeviceSynchronize();
#endif

    start(mapper.timer);

//SEQAN_OMP_PRAGMA(critical(_mapper_mapReads_filter))
//{
    // Instantiate a multiple finder.
    TFinder finder(index);
//    setScoreThreshold(finder, options.errorsPerSeed);

    // Collect seeds from read seqs.
    TSeeds seeds;
    fillSeeds(seeds, readSeqs, options.seedLength, TExecSpace());
    std::cout << "Seeds count:\t\t\t" << length(seeds) << std::endl;

    // Instantiate a pattern object.
    TPattern pattern(seeds);

    // Instantiate an object to save the hits.
    THits hits;

    // Resize space for hits.
    init(hits, pattern);

    // Find hits.
    find(finder, pattern, hits);
//}

#ifdef PLATFORM_CUDA
    cudaDeviceSynchronize();
#endif

    stop(mapper.timer);
    std::cout << "Mapping time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Hits count:\t\t\t" << getCount(hits) << std::endl;
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void runMapper(Mapper<TExecSpace> & mapper, Options const & options)
{
#ifdef _OPENMP
    configureThreads(mapper, options);
#endif

#ifdef ENABLE_GENOME_LOADING
    loadGenome(mapper, options);
#endif

    loadGenomeIndex(mapper, options);

    // Open reads file.
    open(mapper.readsLoader, options.readsFile);

    // Reserve space for reads.
    reserve(mapper.reads, options.mappingBlock);

    // Process reads in blocks.
    while (!atEnd(mapper.readsLoader))
    {
        // Load one block of reads.
        loadReads(mapper, options);

        // Map this block of reads.
        mapReads(mapper, options);

        // Clear mapped reads.
        clear(mapper.reads);
    }

    // Close reads file.
    close(mapper.readsLoader);
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

#if false
template <typename TExecSpace>
void runMapper(Mapper<TExecSpace> & mapper, Options const & options)
{
    typedef Timer<double>                               TTimer;
    typedef Logger<std::ostream>                        TLogger;
    typedef Mapper<TExecSpace>                             TMapper;
    typedef Reads<void, typename TMapper::TReadsConfig>    TReads;

    TTimer timer;
    TLogger cout(std::cout);
    TLogger cerr(std::cerr);

#ifdef _OPENMP
    cout << "Threads count:\t\t\t" << options.threadsCount << std::endl;
#endif

#ifdef _OPENMP
    // Disable nested parallelism.
    omp_set_nested(false);
#endif

#ifdef ENABLE_GENOME_LOADING
    loadGenome(mapper, options);
#endif

    loadGenomeIndex(mapper, options);

    // Open reads file.
    open(mapper.readsLoader, options.readsFile);

    // Process reads in parallel.
    SEQAN_OMP_PRAGMA(parallel firstprivate(timer) num_threads(3))
    {
        // Reserve space for reads.
        TReads reads;
        reserve(reads, options.mappingBlock);

        // Process reads.
        while (true)
        {
            // Load a block of reads.
            SEQAN_OMP_PRAGMA(critical(_mapper_readsLoader_load))
            {
                // No more reads.
                if (!atEnd(mapper.readsLoader))
                {
                    start(timer);
                    setReads(mapper.readsLoader, reads);
                    load(mapper.readsLoader, options.mappingBlock);
                    stop(timer);

                    cout << "Loading reads:\t\t\t" << timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;// <<
//                            "Reads count:\t\t\t" << reads.readsCount << "\t\t\t[" << omp_get_thread_num() << "]" << std::endl;
                }
            }

            // No more reads.
            if (!reads.readsCount) break;

            // Map reads.
            SEQAN_OMP_PRAGMA(critical(_mapper_mapReads))
            {
                #ifdef _OPENMP
                // Enable nested parallelism.
                omp_set_nested(true);
                #endif

                #ifdef _OPENMP
                omp_set_num_threads(options.threadsCount);
                #endif

                start(timer);
                mapReads(mapper, options, getSeqs(reads));
                stop(timer);

                cout << "Mapping reads:\t\t\t" << timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;

                #ifdef _OPENMP
                omp_set_num_threads(1);
                #endif

                #ifdef _OPENMP
                // Disable nested parallelism.
                omp_set_nested(false);
                #endif
            }

            // Writer results.
            SEQAN_OMP_PRAGMA(critical(_mapper_samWriter_write))
            {
                start(timer);
                sleep(reads.readsCount / 1000000.0);
                stop(timer);
                
                cout << "Writing results:\t\t" << timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;
            }

            // Clear mapped reads.
            clear(reads);
        }
    }

    // Close reads file.
    close(mapper.readsLoader);
}
#endif

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_H_
