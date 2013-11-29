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
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Space
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec = void>
struct Space
{
    typedef TObject Type;
};

template <typename TObject>
struct Space<TObject, ExecDevice>
{
    typedef typename Device<TObject>::Type  Type;
};

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

    typedef typename Contigs<TGenome>::Type                         TContigs;
    typedef Index<TContigs, TGenomeIndexSpec>                       THostIndex;
    typedef typename Space<THostIndex, TExecSpace>::Type            TIndex;

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
        index(),
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
#ifdef _OPENMP
    omp_set_num_threads(options.threadsCount);
    std::cout << "Threads count:\t\t\t" << omp_get_max_threads() << std::endl;
#else
    ignoreUnusedVariableWarning(options);
#endif
}

// ----------------------------------------------------------------------------
// Function loadGenome()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void loadGenome(Mapper<TExecSpace> & mapper, Options const & options)
{
#ifdef ENABLE_GENOME_LOADING
    open(mapper.genomeLoader, options.genomeFile);

    std::cout << "Loading genome:\t\t\t" << std::flush;
    start(mapper.timer);
    load(mapper.genomeLoader);
    stop(mapper.timer);
    std::cout << mapper.timer << std::endl;
#else
    ignoreUnusedVariableWarning(mapper);
    ignoreUnusedVariableWarning(options);
#endif
}

// ----------------------------------------------------------------------------
// Function loadGenomeIndex()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void loadGenomeIndex(Mapper<TExecSpace> & mapper, Options const & options)
{
#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif

    std::cout << "Loading genome index:\t\t" << std::flush;
    start(mapper.timer);

    if (!open(mapper.index, toCString(options.genomeIndexFile)))
        throw RuntimeError("Error while opening genome index file.");

    stop(mapper.timer);
    std::cout << mapper.timer << std::endl;

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif
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
    _mapReads(mapper, options, getSeqs(mapper.reads));
}

// ----------------------------------------------------------------------------
// Function _mapReads()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TReadSeqs>
void _mapReads(Mapper<TExecSpace> & mapper, Options const & options, TReadSeqs & readSeqs)
{
    typedef Mapper<TExecSpace>                              TMapper;
    typedef typename TMapper::TIndex                        TIndex;

    typedef String<typename Seed<TReadSeqs>::Type>          TSeedsString;
    typedef typename Space<TSeedsString, TExecSpace>::Type  TSeeds;

//    typedef Multiple<Backtracking<HammingDistance> >    TAlgoSpec;
    typedef Multiple<FinderSTree>                       TAlgoSpec;
    typedef Pattern<TSeeds, TAlgoSpec>                  TPattern;
    typedef Finder2<TIndex, TPattern, TAlgoSpec>        TFinder;
    typedef OccurrencesCounter<TIndex>                  TCounter;

#ifdef PLATFORM_CUDA
    cudaDeviceSynchronize();
#endif

    start(mapper.timer);

//SEQAN_OMP_PRAGMA(critical(_mapper_mapReads_filter))
//{
    // Instantiate a multiple finder.
    TFinder finder(mapper.index);
//    setScoreThreshold(finder, options.errorsPerSeed);

    // Collect seeds from read seqs.
    TSeeds seeds;
    fillSeeds(seeds, readSeqs, options.seedLength, TExecSpace());
    std::cout << "Seeds count:\t\t\t" << length(seeds) << std::endl;

    // Instantiate a pattern object.
    TPattern pattern(seeds);

    // Instantiate an object to save the hits.
    TCounter hits;

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
    configureThreads(mapper, options);

    loadGenome(mapper, options);
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

    loadGenome(mapper, options);
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

// ----------------------------------------------------------------------------
// Function spawnMapper()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void spawnMapper(Options const & options, TExecSpace const & /* tag */)
{
    Mapper<TExecSpace> mapper;
    runMapper(mapper, options);
}

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_H_
