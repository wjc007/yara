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
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString          genomeFile;
    CharString          genomeIndexFile;
    Pair<CharString>    readsFile;

    unsigned            errorRate;
    unsigned            libraryLength;
    unsigned            libraryError;

    unsigned            mappingBlock;
    bool                noCuda;
    unsigned            threadsCount;
    unsigned            hitsThreshold;

    Options() :
        errorRate(5),
        libraryLength(220),
        libraryError(50),
        mappingBlock(200000),
        noCuda(false),
        threadsCount(1),
        hitsThreshold(300)
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

    typedef Index<TFMContigs, TGenomeIndexSpec>                     THostIndex;
    typedef typename Space<THostIndex, TExecSpace>::Type            TIndex;

    typedef FragmentStore<void, CUDAStoreConfig>                    TStore;
    typedef ReadsConfig<False, False, True, True, CUDAStoreConfig>  TReadsConfig;
    typedef Reads<PairedEnd, TReadsConfig>                          TReads;
    typedef ReadsLoader<PairedEnd, TReadsConfig>                    TReadsLoader;
    typedef typename TStore::TReadSeqStore                          THostReadSeqs;
    typedef typename Space<THostReadSeqs, TExecSpace>::Type         TReadSeqs;

    typedef typename Size<TIndex>::Type                             TIndexSize;
    typedef Hits<TIndexSize, Exact>                                 THits;

    typedef Match<void>                                             TMatch;
    typedef String<TMatch>                                          TMatches;

    typedef SeederConfig<Options, TIndex, TReadSeqs, Exact>         TSeederConfig;
    typedef Seeder<TExecSpace, TSeederConfig>                       TSeeder;

    typedef ExtenderConfig<Options, TContigs, TReadSeqs, TSeeder>   TExtenderConfig;
    typedef Extender<TExecSpace, TExtenderConfig>                   TExtender;

    typedef VerifierConfig<Options, TContigs, TReadSeqs>            TVerifierConfig;
    typedef Verifier<TExecSpace, TVerifierConfig>                   TVerifier;

//    typedef WriterConfig<Options, TReadSeqs>                        TWriterConfig;
//    typedef Writer<TExecSpace, TWriterConfig>                       TWriter;

    Timer<double>       timer;
    Options const &     options;

    TGenome             genome;
    TGenomeLoader       genomeLoader;
    TIndex              index;
    TStore              store;
    TReads              reads;
    TReadsLoader        readsLoader;

    THits               hits;
    TMatches            anchors;
    TMatches            mates;

    TSeeder             seeder;
    TExtender           extender;
    TVerifier           verifier;
//    TWriter             writer;

    Mapper(Options const & options) :
        options(options),
        genome(),
        genomeLoader(genome),
        index(),
        store(),
        reads(store),
        readsLoader(reads),
        seeder(options, index, 0u),
        extender(options, contigs(genome), seeder),
        verifier(options, contigs(genome))
//        writer(options, genome)
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
void configureThreads(Mapper<TExecSpace> & mapper)
{
#ifdef _OPENMP
    omp_set_num_threads(mapper.options.threadsCount);
    std::cout << "Threads count:\t\t\t" << omp_get_max_threads() << std::endl;
#else
    ignoreUnusedVariableWarning(mapper);
#endif
}

// ----------------------------------------------------------------------------
// Function loadGenome()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void loadGenome(Mapper<TExecSpace> & mapper)
{
    std::cout << "Loading genome:\t\t\t" << std::flush;
    start(mapper.timer);

    CharString genomeFile = mapper.options.genomeIndexFile;
    append(genomeFile, ".txt");

    if (!open(contigs(mapper.genome), toCString(genomeFile)))
        throw RuntimeError("Error while opening genome file.");

    stop(mapper.timer);
    std::cout << mapper.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function loadGenomeIndex()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void loadGenomeIndex(Mapper<TExecSpace> & mapper)
{
#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif

    std::cout << "Loading genome index:\t\t" << std::flush;
    start(mapper.timer);

    if (!open(mapper.index, toCString(mapper.options.genomeIndexFile)))
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
void loadReads(Mapper<TExecSpace> & mapper)
{
    std::cout << "Loading reads:\t\t\t" << std::flush;
    start(mapper.timer);
    load(mapper.readsLoader, mapper.options.mappingBlock);
    stop(mapper.timer);
    std::cout << mapper.timer << std::endl;

    std::cout << "Reads count:\t\t\t" << mapper.reads.readsCount << std::endl;
}

// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void mapReads(Mapper<TExecSpace> & mapper)
{
//SEQAN_OMP_PRAGMA(critical(_mapper_mapReads_filter))
//{
    _mapReads(mapper, getSeqs(mapper.reads));
//}
}

// ----------------------------------------------------------------------------
// Function filterHits()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TReadSeqs>
inline void filterHits(Mapper<TExecSpace> & mapper, TReadSeqs & readSeqs)
{
    typedef Mapper<TExecSpace>                          TMapper;
    typedef typename TMapper::TSeeder                   TSeeder;
    typedef typename TSeeder::TSeedIds                  TSeedIds;
    typedef typename TMapper::TIndexSize                TIndexSize;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    TReadId pairsCount = length(readSeqs) / 4;

    for (TReadId pairId = 0; pairId < pairsCount; ++pairId)
    {
        // Get mates ids.
        TReadId fwdOneId = pairId;
        TReadId fwdTwoId = pairId + pairsCount;
        TReadId revOneId = pairId + 2 * pairsCount;
        TReadId revTwoId = pairId + 3 * pairsCount;

        // Get seed ids.
        TSeedIds fwdOneSeedIds = getSeedIds(mapper.seeder, fwdOneId);
        TSeedIds fwdTwoSeedIds = getSeedIds(mapper.seeder, fwdTwoId);
        TSeedIds revOneSeedIds = getSeedIds(mapper.seeder, revOneId);
        TSeedIds revTwoSeedIds = getSeedIds(mapper.seeder, revTwoId);

        // Count the hits of each read.
        TIndexSize fwdOneHits = countHits(mapper.hits, fwdOneSeedIds);
        TIndexSize fwdTwoHits = countHits(mapper.hits, fwdTwoSeedIds);
        TIndexSize revOneHits = countHits(mapper.hits, revOneSeedIds);
        TIndexSize revTwoHits = countHits(mapper.hits, revTwoSeedIds);

        // Choose the easiest read as the anchor.
        TIndexSize pairOneTwoHits = std::min(fwdOneHits, revTwoHits);
        TIndexSize pairTwoOneHits = std::min(fwdTwoHits, revOneHits);

        // Clear the hits of the mates.
        TSeedIds mateOneTwoSeedIds = (pairOneTwoHits == fwdOneHits) ? revTwoSeedIds : fwdOneSeedIds;
        TSeedIds mateTwoOneSeedIds = (pairTwoOneHits == fwdTwoHits) ? revOneSeedIds : fwdTwoSeedIds;
        clearHits(mapper.hits, mateOneTwoSeedIds);
        clearHits(mapper.hits, mateTwoOneSeedIds);

        // Clear the hits of the anchor and skip the pair.
        if (pairOneTwoHits + pairTwoOneHits > mapper.options.hitsThreshold)
        {
            TSeedIds anchorOneTwoSeedIds = (pairOneTwoHits == fwdOneHits) ? fwdOneSeedIds : revTwoSeedIds;
            TSeedIds anchorTwoOneSeedIds = (pairTwoOneHits == fwdTwoHits) ? fwdTwoSeedIds : revOneSeedIds;
            clearHits(mapper.hits, anchorOneTwoSeedIds);
            clearHits(mapper.hits, anchorTwoOneSeedIds);
        }
    }
}

// ----------------------------------------------------------------------------
// Function _mapReads()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TReadSeqs>
void _mapReads(Mapper<TExecSpace> & mapper, TReadSeqs & readSeqs)
{
    start(mapper.timer);
    clearHits(mapper.hits);
    fillSeeds(mapper.seeder, readSeqs);
    std::cout << "Seeds count:\t\t\t" << length(mapper.seeder.seeds) << std::endl;
    findSeeds(mapper.seeder, mapper.hits);
    stop(mapper.timer);
    std::cout << "Seeding time:\t\t\t" << mapper.timer << std::endl;

    filterHits(mapper, readSeqs);
    std::cout << "Hits count:\t\t\t" << countHits(mapper.hits) << std::endl;

    start(mapper.timer);
    clear(mapper.anchors);
    reserve(mapper.anchors, countHits(mapper.hits) / 5);
    extendHits(mapper.extender, readSeqs, mapper.hits, indexSA(mapper.index), mapper.anchors);
    stop(mapper.timer);
    std::cout << "Extension time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << length(mapper.anchors) << std::endl;

    start(mapper.timer);
    removeDuplicateMatches(mapper.anchors);
    stop(mapper.timer);
    std::cout << "Compaction time:\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << length(mapper.anchors) << std::endl;

    start(mapper.timer);
    mapper.verifier.matchesCount = 0;
    verifyMatches(mapper.verifier, readSeqs, mapper.anchors, mapper.mates);
    stop(mapper.timer);
    std::cout << "Verification time:\t\t" << mapper.timer << std::endl;
    std::cout << "Mates count:\t\t\t" << mapper.verifier.matchesCount << std::endl;
//    std::cout << "Mates count:\t\t" << length(mapper.mates) << std::endl;

//    start(mapper.timer);
//    runWriter(mapper.writer, readSeqs);
//    stop(mapper.timer);
//    std::cout << "Writing time:\t\t\t" << mapper.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void runMapper(Mapper<TExecSpace> & mapper)
{
    configureThreads(mapper);

    printRuler();

    loadGenome(mapper);
    loadGenomeIndex(mapper);

    // Open reads file.
    open(mapper.readsLoader, mapper.options.readsFile);

    // Reserve space for reads.
    reserve(mapper.reads, mapper.options.mappingBlock);

    // Process reads in blocks.
    while (!atEnd(mapper.readsLoader))
    {
        printRuler();

        // Load one block of reads.
        loadReads(mapper);

        // Map this block of reads.
        mapReads(mapper);

        // Clear mapped reads.
        clear(mapper.reads);
    }

    // Close reads file.
    close(mapper.readsLoader);
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void runMapper(Mapper<TExecSpace> & mapper, Parallel)
{
    typedef Mapper<TExecSpace>                          TMapper;
    typedef Reads<void, typename TMapper::TReadsConfig> TReads;
    typedef Logger<std::ostream>                        TLogger;

    Timer<double> timer;
    TLogger cout(std::cout);
    TLogger cerr(std::cerr);

#ifdef _OPENMP
    cout << "Threads count:\t\t\t" << mapper.options.threadsCount << std::endl;
#endif

#ifdef _OPENMP
    // Disable nested parallelism.
    omp_set_nested(false);
#endif

    loadGenome(mapper);
    loadGenomeIndex(mapper);

    // Open reads file.
    open(mapper.readsLoader);

    // Process reads in parallel.
    SEQAN_OMP_PRAGMA(parallel firstprivate(timer) num_threads(3))
    {
        // Reserve space for reads.
        TReads reads;
        reserve(reads, mapper.options.mappingBlock);

        // Process reads.
        while (true)
        {
            // Load a block of reads.
            SEQAN_OMP_PRAGMA(critical(_mapper_readsLoader_load))
            {
                // No more reads.
                if (!atEnd(mapper.readsLoader))
                {
                    start(mapper.timer);
                    setReads(mapper.readsLoader, reads);
                    load(mapper.readsLoader, mapper.options.mappingBlock);
                    stop(mapper.timer);

                    cout << "Loading reads:\t\t\t" << mapper.timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;// <<
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
                omp_set_num_threads(mapper.options.threadsCount);
                #endif

                start(mapper.timer);
                mapReads(mapper, mapper.options, getSeqs(reads));
                stop(mapper.timer);

                cout << "Mapping reads:\t\t\t" << mapper.timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;

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
                start(mapper.timer);
                sleep(reads.readsCount / 1000000.0);
                stop(mapper.timer);
                
                cout << "Writing results:\t\t" << mapper.timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;
            }

            // Clear mapped reads.
            clear(reads);
        }
    }

    // Close reads file.
    close(mapper.readsLoader);
}

// ----------------------------------------------------------------------------
// Function spawnMapper()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void spawnMapper(Options const & options, TExecSpace const & /* tag */)
{
    Mapper<TExecSpace> mapper(options);
    runMapper(mapper);
}

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_H_
