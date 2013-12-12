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

template <typename TExecSpace, typename TConfig = void>
struct Mapper
{
    typedef Genome<void, TConfig>                                   TGenome;
    typedef GenomeLoader<void, TConfig>                             TGenomeLoader;
    typedef typename Contigs<TGenome>::Type                         TContigs;
    typedef typename Value<TContigs>::Type                          TContig;
    typedef typename StringSetPosition<TContigs>::Type              TContigsPos;

    typedef Index<TFMContigs, TGenomeIndexSpec>                     THostIndex;
    typedef typename Space<THostIndex, TExecSpace>::Type            TIndex;

    typedef FragmentStore<void, TConfig>                            TStore;
    typedef ReadsConfig<False, False, True, True, TConfig>          TReadsConfig;
    typedef Reads<PairedEnd, TReadsConfig>                          TReads;
    typedef ReadsLoader<PairedEnd, TReadsConfig>                    TReadsLoader;
    typedef typename TStore::TReadSeqStore                          THostReadSeqs;
    typedef typename Space<THostReadSeqs, TExecSpace>::Type         TReadSeqs;
    typedef typename Value<TReadSeqs>::Type                         TReadSeq;

    typedef typename Size<TIndex>::Type                             TIndexSize;
    typedef Hits<TIndexSize, Exact>                                 THits;

    typedef Match<void>                                             TMatch;
    typedef String<TMatch>                                          TMatches;

    typedef SeederConfig<Options, TIndex, TReadSeqs, Exact>         TSeederConfig;
    typedef Seeder<TExecSpace, TSeederConfig>                       TSeeder;

    typedef AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_> TMyersSpec;
    typedef Myers<TMyersSpec, True, void>                           TExtenderAlgorithm;
    typedef Extender<TContigs, TReadSeq, TExtenderAlgorithm>        TExtender;

    typedef Verifier<TContigs, TReadSeq, Myers<> >                  TVerifier;

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
        extender(contigs(genome)),
        verifier(contigs(genome))
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

template <typename TExecSpace, typename TConfig>
void configureThreads(Mapper<TExecSpace, TConfig> & mapper)
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

template <typename TExecSpace, typename TConfig>
void loadGenome(Mapper<TExecSpace, TConfig> & mapper)
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

template <typename TExecSpace, typename TConfig>
void loadGenomeIndex(Mapper<TExecSpace, TConfig> & mapper)
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

template <typename TExecSpace, typename TConfig>
void loadReads(Mapper<TExecSpace, TConfig> & mapper)
{
    std::cout << "Loading reads:\t\t\t" << std::flush;
    start(mapper.timer);
    load(mapper.readsLoader, mapper.options.mappingBlock);
    stop(mapper.timer);
    std::cout << mapper.timer << std::endl;

    std::cout << "Reads count:\t\t\t" << mapper.reads.readsCount << std::endl;
}

// ----------------------------------------------------------------------------
// Function filterHits()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs>
inline void filterHits(Mapper<TExecSpace, TConfig> & mapper, TReadSeqs & readSeqs)
{
    typedef Mapper<TExecSpace, TConfig>                 TMapper;
    typedef typename TMapper::TSeeder                   TSeeder;
    typedef typename TSeeder::TSeedIds                  TSeedIds;
    typedef typename TMapper::TIndexSize                TIndexSize;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    TReadId pairsCount = getPairsCount(readSeqs);

    for (TReadId pairId = 0; pairId < pairsCount; ++pairId)
    {
        // Get mates ids.
        TReadId firstMateFwdSeqId = getFirstMateFwdSeqId(readSeqs, pairId);
        TReadId secondMateFwdSeqId = getSecondMateFwdSeqId(readSeqs, pairId);
        TReadId firstMateRevSeqId = getFirstMateRevSeqId(readSeqs, pairId);
        TReadId secondMateRevSeqId = getSecondMateRevSeqId(readSeqs, pairId);

        // Get seed ids.
        TSeedIds firstMateFwdSeedIds = getSeedIds(mapper.seeder, firstMateFwdSeqId);
        TSeedIds secondMateFwdSeedIds = getSeedIds(mapper.seeder, secondMateFwdSeqId);
        TSeedIds firstMateRevSeedIds = getSeedIds(mapper.seeder, firstMateRevSeqId);
        TSeedIds secondMateRevSeedIds = getSeedIds(mapper.seeder, secondMateRevSeqId);

        // Count the hits of each read.
        TIndexSize firstMateFwdHits = countHits(mapper.hits, firstMateFwdSeedIds);
        TIndexSize secondMateFwdHits = countHits(mapper.hits, secondMateFwdSeedIds);
        TIndexSize firstMateRevHits = countHits(mapper.hits, firstMateRevSeedIds);
        TIndexSize secondMateRevHits = countHits(mapper.hits, secondMateRevSeedIds);

        // Choose the easiest read as the anchor.
        TIndexSize firstFwdSecondRevHits = std::min(firstMateFwdHits, secondMateRevHits);
        TIndexSize secondFwdFirstRevHits = std::min(secondMateFwdHits, firstMateRevHits);

        // Clear the hits of the mates.
        TSeedIds mateFirstFwdSecondRevSeedIds = (firstFwdSecondRevHits == firstMateFwdHits) ? secondMateRevSeedIds : firstMateFwdSeedIds;
        TSeedIds mateSecondFwdFirstRevSeedIds = (secondFwdFirstRevHits == secondMateFwdHits) ? firstMateRevSeedIds : secondMateFwdSeedIds;
        clearHits(mapper.hits, mateFirstFwdSecondRevSeedIds);
        clearHits(mapper.hits, mateSecondFwdFirstRevSeedIds);

        // Clear the hits of the anchor and skip the pair.
        if (firstFwdSecondRevHits + secondFwdFirstRevHits > mapper.options.hitsThreshold)
        {
            TSeedIds anchorFirstFwdSecondRevSeedIds = (firstFwdSecondRevHits == firstMateFwdHits) ? firstMateFwdSeedIds : secondMateRevSeedIds;
            TSeedIds anchorSecondFwdFirstRevSeedIds = (secondFwdFirstRevHits == secondMateFwdHits) ? secondMateFwdSeedIds : firstMateRevSeedIds;
            clearHits(mapper.hits, anchorFirstFwdSecondRevSeedIds);
            clearHits(mapper.hits, anchorSecondFwdFirstRevSeedIds);
        }
    }
}

// ----------------------------------------------------------------------------
// Function extendHits()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs>
inline void extendHits(Mapper<TExecSpace, TConfig> & mapper, TReadSeqs & readSeqs)
{
    typedef Mapper<TExecSpace, TConfig>                 TMapper;
    typedef typename TMapper::TSeeder                   TSeeder;

    typedef typename TMapper::TContigs                  TContigs;
    typedef typename TMapper::TContig                   TContig;
    typedef typename Size<TContigs>::Type               TContigId;
    typedef typename Size<TContig>::Type                TContigPos;
    typedef typename TMapper::TContigsPos               TContigsPos;

    typedef typename TMapper::TReadSeq                  TReadSeq;
    typedef typename Size<TReadSeqs>::Type              TReadId;
    typedef typename TSeeder::TReadPos                  TReadPos;
    typedef typename TSeeder::TReadSeqSize              TReadSeqSize;

    typedef typename TSeeder::TSeedIds                  TSeedIds;
    typedef typename TSeeder::TSeedId                   TSeedId;

    typedef typename TMapper::THits                     THits;
    typedef typename THits::THitId                      THitId;
    typedef typename THits::THitRange                   THitRange;
    typedef typename THits::THitErrors                  THitErrors;

    typedef typename TMapper::TMatches                  TMatches;
    typedef typename Value<TMatches>::Type              TMatch;

    typedef typename TMapper::TIndex                    TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type       TSA;
    typedef typename Size<TSA>::Type                    TSAPos;
    typedef typename Value<TSA>::Type                   TSAValue;

    TSA & sa = indexSA(mapper.index);

    TReadId readSeqsCount = getReadSeqsCount(readSeqs);
    for (TReadId readSeqId = 0; readSeqId < readSeqsCount; ++readSeqId)
    {
        TReadSeq readSeq = readSeqs[readSeqId];

        TSeedIds seedIds = getSeedIds(mapper.seeder, readSeqId);

        for (TSeedId seedId = getValueI1(seedIds); seedId < getValueI2(seedIds); ++seedId)
        {
            // Get position in read.
            TReadPos readPos = getPosInRead(mapper.seeder, seedId);
            TReadSeqSize seedLength = getValueI2(readPos) - getValueI1(readPos);

            // TODO(esiragusa): iterate over all hits of the seed.
            // THitIds hitIds = getHitIds(mapper.seeder, seedId);
            {
                THitId hitId = seedId;

                THitRange hitRange = getHitRange(mapper.hits, hitId);
                THitErrors hitErrors = getHitErrors(mapper.hits, hitId);

                for (TSAPos saPos = getValueI1(hitRange); saPos < getValueI2(hitRange); ++saPos)
                {
                    // Invert SA value.
                    TSAValue saValue = sa[saPos];
                    setSeqOffset(saValue, suffixLength(saValue, contigs(mapper.genome)) - seedLength);

                    // Compute position in contig.
                    TContigsPos contigBegin = TContigsPos(getValueI1(saValue), getValueI2(saValue));
                    TContigsPos contigEnd = contigBegin;
                    posAdd(contigEnd, seedLength);

                    extend(mapper.extender,
                           readSeq,
                           contigBegin, contigEnd,
                           readPos.i1, readPos.i2,
                           hitErrors, mapper.options.errorRate,
                           mapper.anchors);
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function _getMateContigPos()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TContigPos, typename TMatch>
inline void _getMateContigPos(Mapper<TExecSpace, TConfig> & mapper,
                              TContigPos & contigBegin,
                              TContigPos & contigEnd,
                              TMatch const & anchor,
                              RightMate)
{
    typedef Mapper<TExecSpace, TConfig>             TMapper;
    typedef typename TMapper::TContig               TContig;
    typedef typename Size<TContig>::Type            TContigSize;

    TContigSize contigLength = length(contigs(mapper.genome)[anchor.contigId]);

    setValueI1(contigBegin, anchor.contigId);
    setValueI1(contigEnd, anchor.contigId);

    contigBegin.i2 = 0;
    if (anchor.contigBegin + mapper.options.libraryLength > mapper.options.libraryError)
        contigBegin.i2 = anchor.contigBegin + mapper.options.libraryLength - mapper.options.libraryError;
    contigBegin.i2 = _min(contigBegin.i2, contigLength);

    contigEnd.i2 = _min(anchor.contigBegin + mapper.options.libraryLength + mapper.options.libraryError, contigLength);

    SEQAN_ASSERT_LEQ(getValueI2(contigBegin), getValueI2(contigEnd));
    SEQAN_ASSERT_LEQ(getValueI2(contigEnd) - getValueI2(contigBegin), 2 * mapper.options.libraryError);
}

template <typename TExecSpace, typename TConfig, typename TContigPos, typename TMatch>
inline void _getMateContigPos(Mapper<TExecSpace, TConfig> & mapper,
                              TContigPos & contigBegin,
                              TContigPos & contigEnd,
                              TMatch const & anchor,
                              LeftMate)
{
    setValueI1(contigBegin, anchor.contigId);
    setValueI1(contigEnd, anchor.contigId);

    contigBegin.i2 = 0;
    if (anchor.contigEnd > mapper.options.libraryLength + mapper.options.libraryError)
        contigBegin.i2 = anchor.contigEnd - mapper.options.libraryLength - mapper.options.libraryError;

    contigEnd.i2 = 0;
    if (anchor.contigEnd + mapper.options.libraryError > mapper.options.libraryLength)
        contigEnd.i2 = anchor.contigEnd - mapper.options.libraryLength + mapper.options.libraryError;

    SEQAN_ASSERT_LEQ(getValueI2(contigBegin), getValueI2(contigEnd));
    SEQAN_ASSERT_LEQ(getValueI2(contigEnd) - getValueI2(contigBegin), 2 * mapper.options.libraryError);
}

// ----------------------------------------------------------------------------
// Function findMates()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename TMatches>
inline void findMates(Mapper<TExecSpace, TConfig> & mapper,
                      TReadSeqs & readSeqs,
                      TMatches const & anchors,
                      TMatches & mates)
{
    typedef Mapper<TExecSpace, TConfig>                 TMapper;
    typedef typename TMapper::TContigsPos               TContigsPos;

    typedef typename Size<TMatches>::Type               TMatchId;
    typedef typename Value<TMatches>::Type              TMatch;

    typedef typename Size<TReadSeqs>::Type              TReadId;
    typedef typename Value<TReadSeqs>::Type             TReadSeq;

    TMatchId matchesCount = length(anchors);

    for (TMatchId matchId = 0; matchId < matchesCount; ++matchId)
    {
        TMatch const & match = anchors[matchId];
        TReadId mateId = getMateSeqId(readSeqs, match.readId);
        TReadSeq mateSeq = readSeqs[mateId];

        TContigsPos contigBegin;
        TContigsPos contigEnd;

        if (isRevReadSeq(readSeqs, mateId))
            _getMateContigPos(mapper, contigBegin, contigEnd, match, RightMate());
        else
            _getMateContigPos(mapper, contigBegin, contigEnd, match, LeftMate());

        // TODO(esiragusa): convert this to errorsPerRead.
        verify(mapper.verifier, mateSeq, contigBegin, contigEnd, mapper.options.errorRate, mates);
    }
}

// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig>
void mapReads(Mapper<TExecSpace, TConfig> & mapper)
{
//SEQAN_OMP_PRAGMA(critical(_mapper_mapReads_filter))
//{
    _mapReads(mapper, getSeqs(mapper.reads));
//}
}

// ----------------------------------------------------------------------------
// Function _mapReads()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs>
void _mapReads(Mapper<TExecSpace, TConfig> & mapper, TReadSeqs & readSeqs)
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
    extendHits(mapper, readSeqs);
    stop(mapper.timer);
    std::cout << "Extension time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << length(mapper.anchors) << std::endl;

    start(mapper.timer);
    removeDuplicateMatches(mapper.anchors);
    stop(mapper.timer);
    std::cout << "Compaction time:\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << length(mapper.anchors) << std::endl;
    std::cout << "Anchored pairs:\t\t\t" << countMatches(readSeqs, mapper.anchors, PairedEnd()) << std::endl;

    start(mapper.timer);
    clear(mapper.mates);
    reserve(mapper.mates, length(mapper.anchors), Exact());
    findMates(mapper, readSeqs, mapper.anchors, mapper.mates);
    stop(mapper.timer);
    std::cout << "Verification time:\t\t" << mapper.timer << std::endl;
    std::cout << "Mates count:\t\t\t" << length(mapper.mates) << std::endl;
    std::cout << "Mapped pairs:\t\t\t" << countMatches(readSeqs, mapper.mates, PairedEnd()) << std::endl;

//    start(mapper.timer);
//    runWriter(mapper.writer, readSeqs);
//    stop(mapper.timer);
//    std::cout << "Writing time:\t\t\t" << mapper.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig>
void runMapper(Mapper<TExecSpace, TConfig> & mapper)
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

template <typename TExecSpace, typename TConfig>
void runMapper(Mapper<TExecSpace, TConfig> & mapper, Parallel)
{
    typedef Mapper<TExecSpace, TConfig>                          TMapper;
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
    Mapper<TExecSpace, CUDAStoreConfig> mapper(options);
    runMapper(mapper);
}

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_H_
