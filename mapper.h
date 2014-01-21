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

    typedef StringSet<TReadSeqs, Segment<TReadSeqs> >               TSeedsSet;
    typedef Tuple<TSeedsSet, 3>                                     TSeeds;

    typedef Multiple<FinderSTree>                                   TSeedingExt;
    typedef Multiple<Backtracking<HammingDistance> >                TSeedingApx;
    typedef Pattern<TSeedsSet, TSeedingExt>                         TSeedsExt;
    typedef Pattern<TSeedsSet, TSeedingApx>                         TSeedsApx;
    typedef Finder2<TIndex, TSeedsExt, TSeedingExt>                 TSeederExt;
    typedef Finder2<TIndex, TSeedsApx, TSeedingApx>                 TSeederApx;

    typedef typename Size<TIndex>::Type                             TIndexSize;
    typedef Hit<TIndexSize, HammingDistance>                        THit;
    typedef String<THit>                                            THitsString;
    typedef Tuple<THitsString, 3>                                   THits;

    typedef Match<void>                                             TMatch;
    typedef String<TMatch>                                          TMatches;

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

    TSeeds              seeds;
    THits               hits;
    TMatches            anchors;
    TMatches            mates;

    TSeederExt          seederExt;
    TSeederApx          seederApx;
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
//        seeds(getSeqs(reads)),
        hits(),
        anchors(),
        mates(),
        seederExt(index),
        seederApx(index),
        extender(contigs(genome)),
        verifier(contigs(genome))
//        writer(options, genome)
    {
        setHost(seeds[0], getSeqs(reads));
        // Set the error threshold.
    //    setScoreThreshold(seederApx, errorsPerSeed);
    };
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
// Function selectSeeds()
// ----------------------------------------------------------------------------

template <typename TReadSeqs, typename TReadSeqId, typename TDelegate>
inline void _selectSeeds(TReadSeqs const & readSeqs, TReadSeqId readSeqId, TDelegate & delegate)
{
    typedef typename StringSetPosition<TReadSeqs>::Type     TPos;
    typedef typename Value<TReadSeqs>::Type                 TReadSeq;
    typedef typename Size<TReadSeq>::Type                   TSize;

    TSize readLength = length(readSeqs[readSeqId]);
//    TSize errorsPerRead = std::ceil(readsLength * (options.errorRate / 100.0));
    TSize seedsPerRead = 6;//errorsPerRead + 1;
    TSize seedsLength = readLength / seedsPerRead;

    for (TSize seedId = 0; seedId < seedsPerRead; ++seedId)
        delegate(TPos(readSeqId, seedId * seedsLength), seedsLength);
}

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename TSpec>
inline void selectSeeds(Mapper<TExecSpace, TConfig> & /* mapper */, StringSet<TReadSeqs, Segment<TSpec> > & seeds, TReadSeqs & readSeqs)
{
    typedef StringSet<TReadSeqs, Segment<TSpec> >       TSeedsSet;
    typedef typename Value<TReadSeqs>::Type             TReadSeq;
    typedef typename Id<TReadSeqs>::Type                TId;
    typedef typename Size<TReadSeq>::Type               TSize;
    typedef SeedsCounter<TSize>                         TCounter;
    typedef SeedsManager<TSeedsSet, String<TSize> >     TManager;

    TId readsCount = length(readSeqs);

    TCounter counter(readsCount);
    for (TId readSeqId = 0; readSeqId < readsCount; ++readSeqId)
        _selectSeeds(readSeqs, readSeqId, counter);

    TManager manager(seeds, counter.seedsPerRead);
    for (TId readSeqId = 0; readSeqId < readsCount; ++readSeqId)
        _selectSeeds(readSeqs, readSeqId, manager);
}

// ----------------------------------------------------------------------------
// Function findSeeds()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename THitsString, typename TSeeder, typename TSeeds>
inline void findSeeds(Mapper<TExecSpace, TConfig> & /* mapper */, THitsString & hits, TSeeder & seeder, TSeeds const & seeds)
{
    typedef Mapper<TExecSpace, TConfig>                 TMapper;
//    typedef typename TMapper::THitsString               THitsString;
    typedef typename TMapper::TSeedsExt                 TSeedsPattern;

    HitsManager<THitsString> manager(hits);

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif

    // Instantiate a pattern object.
    TSeedsPattern pattern(seeds);

    // Initialize the delegate.
    init(manager, pattern);

    // Find hits.
    find(seeder, pattern, manager);

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif
}

// ----------------------------------------------------------------------------
// Function filterHits()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename THits, typename TSeeds>
inline void filterHits(Mapper<TExecSpace, TConfig> & mapper, TReadSeqs & readSeqs, THits & hits, TSeeds & seeds)
{
    typedef Mapper<TExecSpace, TConfig>                 TMapper;
    typedef typename TMapper::TSeedsSet                 TSeedsSet;
    typedef typename Id<TSeedsSet>::Type                TSeedId;
    typedef Pair<TSeedId>                               TSeedIds;
    typedef typename TMapper::THit                      THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef Pair<THitId>                                THitIds;
    typedef typename Size<THit>::Type                   THitSize;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    TReadId readsCount = getReadsCount(readSeqs);

    for (TReadId readSeqId = 0; readSeqId < readsCount; ++readSeqId)
    {
        // Get mate id.
        TReadId mateSeqId = getMateSeqId(readSeqs, readSeqId);

        // Get seed ids.
        TSeedIds readSeedIds = getSeedIds(seeds, readSeqId);
        TSeedIds mateSeedIds = getSeedIds(seeds, mateSeqId);

        // Get hit ids.
        THitIds readHitIds = getHitIds(hits, readSeedIds);
        THitIds mateHitIds = getHitIds(hits, mateSeedIds);

        // Count the hits of each read.
        THitSize readHits = countHits<THitSize>(hits, readHitIds);
        THitSize mateHits = countHits<THitSize>(hits, mateHitIds);

        // Choose the easiest read as the anchor.
        THitSize anchorHits = std::min(readHits, mateHits);

        // Clear the hits of the hardest.
        THitIds otherHitIds = (anchorHits == readHits) ? mateHitIds : readHitIds;
        clearHits(hits, otherHitIds);

        // Clear also the hits of the anchor and skip the pair.
        if (anchorHits > mapper.options.hitsThreshold)
        {
            THitIds anchorHitIds = (anchorHits == readHits) ? readHitIds : mateHitIds;
            clearHits(hits, anchorHitIds);
        }
    }
}

// ----------------------------------------------------------------------------
// Function extendHits()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TConfig, typename TReadSeqs, typename THits, typename TSeeds>
inline void extendHits(Mapper<TExecSpace, TConfig> & mapper, TReadSeqs & readSeqs, THits & hits, TSeeds & seeds)
{
    typedef Mapper<TExecSpace, TConfig>                 TMapper;

    typedef typename TMapper::TContigs                  TContigs;
    typedef typename Size<TContigs>::Type               TContigId;
    typedef typename TMapper::TContigsPos               TContigsPos;

    typedef typename TMapper::TReadSeq                  TReadSeq;
    typedef typename Id<TReadSeqs>::Type                TReadId;
    typedef Pair<typename Position<TReadSeq>::Type>     TReadPos;
    typedef typename Size<TReadSeq>::Type               TReadSeqSize;

    typedef typename TMapper::TSeedsSet                 TSeedsSet;
    typedef typename Id<TSeedsSet>::Type                TSeedId;

    typedef typename TMapper::THit                      THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef typename Position<THit>::Type               THitRange;
    typedef unsigned char                               THitErrors;

    typedef typename TMapper::TMatches                  TMatches;
    typedef typename Value<TMatches>::Type              TMatch;

    typedef typename TMapper::TIndex                    TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type       TSA;
    typedef typename Size<TSA>::Type                    TSAPos;
    typedef typename Value<TSA>::Type                   TSAValue;

    typedef MatchesManager<TMatches>                    TManager;

    TSA & sa = indexSA(mapper.index);

    TManager anchorsManager(mapper.anchors);

    THitId hitsCount = length(hits);

    for (THitId hitId = 0; hitId < hitsCount; ++hitId)
    {
        // Extract hit info.
        TSeedId seedId = getSeedId(hits, hitId);
        THitRange hitRange = getRange(hits, hitId);
        THitErrors hitErrors = getErrors(hits, hitId);

        // Get read.
        TReadId readSeqId = getReadSeqId(seeds, seedId);
        TReadSeq readSeq = readSeqs[readSeqId];

        // Fill readSeqId.
        anchorsManager.prototype.readId = readSeqId;

        // Get position in read.
        TReadPos readPos = getPosInRead(seeds, seedId);
        TReadSeqSize seedLength = getValueI2(readPos) - getValueI1(readPos);

        for (TSAPos saPos = getValueI1(hitRange); saPos < getValueI2(hitRange); ++saPos)
        {
            // Invert SA value.
            TSAValue saValue = sa[saPos];
            setSeqOffset(saValue, suffixLength(saValue, contigs(mapper.genome)) - seedLength);

            // Compute position in contig.
            TContigsPos contigBegin = saValue;
            TContigsPos contigEnd = posAdd(contigBegin, seedLength);

            extend(mapper.extender,
                   readSeq,
                   contigBegin, contigEnd,
                   readPos.i1, readPos.i2,
                   hitErrors, mapper.options.errorRate,
                   anchorsManager);
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
    typedef MatchesManager<TMatches>                    TManager;

    typedef typename Size<TReadSeqs>::Type              TReadId;
    typedef typename Value<TReadSeqs>::Type             TReadSeq;

    TManager matesManager(mates);

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

        // Fill readId.
        matesManager.prototype.readId = mateId;

        // TODO(esiragusa): convert this to errorsPerRead.
        verify(mapper.verifier, mateSeq, contigBegin, contigEnd, mapper.options.errorRate, matesManager);
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
    clear(mapper.hits);

    selectSeeds(mapper, mapper.seeds[0], readSeqs);
    std::cout << "Seeds count:\t\t\t" << length(mapper.seeds[0]) << std::endl;
    findSeeds(mapper, mapper.hits[0], mapper.seederExt, mapper.seeds[0]);
    stop(mapper.timer);
    std::cout << "Seeding time:\t\t\t" << mapper.timer << std::endl;

//#ifdef _OPENMP
//    sortHits(mapper.hits[0]);
//#endif
    filterHits(mapper, readSeqs, mapper.hits[0], mapper.seeds[0]);
    std::cout << "Hits count:\t\t\t" << countHits<unsigned long>(mapper.hits[0]) << std::endl;

    start(mapper.timer);
    clear(mapper.anchors);
    reserve(mapper.anchors, countHits<unsigned long>(mapper.hits[0]) / 5);
    extendHits(mapper, readSeqs, mapper.hits[0], mapper.seeds[0]);
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
