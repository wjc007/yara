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
    typedef std::string             TString;
    typedef std::vector<TString>    TList;

    enum MappingMode
    {
        ALL, ALL_BEST, ANY_BEST
    };

    CharString          genomeFile;
    CharString          genomeIndexFile;
    Pair<CharString>    readsFile;

    TList               mappingModeList;
    MappingMode         mappingMode;

    unsigned            errorRate;
    bool                singleEnd;
    bool                anchorOne;
    unsigned            libraryLength;
    unsigned            libraryError;

    unsigned            mappingBlock;
    bool                noCuda;
    unsigned            threadsCount;
    unsigned            hitsThreshold;
    bool                verbose;

    Options() :
        mappingMode(ALL),
        errorRate(5),
        singleEnd(true),
        anchorOne(false),
        libraryLength(220),
        libraryError(50),
        mappingBlock(200000),
        noCuda(false),
        threadsCount(1),
        hitsThreshold(300),
        verbose(true)
    {
        mappingModeList.push_back("all");
        mappingModeList.push_back("all-best");
        mappingModeList.push_back("any-best");
    }
};

// ----------------------------------------------------------------------------
// Mapper Configuration
// ----------------------------------------------------------------------------

template <typename TExecSpace_  = ExecHost,
          typename TSequencing_ = SingleEnd,
          typename TStrategy_   = AnyBest,
          typename TAnchoring_  = Nothing,
          unsigned BUCKETS_     = 3>
struct ReadMapperConfig : public CUDAStoreConfig
{
    typedef TExecSpace_     TExecSpace;
    typedef TSequencing_    TSequencing;
    typedef TStrategy_      TStrategy;
    typedef TAnchoring_     TAnchoring;

    static const unsigned BUCKETS = BUCKETS_;
};

// ----------------------------------------------------------------------------
// Class Mapper
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig = void>
struct Mapper
{
    typedef Genome<void, TConfig>                                   TGenome;
//    typedef GenomeLoader<void, TConfig>                             TGenomeLoader;
    typedef typename Contigs<TGenome>::Type                         TContigs;
    typedef typename Value<TContigs>::Type                          TContig;
    typedef typename StringSetPosition<TContigs>::Type              TContigsPos;

    typedef typename TConfig::TExecSpace                            TExecSpace;
    typedef Index<TFMContigs, TGenomeIndexSpec>                     THostIndex;
    typedef typename Space<THostIndex, TExecSpace>::Type            TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type                   TSA;

    typedef typename TConfig::TSequencing                           TSequencing;
    typedef FragmentStore<void, TConfig>                            TStore;
    typedef ReadsConfig<False, False, True, True, TConfig>          TReadsConfig;
    typedef Reads<TSequencing, TReadsConfig>                        TReads;
    typedef ReadsLoader<TSequencing, TReadsConfig>                  TReadsLoader;
    typedef typename TStore::TReadSeqStore                          THostReadSeqs;
    typedef typename Space<THostReadSeqs, TExecSpace>::Type         TReadSeqs;
    typedef typename Value<TReadSeqs>::Type                         TReadSeq;

    typedef ReadContext<TSpec, TConfig>                             TReadContext;
    typedef String<TReadContext>                                    TReadsContext;

    typedef StringSet<TReadSeqs, Segment<TReadSeqs> >               TSeedsSet;
    typedef Tuple<TSeedsSet, TConfig::BUCKETS>                      TSeeds;

    typedef typename Size<TIndex>::Type                             TIndexSize;
    typedef Hit<TIndexSize, HammingDistance>                        THit;
    typedef String<THit>                                            THitsString;
    typedef Tuple<THitsString, TConfig::BUCKETS>                    THits;

    typedef Match<void>                                             TMatch;
    typedef String<TMatch>                                          TMatches;

    typedef Multiple<FinderSTree>                                   TAlgorithmExt;
    typedef Multiple<Backtracking<HammingDistance> >                TAlgorithmApx;
    typedef Pattern<TSeedsSet, TAlgorithmExt>                       TPatternExt;
    typedef Pattern<TSeedsSet, TAlgorithmApx>                       TPatternApx;
    typedef Finder2<TIndex, TPatternExt, TAlgorithmExt>             TFinderExt;
    typedef Finder2<TIndex, TPatternApx, TAlgorithmApx>             TFinderApx;

//    typedef WriterConfig<Options, TReadSeqs>                        TWriterConfig;
//    typedef Writer<TSpec, TWriterConfig>                            TWriter;

    struct Traits
    {
        typedef TContigs            TContigs;
        typedef TContigsPos         TContigsPos;
        typedef TReadSeqs           TReadSeqs;
        typedef TReadSeq            TReadSeq;
        typedef TReadsContext       TReadsContext;
        typedef TMatches            TMatches;
        typedef TMatch              TMatch;
        typedef TSeedsSet           TSeeds;
        typedef THitsString         THits;
        typedef TSA                 TSA;
    };

    Timer<double>       timer;
    Options const &     options;

    TGenome             genome;
//    TGenomeLoader       genomeLoader;
    TIndex              index;
    TStore              store;
    TReads              reads;
    TReadsLoader        readsLoader;

    TReadsContext       ctx;
    TSeeds              seeds;
    THits               hits;
    TMatches            anchors;
    TMatches            mates;

    TFinderExt          finderExt;
    TFinderApx          finderApx;
//    TWriter             writer;

    Mapper(Options const & options) :
        options(options),
        genome(),
//        genomeLoader(genome),
        index(),
        store(),
        reads(store),
        readsLoader(reads),
        ctx(),
        seeds(),
        hits(),
        anchors(),
        mates(),
        finderExt(index),
        finderApx(index)
//        writer(options, genome)
    {};
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getReadErrors()
// ----------------------------------------------------------------------------
// Returns the absolute number of errors for a given read sequence.

template <typename TReadSeqSize>
inline TReadSeqSize getReadErrors(Options const & options, TReadSeqSize readSeqLength)
{
    return std::ceil(readSeqLength * (options.errorRate / 100.0));
}

// ----------------------------------------------------------------------------
// Function configureThreads()
// ----------------------------------------------------------------------------
// Sets the number of threads that OpenMP can spawn.

template <typename TSpec, typename TConfig>
inline void configureThreads(Mapper<TSpec, TConfig> & mapper)
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

template <typename TSpec, typename TConfig>
inline void loadGenome(Mapper<TSpec, TConfig> & mapper)
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

template <typename TSpec, typename TConfig>
inline void loadGenomeIndex(Mapper<TSpec, TConfig> & mapper)
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
// Function openReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void openReads(Mapper<TSpec, TConfig> & mapper)
{
    _openReadsImpl(mapper, typename TConfig::TSequencing());
}

template <typename TSpec, typename TConfig, typename TSequencing>
inline void _openReadsImpl(Mapper<TSpec, TConfig> & mapper, TSequencing const & /*tag */)
{
    open(mapper.readsLoader, mapper.options.readsFile);
}

template <typename TSpec, typename TConfig>
inline void _openReadsImpl(Mapper<TSpec, TConfig> & mapper, SingleEnd const & /* tag */)
{
    open(mapper.readsLoader, mapper.options.readsFile.i1);
}

// ----------------------------------------------------------------------------
// Function loadReads()
// ----------------------------------------------------------------------------
// Loads one block of reads.

template <typename TSpec, typename TConfig>
inline void loadReads(Mapper<TSpec, TConfig> & mapper)
{
    std::cout << "Loading reads:\t\t\t" << std::flush;
    start(mapper.timer);
    clear(mapper.reads);
    load(mapper.readsLoader, mapper.options.mappingBlock);
    stop(mapper.timer);
    std::cout << mapper.timer << std::endl;

    std::cout << "Reads count:\t\t\t" << mapper.reads.readsCount << std::endl;
}

// ----------------------------------------------------------------------------
// Function clearReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void clearReads(Mapper<TSpec, TConfig> & mapper)
{
    clear(mapper.reads);
}

// ----------------------------------------------------------------------------
// Function initSeeds()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void initSeeds(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs)
{
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
        setHost(mapper.seeds[bucketId], readSeqs);
}

// ----------------------------------------------------------------------------
// Function initReadsContext()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void initReadsContext(Mapper<TSpec, TConfig> & mapper, TReadSeqs const & readSeqs)
{
    clear(mapper.ctx);
    resize(mapper.ctx, getReadSeqsCount(readSeqs));
}

// ----------------------------------------------------------------------------
// Function seedReads()
// ----------------------------------------------------------------------------
// Samples seeds for all reads.

template <typename TSpec, typename TConfig, typename TReadSeqs, typename TErrors>
inline void seedReads(Mapper<TSpec, TConfig> & mapper, TReadSeqs const & readSeqs, Pair<TErrors> seedErrors)
{
    typedef Mapper<TSpec, TConfig>                      TMapper;
    typedef typename Value<TReadSeqs>::Type             TReadSeq;
    typedef typename Id<TReadSeqs>::Type                TId;
    typedef typename Size<TReadSeq>::Type               TSize;
    typedef typename TMapper::TSeeds                    TSeeds;
    typedef typename Value<TSeeds>::Type                TSeedsSet;
    typedef SeedsCounter<TSize>                         TCounter;
    typedef SeedsManager<TSeedsSet, String<TSize> >     TManager;
    typedef Tuple<TCounter, TConfig::BUCKETS>           TCounters;
    typedef Tuple<TManager, TConfig::BUCKETS>           TManagers;

    TId readsCount = getReadSeqsCount(readSeqs);

    // Initialize counters.
    TCounters counters;
    for (TErrors errors = getValueI1(seedErrors); errors <= getValueI2(seedErrors); ++errors)
        resize(counters[errors], readsCount);

    // Count seeds.
    // Counters(ctx)
//    iterate(readSeqs, counters, Rooted(), Parallel());

    for (TId readSeqId = 0; readSeqId < readsCount; ++readSeqId)
    {
        if (getStatus(mapper.ctx, readSeqId) == STATUS_UNSEEDED)
        {
            unsigned char errors = getSeedErrors(mapper.ctx, readSeqId);
            SEQAN_ASSERT_GEQ(errors, getValueI1(seedErrors));
            SEQAN_ASSERT_LEQ(errors, getValueI2(seedErrors));
            _getSeedsPerRead(mapper, readSeqs, readSeqId, counters[errors]);
        }
    }

    // Initialize managers.
    TManagers managers;
    for (TErrors errors = getValueI1(seedErrors); errors <= getValueI2(seedErrors); ++errors)
        init(managers[errors], mapper.seeds[errors], counters[errors].seedsPerRead);

    // Select seeds.
//    iterate(readSeqs, managers, Rooted(), Parallel());

    for (TId readSeqId = 0; readSeqId < readsCount; ++readSeqId)
    {
        if (getStatus(mapper.ctx, readSeqId) == STATUS_UNSEEDED)
        {
            unsigned char errors = getSeedErrors(mapper.ctx, readSeqId);
            SEQAN_ASSERT_GEQ(errors, getValueI1(seedErrors));
            SEQAN_ASSERT_LEQ(errors, getValueI2(seedErrors));
            _getSeedsPerRead(mapper, readSeqs, readSeqId, managers[errors]);
            setStatus(mapper.ctx, readSeqId, STATUS_SEEDED);
        }
    }
}

// ----------------------------------------------------------------------------
// Function _getSeedsPerRead()
// ----------------------------------------------------------------------------
// Enumerates the seeds for a given read sequence.

template <typename TSpec, typename TConfig, typename TReadSeqs, typename TReadSeqId, typename TDelegate>
inline void _getSeedsPerRead(Mapper<TSpec, TConfig> & mapper, TReadSeqs const & readSeqs, TReadSeqId readSeqId, TDelegate & delegate)
{
    typedef typename StringSetPosition<TReadSeqs>::Type     TPos;
    typedef typename Value<TReadSeqs>::Type                 TReadSeq;
    typedef typename Size<TReadSeq>::Type                   TSize;

    TSize readLength = length(readSeqs[readSeqId]);
    TSize readErrors = getReadErrors(mapper.options, readLength);
    TSize seedErrors = getSeedErrors(mapper.ctx, readSeqId);
    TSize seedsCount = std::ceil((readErrors + 1) / (seedErrors + 1.0));
    TSize seedsLength = readLength / seedsCount;

    for (TSize seedId = 0; seedId < seedsCount; ++seedId)
        delegate(TPos(readSeqId, seedId * seedsLength), seedsLength);
}

// ----------------------------------------------------------------------------
// Function findSeeds()
// ----------------------------------------------------------------------------

template <unsigned ERRORS, typename TSpec, typename TConfig, typename TBucketId>
inline void findSeeds(Mapper<TSpec, TConfig> & mapper, TBucketId bucketId)
{
    typedef Mapper<TSpec, TConfig>          TMapper;
    typedef typename TMapper::TPatternExt   TPatternExt;
    typedef typename TMapper::TPatternApx   TPatternApx;

    start(mapper.timer);
    std::cout << "Seeds count:\t\t\t" << length(mapper.seeds[bucketId]) << std::endl;

    if (ERRORS > 0)
    {
        setScoreThreshold(mapper.finderApx, ERRORS);
        _findSeedsImpl(mapper, mapper.hits[bucketId], mapper.seeds[bucketId], mapper.finderApx, TPatternApx());
    }
    else
    {
        _findSeedsImpl(mapper, mapper.hits[bucketId], mapper.seeds[bucketId], mapper.finderExt, TPatternExt());
    }

    stop(mapper.timer);
    std::cout << "Seeding time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Hits count:\t\t\t" << countHits<unsigned long>(mapper.hits[bucketId]) << std::endl;
//    writeHits(mapper, readSeqs, mapper.hits[bucketId], mapper.seeds[bucketId], "hits.csv");
}

// ----------------------------------------------------------------------------
// Function _findSeedsImpl()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename THitsString, typename TSeedsSet, typename TFinder, typename TPattern>
inline void _findSeedsImpl(Mapper<TSpec, TConfig> & /* mapper */, THitsString & hits, TSeedsSet & seeds, TFinder & finder, TPattern)
{
    HitsManager<THitsString> manager(hits);

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif

    TPattern pattern(seeds);

    // Initialize the delegate.
    init(manager, pattern);

    // Find hits.
    find(finder, pattern, manager);

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif

#ifdef _OPENMP
    // Sort the hits by seedId.
    sortHits(hits);
#endif
}

// ----------------------------------------------------------------------------
// Function classifyReads()
// ----------------------------------------------------------------------------
// Classifies the reads by hardness.

template <typename TSpec, typename TConfig>
inline void classifyReads(Mapper<TSpec, TConfig> & mapper)
{
    typedef Mapper<TSpec, TConfig>                                  TMapper;
    typedef typename TMapper::TReadsContext                         TReadsContext;
    typedef typename TMapper::THitsString                           THits;
    typedef typename TMapper::TSeedsSet                             TSeeds;
    typedef ReadsClassifier<TReadsContext, THits, TSeeds, TConfig>  TClassifier;

    TClassifier classifier(mapper.ctx, mapper.hits[0], mapper.seeds[0], mapper.options);

    std::cout << "Hits count:\t\t\t" << countHits<unsigned long>(mapper.hits[0]) << std::endl;
}

// ----------------------------------------------------------------------------
// Function clearHits()
// ----------------------------------------------------------------------------
// Clears the hits in all buckets.

template <typename TSpec, typename TConfig>
inline void clearHits(Mapper<TSpec, TConfig> & mapper)
{
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
        clear(mapper.hits[bucketId]);
}

// ----------------------------------------------------------------------------
// Function countHits()
// ----------------------------------------------------------------------------
// Counts the hits in all buckets.

template <typename TSpec, typename TConfig>
inline unsigned long countHits(Mapper<TSpec, TConfig> const & mapper)
{
    unsigned long hitsCount = 0;

    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
        hitsCount += countHits<unsigned long>(mapper.hits[bucketId]);

    return hitsCount;
}

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class HitsExtender
// ----------------------------------------------------------------------------
// One instance per thread.

template <typename TSpec, typename TConfig>
struct HitsExtender
{
    typedef typename TConfig::TContigs          TContigs;
    typedef typename TConfig::TContigsPos       TContigsPos;
    typedef typename TConfig::TReadSeqs         TReadSeqs;
    typedef typename TConfig::TReadSeq          TReadSeq;
    typedef typename TConfig::TReadsContext     TReadsContext;
    typedef typename TConfig::TMatches          TMatches;
    typedef typename TConfig::TMatch            TMatch;
    typedef typename TConfig::TSeeds            TSeeds;
    typedef typename TConfig::THits             THits;
    typedef typename TConfig::TSA               TSA;

    typedef AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_> TMyersSpec;
    typedef Myers<TMyersSpec, True, void>               TAlgorithm;
    typedef Extender<TContigs, TReadSeq, TAlgorithm>    TExtender;

    // Thread-private data.
    TExtender           extender;
    TMatch              prototype;

    // Shared-memory read-write data.
    TReadsContext &     ctx;
    TMatches &          matches;

    // Shared-memory read-only data.
    TContigs const &    contigs;
    TReadSeqs const &   readSeqs;
    TSeeds const &      seeds;
    THits const &       hits;
    TSA const &         sa;
    Options const &     options;

    HitsExtender(TReadsContext & ctx,
                 TMatches & matches,
                 TContigs const & contigs,
                 TSeeds const & seeds,
                 THits const & hits,
                 TSA const & sa,
                 Options const & options) :
        extender(contigs),
        prototype(),
        ctx(ctx),
        matches(matches),
        contigs(contigs),
        readSeqs(host(seeds)),
        seeds(seeds),
        hits(hits),
        sa(sa),
        options(options)
    {
        // Iterate over all hits.
        iterate(hits, *this, Rooted(), Serial());
    }

    template <typename THitsIterator>
    void operator() (THitsIterator const & hitsIt)
    {
        _extendHitImpl(*this, hitsIt);
    }

    void operator() (TExtender const & /* extender */)
    {
//        _addMatchImpl(*this, extender);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _extendHitImpl()
// ----------------------------------------------------------------------------
// Extends one hit.

template <typename TSpec, typename TConfig, typename THitsIterator>
inline void _extendHitImpl(HitsExtender<TSpec, TConfig> & me, THitsIterator const & hitsIt)
{
    typedef HitsExtender<TSpec, TConfig>                THitsExtender;

    typedef typename THitsExtender::TContigs            TContigs;
    typedef typename Size<TContigs>::Type               TContigId;
    typedef typename THitsExtender::TContigsPos         TContigsPos;

    typedef typename THitsExtender::TReadSeq            TReadSeq;
    typedef typename Id<TReadSeqs>::Type                TReadId;
    typedef Pair<typename Position<TReadSeq>::Type>     TReadPos;
    typedef typename Size<TReadSeq>::Type               TReadSeqSize;
    typedef typename Size<TReadSeq>::Type               TErrors;

    typedef typename THitsExtender::TSeeds              TSeedsSet;
    typedef typename Id<TSeedsSet>::Type                TSeedId;

    typedef typename THitsExtender::THits               THits;
    typedef typename Value<THits>::Type                 THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef typename Position<THit>::Type               THitRange;
    typedef unsigned char                               THitErrors;

    typedef typename THitsExtender::TMatches            TMatches;
    typedef typename Value<TMatches>::Type              TMatch;

    typedef typename THitsExtender::TSA                 TSA;
    typedef typename Size<TSA>::Type                    TSAPos;
    typedef typename Value<TSA>::Type                   TSAValue;

    typedef typename THitsExtender::TReadsContext       TReadsContext;

    // Get hit id.
    THitId hitId = position(hitsIt);

    // Extract hit info.
    TSeedId seedId = getSeedId(me.hits, hitId);
    THitRange hitRange = getRange(me.hits, hitId);
    THitErrors hitErrors = getErrors(me.hits, hitId);

    // Get read.
    TReadId readSeqId = getReadSeqId(me.seeds, seedId);
    TReadSeq readSeq;// = me.readSeqs[readSeqId];

    // Skip unseeded and mapped reads.
    if (getStatus(me.ctx, readSeqId) != STATUS_SEEDED) return;

    // Fill readSeqId.
    me.prototype.readId = readSeqId;

    // Get position in read.
    TReadPos readPos = getPosInRead(me.seeds, seedId);
    TReadSeqSize seedLength = getValueI2(readPos) - getValueI1(readPos);

    for (TSAPos saPos = getValueI1(hitRange); saPos < getValueI2(hitRange); ++saPos)
    {
        // Invert SA value.
        TSAValue saValue = me.sa[saPos];
        setSeqOffset(saValue, suffixLength(saValue, me.contigs) - seedLength);

        // Compute position in contig.
        TContigsPos contigBegin = saValue;
        TContigsPos contigEnd = posAdd(contigBegin, seedLength);

        // Get absolute number of errors.
        TErrors maxErrors = getReadErrors(me.options, length(readSeq));

        extend(me.extender,
               readSeq,
               contigBegin, contigEnd,
               readPos.i1, readPos.i2,
               hitErrors, maxErrors,
               me);

        // Stop when the read has been mapped.
        if (getStatus(me.ctx, readSeqId) == STATUS_MAPPED) break;
    }

    // Full stratum analyzed.
    // TODO(esiragusa): generalize stratum for approximate seeds, where one seed can have multiple hits.
    // TODO(esiragusa): extend any fwd + any rev seed consecutively, then increase the stratum of both fwd & rev read.
    incStratum(me.ctx, readSeqId);
}

//template <typename TSpec, typename TConfig, typename TExtender>
//inline void _addMatchImpl(HitsExtender<TSpec, TConfig> & me, TExtender const & extender)
//{
//    me.prototype.contigId = getValueI1(matchBegin);
//    me.prototype.contigBegin = getValueI2(matchBegin);
//    me.prototype.contigEnd = getValueI2(matchEnd);
//    me.prototype.errors = errors;
//    appendValue(me.matches, me.prototype); // typename TConfig::TParallel()
//}

// ----------------------------------------------------------------------------
// Function extendHits()
// ----------------------------------------------------------------------------
// Extends the hits in all buckets.

template <typename TSpec, typename TConfig>
inline void extendHits(Mapper<TSpec, TConfig> & mapper)
{
    typedef Mapper<TSpec, TConfig>          TMapper;
    typedef typename TMapper::Traits        TTraits;
    typedef HitsExtender<TSpec, TTraits>    THitsExtender;

    // TODO(esiragusa): guess the number of matches.
    clear(mapper.anchors);
    reserve(mapper.anchors, countHits(mapper) / 5);

    start(mapper.timer);
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
    {
        THitsExtender extender(mapper.ctx, mapper.anchors, contigs(mapper.genome),
                               mapper.seeds[bucketId], mapper.hits[bucketId],
                               indexSA(mapper.index), mapper.options);
    }
    stop(mapper.timer);

    std::cout << "Extension time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << length(mapper.anchors) << std::endl;
}

// ----------------------------------------------------------------------------
// Function _writeHitsImpl()
// ----------------------------------------------------------------------------
// Debug stuff.

template <typename TSpec, typename TConfig, typename TReadSeqs, typename THits, typename TSeedsSet, typename TFilename>
inline void _writeHitsImpl(Mapper<TSpec, TConfig> const & mapper, TReadSeqs const & readSeqs, THits const & hits, TSeedsSet const & seeds, TFilename const & filename)
{
    typedef Mapper<TSpec, TConfig>                      TMapper;
    typedef typename Id<TSeedsSet>::Type                TSeedId;
    typedef Pair<TSeedId>                               TSeedIds;
    typedef typename TMapper::THit                      THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef typename Id<THit>::Type                     THitId;
    typedef Pair<THitId>                                THitIds;
    typedef typename Position<THit>::Type               THitRange;
    typedef typename Size<THit>::Type                   THitSize;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    std::ofstream file;
    file.open(filename);

    TReadId readSeqsCount = getReadSeqsCount(readSeqs);

    for (TReadId readSeqId = 0; readSeqId < readSeqsCount; ++readSeqId)
    {
        TSeedIds readSeedIds = getSeedIds(seeds, readSeqId);

        if (getStatus(mapper.ctx, readSeqId) != STATUS_SEEDED) continue;
        if (getValueI2(readSeedIds) <= getValueI1(readSeedIds)) continue;

        for (TSeedId seedId = getValueI1(readSeedIds); seedId < getValueI2(readSeedIds); ++seedId)
        {
            THitIds seedHitIds = getHitIds(hits, seedId);
            file << countHits<THitSize>(hits, seedHitIds) << "\t";
        }
        file << "\n";
    }

    file.close();
}

// ----------------------------------------------------------------------------
// Function verifyMates()
// ----------------------------------------------------------------------------
// Verifies all mates in within the insert window of their anchors.

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void verifyMates(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs)
{
    start(mapper.timer);
    _verifyMatesImpl(mapper, readSeqs, typename TConfig::TAnchoring());
    stop(mapper.timer);
    std::cout << "Verification time:\t\t" << mapper.timer << std::endl;
    std::cout << "Mates count:\t\t\t" << length(mapper.mates) << std::endl;
    std::cout << "Mapped pairs:\t\t\t" << countMatches(readSeqs, mapper.mates, typename TConfig::TSequencing()) << std::endl;
}

template <typename TSpec, typename TConfig, typename TReadSeqs, typename TAnchoring>
inline void _verifyMatesImpl(Mapper<TSpec, TConfig> & /* mapper */, TReadSeqs & /* readSeqs */, TAnchoring)
{}

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _verifyMatesImpl(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs, AnchorOne)
{
    typedef Mapper<TSpec, TConfig>                      TMapper;

    typedef typename TMapper::TContigs                  TContigs;
    typedef typename TMapper::TContigsPos               TContigsPos;

    typedef typename Size<TReadSeqs>::Type              TReadId;
    typedef typename Value<TReadSeqs>::Type             TReadSeq;
    typedef typename Size<TReadSeq>::Type               TErrors;

    typedef typename TMapper::TMatches                  TMatches;
    typedef typename Size<TMatches>::Type               TMatchId;
    typedef typename Value<TMatches>::Type              TMatch;

    typedef typename TMapper::TReadsContext             TReadsContext;
    typedef MatchesManager<TReadSeqs, TReadsContext, TMatches>  TManager;

    typedef Myers<>                                     TAlgorithm;
//    typedef Filter<MultipleShiftAnd>                    TAlgorithm;
    typedef Verifier<TContigs, TReadSeq, TAlgorithm>    TVerifier;

    TVerifier verifier(contigs(mapper.genome));

    TManager matesManager(readSeqs, mapper.ctx, mapper.mates);
    clear(mapper.mates);
    reserve(mapper.mates, length(mapper.anchors), Exact());

//    TVerifier verifier(contigs(mapper.genome), ctx, mates);
//    iterate(anchors, extender, Rooted(), Parallel());

    TMatchId matchesCount = length(mapper.anchors);
    for (TMatchId matchId = 0; matchId < matchesCount; ++matchId)
    {
        TMatch const & match = mapper.anchors[matchId];
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

        // Get absolute number of errors.
        TErrors maxErrors = getReadErrors(mapper.options, length(mateSeq));

        verify(verifier, mateSeq, contigBegin, contigEnd, maxErrors, matesManager);
    }
}

// ----------------------------------------------------------------------------
// Function _getMateContigPos()
// ----------------------------------------------------------------------------
// Computes the insert window.

template <typename TSpec, typename TConfig, typename TContigPos, typename TMatch>
inline void _getMateContigPos(Mapper<TSpec, TConfig> const & mapper,
                              TContigPos & contigBegin,
                              TContigPos & contigEnd,
                              TMatch const & anchor,
                              RightMate)
{
    typedef Mapper<TSpec, TConfig>                  TMapper;
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

template <typename TSpec, typename TConfig, typename TContigPos, typename TMatch>
inline void _getMateContigPos(Mapper<TSpec, TConfig> const & mapper,
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
// Function removeDuplicates()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void removeDuplicates(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs)
{
    start(mapper.timer);
    removeDuplicateMatches(mapper.anchors);
    stop(mapper.timer);
    std::cout << "Compaction time:\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << length(mapper.anchors) << std::endl;
    std::cout << "Anchored pairs:\t\t\t" << countMatches(readSeqs, mapper.anchors, typename TConfig::TSequencing()) << std::endl;
}

// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void mapReads(Mapper<TSpec, TConfig> & mapper)
{
//SEQAN_OMP_PRAGMA(critical(_mapper_mapReads_filter))
//{
    _mapReadsImpl(mapper, getSeqs(mapper.reads), typename TConfig::TSequencing(), typename TConfig::TStrategy());
//}
}

// ----------------------------------------------------------------------------
// Function _mapReadsImpl(); SingleEnd or PairedEnd, All
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs, typename TSequencing>
inline void _mapReadsImpl(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs, TSequencing, All)
{
    initReadsContext(mapper, readSeqs);
    initSeeds(mapper, readSeqs);
    clearHits(mapper);

    seedReads(mapper, readSeqs, Pair<unsigned>(0, 0));
    findSeeds<0>(mapper, 0);
    classifyReads(mapper);
    seedReads(mapper, readSeqs, Pair<unsigned>(1, 2));
    findSeeds<1>(mapper, 1);
    findSeeds<2>(mapper, 2);
    extendHits(mapper);
    removeDuplicates(mapper, readSeqs);
    verifyMates(mapper, readSeqs);
}

// ----------------------------------------------------------------------------
// Function _mapReadsImpl(); SingleEnd, AnyBest or AllBest
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _mapReadsImpl(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs, SingleEnd, AnyBest)
{
    _mapReadsByStrata(mapper, readSeqs);
}

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _mapReadsImpl(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs, SingleEnd, AllBest)
{
    _mapReadsByStrata(mapper, readSeqs);
}

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _mapReadsByStrata(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs)
{
    typedef Mapper<TSpec, TConfig>          TMapper;
    typedef typename TMapper::THit          THit;

    initReadsContext(mapper, readSeqs);
    initSeeds(mapper, readSeqs);
    clearHits(mapper);

    start(mapper.timer);
    seedReads(mapper, readSeqs, Pair<unsigned>(0, 0));
    findSeeds<0>(mapper, 0);
    classifyReads(mapper);
    seedReads(mapper, readSeqs, Pair<unsigned>(1, 2));
    findSeeds<0>(mapper, 1);
    findSeeds<0>(mapper, 2);

    // Sort hits by range size.
    // TODO(esiragusa): generalize sorting for approximate seeds, where one seed can have multiple hits.
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; ++bucketId)
        std::stable_sort(begin(mapper.hits[bucketId], Standard()), end(mapper.hits[bucketId], Standard()), HitsSorterByCount<THit>());

    extendHits(mapper);

    std::cout << "Mapped reads:\t\t\t" << countMapped(mapper.ctx) << std::endl;

//    clearHits(mapper);
//    seedReads(mapper, readSeqs, Pair<unsigned>(1, 2));
//    findSeeds<1>(mapper, 1);
//    findSeeds<1>(mapper, 2);
////    sortHits(mapper);
//    extendHits(mapper);
//
//    clearHits(mapper);
//    seedReads(mapper, readSeqs, Pair<unsigned>(2, 2));
//    findSeeds<2>(mapper, 2);
////    sortHits(mapper);
//    extendHits(mapper);
//
//    extendHits(mapper);
//    removeDuplicates(mapper, readSeqs);
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void runMapper(Mapper<TSpec, TConfig> & mapper)
{
    configureThreads(mapper);

    printRuler();

    loadGenome(mapper);
    loadGenomeIndex(mapper);

    // Open reads file.
    openReads(mapper);

    // Process reads in blocks.
    while (!atEnd(mapper.readsLoader))
    {
        printRuler();

        loadReads(mapper);
        mapReads(mapper);
        clearReads(mapper);
    }

    // Close reads file.
    close(mapper.readsLoader);
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

/*
template <typename TSpec, typename TConfig>
inline void runMapper(Mapper<TSpec, TConfig> & mapper, Parallel)
{
    typedef Mapper<TSpec, TConfig>                      TMapper;
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
*/

// ----------------------------------------------------------------------------
// Function spawnMapper()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TSequencing, typename TStrategy, typename TAnchoring>
inline void spawnMapper(Options const & options,
                        TExecSpace const & /* tag */,
                        TSequencing const & /* tag */,
                        TStrategy const & /* tag */,
                        TAnchoring const & /* tag */)
{
    typedef ReadMapperConfig<TExecSpace, TSequencing, TStrategy, TAnchoring> TConfig;

    Mapper<void, TConfig> mapper(options);
    runMapper(mapper);
}

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_H_
