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

    enum OutputFormat
    {
        SAM, BAM
    };

    CharString          genomeFile;
    CharString          genomeIndexFile;
    Pair<CharString>    readsFile;
    CharString          outputFile;
    OutputFormat        outputFormat;

    TList               mappingModeList;
    MappingMode         mappingMode;

    TList               outputFormatList;
    TList               outputFormatExtensions;

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
        outputFormat(SAM),
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

        outputFormatList.push_back("sam");
        outputFormatExtensions.push_back("sam");
#ifdef SEQAN_HAS_ZLIB
        outputFormatList.push_back("bam");
        outputFormatExtensions.push_back("bam");
#endif
    }
};

// ----------------------------------------------------------------------------
// Mapper Configuration
// ----------------------------------------------------------------------------

template <typename TExecSpace_      = ExecHost,
          typename TThreading_      = Parallel,
          typename TOutputFormat_   = Sam,
          typename TSequencing_     = SingleEnd,
          typename TStrategy_       = AnyBest,
          typename TAnchoring_      = Nothing,
          unsigned BUCKETS_         = 3>
struct ReadMapperConfig : public ContigsConfig<void>, public ReadsConfig<void>
{
    typedef TExecSpace_     TExecSpace;
    typedef TThreading_     TThreading;
    typedef TOutputFormat_  TOutputFormat;
    typedef TSequencing_    TSequencing;
    typedef TStrategy_      TStrategy;
    typedef TAnchoring_     TAnchoring;

    static const unsigned BUCKETS = BUCKETS_;
};

// ----------------------------------------------------------------------------
// Mapper Traits
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
struct MapperTraits
{
    typedef typename TConfig::TExecSpace                            TExecSpace;
    typedef typename TConfig::TThreading                            TThreading;
    typedef typename TConfig::TOutputFormat                         TOutputFormat;
    typedef typename TConfig::TSequencing                           TSequencing;
    typedef typename TConfig::TStrategy                             TStrategy;
    typedef typename TConfig::TAnchoring                            TAnchoring;

    typedef Contigs<void, TConfig>                                  TContigs;
    typedef typename TContigs::TContigSeqs                          TContigSeqs;
    typedef typename Value<TContigSeqs>::Type                       TContig;
    typedef typename StringSetPosition<TContigSeqs>::Type           TContigsPos;

    typedef Index<TFMContigs, TGenomeIndexSpec>                     THostIndex;
    typedef typename Space<THostIndex, TExecSpace>::Type            TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type                   TSA;

    typedef Reads<TSequencing, TConfig>                             TReads;
    typedef ReadsLoader<TSequencing, TConfig>                       TReadsLoader;
    typedef typename TReads::TReadSeqs                              THostReadSeqs;
    typedef typename Space<THostReadSeqs, TExecSpace>::Type         TReadSeqs;
    typedef typename Value<TReadSeqs>::Type                         TReadSeq;
    typedef typename Size<TReadSeqs>::Type                          TReadSeqSize;
    typedef String<TReadSeqSize>                                    TSeedsCount;

    typedef typename TContigs::TContigNames                         TContigNames;
    typedef typename TContigs::TContigNamesCache                    TContigNamesCache;
    typedef Stream<FileStream<File<> > >                            TOutputStream;
    typedef BamIOContext<TContigNames, TContigNamesCache>           TOutputContext;

    typedef ReadContext<TSpec, TConfig>                             TReadContext;
    typedef String<TReadContext>                                    TReadsContext;

    typedef StringSet<TReadSeqs, Segment<TReadSeqs> >               TSeeds;
    typedef Tuple<TSeeds, TConfig::BUCKETS>                         TSeedsBuckets;

    typedef typename Size<TIndex>::Type                             TIndexSize;
    typedef Hit<TIndexSize, HammingDistance>                        THit;
    typedef String<THit>                                            THits;
    typedef Tuple<THits, TConfig::BUCKETS>                          THitsBuckets;

    typedef Match<void>                                             TMatch;
    typedef String<TMatch>                                          TMatches;
    typedef StringSet<TMatches, Segment<TMatches> >                 TMatchesSet;

    typedef Multiple<FinderSTree>                                   TAlgorithmExt;
    typedef Multiple<Backtracking<HammingDistance> >                TAlgorithmApx;
    typedef Pattern<TSeeds, TAlgorithmExt>                          TPatternExt;
    typedef Pattern<TSeeds, TAlgorithmApx>                          TPatternApx;
    typedef Finder2<TIndex, TPatternExt, TAlgorithmExt>             TFinderExt;
    typedef Finder2<TIndex, TPatternApx, TAlgorithmApx>             TFinderApx;
};

// ----------------------------------------------------------------------------
// Class Mapper
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig = void>
struct Mapper
{
    typedef MapperTraits<TSpec, TConfig>    Traits;

    Timer<double>                       timer;
    Options const &                     options;

    typename Traits::TContigs           contigs;
    typename Traits::TIndex             index;
    typename Traits::TReads             reads;
    typename Traits::TReadsLoader       readsLoader;
    typename Traits::TOutputStream      outputStream;
    typename Traits::TOutputContext     outputCtx;

    typename Traits::TReadsContext      ctx;
    typename Traits::TSeedsBuckets      seeds;
    typename Traits::THitsBuckets       hits;
    typename Traits::TMatches           anchors;
    typename Traits::TMatches           mates;
    typename Traits::TMatchesSet        anchorsSet;
    typename Traits::TMatchesSet        matesSet;

    typename Traits::TFinderExt         finderExt;
    typename Traits::TFinderApx         finderApx;

    Mapper(Options const & options) :
        options(options),
        contigs(),
        index(),
        reads(),
        readsLoader(),
        outputStream(),
        outputCtx(contigs._contigNames, contigs._contigNamesCache),
        ctx(),
        seeds(),
        hits(),
        anchors(),
        mates(),
        anchorsSet(),
        matesSet(),
        finderExt(index),
        finderApx(index)
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

    if (!open(mapper.contigs, toCString(mapper.options.genomeIndexFile)))
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
    load(mapper.reads, mapper.readsLoader, mapper.options.mappingBlock);
    stop(mapper.timer);
    std::cout << mapper.timer << std::endl;
    std::cout << "Reads count:\t\t\t" << getReadsCount(mapper.reads._readSeqs) << std::endl;
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
// Function initOutput()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void initOutput(Mapper<TSpec, TConfig> & mapper)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;

    BamHeader header;

    if (!open(mapper.outputStream, toCString(mapper.options.outputFile), OPEN_RDWR | OPEN_CREATE))
        throw RuntimeError("Error while opening output file.");

    // Fill header.
//    _fillHeader(header, mapper.contigs);

    // Write header to stream.
    write2(mapper.outputStream, header, mapper.outputCtx, typename TTraits::TOutputFormat());
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
// Function collectSeeds()
// ----------------------------------------------------------------------------
// Collects seeds from all reads.

template <unsigned ERRORS, typename TSpec, typename TConfig, typename TReadSeqs>
inline void collectSeeds(Mapper<TSpec, TConfig> & mapper, TReadSeqs const & readSeqs)
{
    typedef MapperTraits<TSpec, TConfig>                TTraits;
    typedef SeedsCollector<Counter, TTraits>            TCounter;
    typedef SeedsCollector<void, TTraits>               TFiller;

    typename TTraits::TSeedsCount seedsCounts;

    TCounter counter(mapper.ctx, mapper.seeds[ERRORS], seedsCounts, readSeqs, mapper.options, ERRORS);
    TFiller filler(mapper.ctx, mapper.seeds[ERRORS], seedsCounts, readSeqs, mapper.options, ERRORS);
}

// ----------------------------------------------------------------------------
// Function findSeeds()
// ----------------------------------------------------------------------------

template <unsigned ERRORS, typename TSpec, typename TConfig, typename TBucketId>
inline void findSeeds(Mapper<TSpec, TConfig> & mapper, TBucketId bucketId)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;
    typedef typename TTraits::TPatternExt           TPatternExt;
    typedef typename TTraits::TPatternApx           TPatternApx;

    start(mapper.timer);
    std::cout << "Seeds count:\t\t\t" << length(mapper.seeds[bucketId]) << std::endl;

    if (ERRORS > 0)
    {
        setScoreThreshold(mapper.finderApx, ERRORS);
        reserve(mapper.hits[bucketId], lengthSum(mapper.seeds[bucketId]) * Power<ERRORS, 2>::VALUE, Exact());
        _findSeedsImpl(mapper, mapper.hits[bucketId], mapper.seeds[bucketId], mapper.finderApx, TPatternApx());
    }
    else
    {
        reserve(mapper.hits[bucketId], length(mapper.seeds[bucketId]), Exact());
        _findSeedsImpl(mapper, mapper.hits[bucketId], mapper.seeds[bucketId], mapper.finderExt, TPatternExt());
    }

    stop(mapper.timer);
    std::cout << "Seeding time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Hits count:\t\t\t" << countHits<unsigned long>(mapper.hits[bucketId], typename TConfig::TThreading()) << std::endl;
//    writeHits(mapper, readSeqs, mapper.hits[bucketId], mapper.seeds[bucketId], "hits.csv");
}

template <typename TSpec, typename TConfig, typename THits, typename TSeeds, typename TFinder, typename TPattern>
inline void _findSeedsImpl(Mapper<TSpec, TConfig> & /* mapper */, THits & hits, TSeeds & seeds, TFinder & finder, TPattern)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;
    typedef FilterDelegate<TSpec, TTraits>          TDelegate;

    TDelegate delegate(hits);
    TPattern pattern(seeds);

    // Find hits.
    find(finder, pattern, delegate);

    // Sort the hits by seedId.
    if (IsSameType<typename TConfig::TThreading, Parallel>::VALUE)
        sortHits(hits, typename TConfig::TThreading());
}

// ----------------------------------------------------------------------------
// Function classifyReads()
// ----------------------------------------------------------------------------
// Classifies the reads by hardness.

template <typename TSpec, typename TConfig>
inline void classifyReads(Mapper<TSpec, TConfig> & mapper)
{
    typedef MapperTraits<TSpec, TConfig>                TTraits;
    typedef ReadsClassifier<TSpec, TTraits>             TClassifier;

    TClassifier classifier(mapper.ctx, mapper.hits[0], mapper.seeds[0], mapper.options);

    std::cout << "Hits count:\t\t\t" << countHits<unsigned long>(mapper.hits[0], typename TConfig::TThreading()) << std::endl;
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
        hitsCount += countHits<unsigned long>(mapper.hits[bucketId], typename TConfig::TThreading());

    return hitsCount;
}

// ----------------------------------------------------------------------------
// Function extendHits()
// ----------------------------------------------------------------------------
// Extends the hits in all buckets.

template <typename TSpec, typename TConfig>
inline void extendHits(Mapper<TSpec, TConfig> & mapper)
{
    typedef MapperTraits<TSpec, TConfig>    TTraits;
    typedef HitsExtender<TSpec, TTraits>    THitsExtender;

    start(mapper.timer);

    // TODO(esiragusa): guess the number of matches.
    clear(mapper.anchors);
    reserve(mapper.anchors, countHits(mapper) / 3);

    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
    {
        THitsExtender extender(mapper.ctx, mapper.anchors, mapper.contigs._contigSeqs,
                               mapper.seeds[bucketId], mapper.hits[bucketId],
                               indexSA(mapper.index), mapper.options);
    }

    stop(mapper.timer);

    std::cout << "Extension time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << length(mapper.anchors) << std::endl;
}

// --------------------------------------------------------------------------
// Function aggregate()
// --------------------------------------------------------------------------
// TODO(esiragusa): Parallelize and move into Segment StringSet.

template <typename THost, typename TSpec, typename TKey>
inline void aggregate(StringSet<THost, Segment<TSpec> > & me, TKey const & key)
{
    typedef typename Iterator<THost, Standard>::Type    THostIter;

    THostIter beginIt = begin(host(me), Standard());
    THostIter endIt = end(host(me), Standard());
    THostIter firstIt = beginIt;
    THostIter lastIt = firstIt;

    clear(me);

    while (firstIt != endIt)
    {
        while (lastIt != endIt && key(value(firstIt)) == key(value(lastIt))) ++lastIt;

        appendInfixWithLength(me, firstIt - beginIt, lastIt - firstIt, Generous());

        firstIt = lastIt;
    }
}

// ----------------------------------------------------------------------------
// Function aggregateAnchors()
// ----------------------------------------------------------------------------
// Aggregate anchors by readId.

template <typename TSpec, typename TConfig>
inline void aggregateAnchors(Mapper<TSpec, TConfig> & mapper)
{
    typedef MapperTraits<TSpec, TConfig>    TTraits;
    typedef typename TTraits::TMatch        TMatch;

    start(mapper.timer);

    // Sort anchors by readId.
//    if (IsSameType<typename TConfig::TThreading, Parallel>::VALUE)
    sort(mapper.anchors, MatchSorter<TMatch, SortReadId>(), typename TConfig::TThreading());

    setHost(mapper.anchorsSet, mapper.anchors);
    aggregate(mapper.anchorsSet, MatchReadId<TMatch>());

    stop(mapper.timer);

    std::cout << "Sorting time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Mapped anchors:\t\t\t" << length(mapper.anchorsSet) << std::endl;
}

// ----------------------------------------------------------------------------
// Function compactAnchors()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void compactAnchors(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs)
{
    typedef MapperTraits<TSpec, TConfig>    TTraits;

    start(mapper.timer);

    removeDuplicates(mapper.anchorsSet, typename TConfig::TThreading());

    stop(mapper.timer);
    std::cout << "Compaction time:\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << lengthSum(mapper.anchorsSet) << std::endl;
    std::cout << "Anchored pairs:\t\t\t" << countMatches(readSeqs, concat(mapper.anchorsSet),
                                                         typename TConfig::TSequencing(),
                                                         typename TConfig::TThreading()) << std::endl;
}

// ----------------------------------------------------------------------------
// Function verifyAnchors()
// ----------------------------------------------------------------------------
// Verifies all mates in within the insert window of their anchors.

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void verifyAnchors(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs)
{
    _verifyAnchorsImpl(mapper, readSeqs, typename TConfig::TAnchoring());
}

template <typename TSpec, typename TConfig, typename TReadSeqs, typename TAnchoring>
inline void _verifyAnchorsImpl(Mapper<TSpec, TConfig> & /* mapper */, TReadSeqs & /* readSeqs */, TAnchoring) {}

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _verifyAnchorsImpl(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs, AnchorOne)
{
    typedef MapperTraits<TSpec, TConfig>    TTraits;
    typedef AnchorsVerifier<TSpec, TTraits> TAnchorsVerifier;

    start(mapper.timer);

    // TODO(esiragusa): guess the number of mates.
    clear(mapper.mates);
    reserve(mapper.mates, lengthSum(mapper.anchorsSet));

    TAnchorsVerifier verifier(mapper.ctx, mapper.mates,
                              mapper.contigs._contigSeqs, readSeqs,
                              concat(mapper.anchorsSet), mapper.options);

    // Sort mates by readId.
//    if (IsSameType<typename TConfig::TThreading, Parallel>::VALUE)
//    sort(mapper.mates, MatchSorter<typename TTraits::TMatch, SortReadId>(), typename TConfig::TThreading());

    stop(mapper.timer);

    std::cout << "Verification time:\t\t" << mapper.timer << std::endl;
    std::cout << "Mates count:\t\t\t" << length(mapper.mates) << std::endl;
    std::cout << "Mapped pairs:\t\t\t" << countMatches(readSeqs, mapper.mates,
                                                       typename TConfig::TSequencing(),
                                                       typename TConfig::TThreading()) << std::endl;
}

// ----------------------------------------------------------------------------
// Function writeMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void writeMatches(Mapper<TSpec, TConfig> & mapper)
{
    typedef MapperTraits<TSpec, TConfig>        TTraits;
    typedef MatchesWriter<TSpec, TTraits>       TMatchesWriter;
    typedef typename TTraits::TMatchesSet       TMatchesSet;
//    typedef typename Value<TMatchesSet>::Type   TMatches;
    typedef typename Iterator<TMatchesSet, Standard>::Type  TMatchesIt;

    start(mapper.timer);

    // Sort each set of matches by errors.
//    forEach(mapper.anchorsSet, sortMatches<TMatches, SortErrors>, typename TConfig::TThreading());
    iterate(mapper.anchorsSet, sortMatches<TMatchesIt, SortErrors>, Standard(), typename TConfig::TThreading());

//    TMatchesWriter writer(mapper.outputStream, mapper.outputCtx,
//                          mapper.ctx, mapper.contigs, mapper.reads,
//                          mapper.anchorsSet, mapper.options);

    stop(mapper.timer);
    std::cout << "Output time:\t\t\t" << mapper.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function writeHits()
// ----------------------------------------------------------------------------
// Debug stuff.

template <typename TSpec, typename TConfig, typename TReadSeqs, typename THits, typename TSeeds, typename TFilename>
inline void writeHits(Mapper<TSpec, TConfig> const & mapper, TReadSeqs const & readSeqs,
                      THits const & hits, TSeeds const & seeds, TFilename const & filename)
{
    typedef MapperTraits<TSpec, TConfig>                TTraits;
    typedef typename Id<TSeeds>::Type                   TSeedId;
    typedef Pair<TSeedId>                               TSeedIds;
    typedef typename TTraits::THit                      THit;
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
            file << countHits<THitSize>(hits, seedHitIds, Serial()) << "\t";
        }
        file << "\n";
    }

    file.close();
}

// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void mapReads(Mapper<TSpec, TConfig> & mapper)
{
    _mapReadsImpl(mapper, mapper.reads._readSeqs, typename TConfig::TSequencing(), typename TConfig::TStrategy());
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

    collectSeeds<0>(mapper, readSeqs);
    findSeeds<0>(mapper, 0);
    classifyReads(mapper);
    collectSeeds<1>(mapper, readSeqs);
    collectSeeds<2>(mapper, readSeqs);
    findSeeds<1>(mapper, 1);
    findSeeds<2>(mapper, 2);
    extendHits(mapper);
    aggregateAnchors(mapper);
    compactAnchors(mapper, readSeqs);
    verifyAnchors(mapper, readSeqs);
    writeMatches(mapper);
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
    typedef MapperTraits<TSpec, TConfig>    TTraits;
    typedef typename TTraits::THit          THit;

    initReadsContext(mapper, readSeqs);
    initSeeds(mapper, readSeqs);
    clearHits(mapper);

    start(mapper.timer);
    collectSeeds<0>(mapper, readSeqs);
    findSeeds<0>(mapper, 0);
    classifyReads(mapper);
    collectSeeds<1>(mapper, readSeqs);
    collectSeeds<2>(mapper, readSeqs);
    findSeeds<0>(mapper, 1);
    findSeeds<0>(mapper, 2);

    // Sort hits by range size.
    // TODO(esiragusa): generalize sorting for approximate seeds, where one seed can have multiple hits.
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; ++bucketId)
        stableSort(mapper.hits[bucketId], HitsSorterByCount<THit>());

    extendHits(mapper);
    aggregateAnchors(mapper);

//    std::cout << "Mapped reads:\t\t\t" << countMapped(mapper.ctx, typename TConfig::TThreading()) << std::endl;

//    clearHits(mapper);
//    collectSeeds<1>(mapper, readSeqs);
//    collectSeeds<2>(mapper, readSeqs);
//    findSeeds<1>(mapper, 1);
//    findSeeds<1>(mapper, 2);
////    sortHits(mapper);
//    extendHits(mapper);
//
//    clearHits(mapper);
//    collectSeeds<2>(mapper, readSeqs);
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

    // Open and init output file.
    initOutput(mapper);

    // Process reads in blocks.
    while (!atEnd(mapper.readsLoader))
    {
        printRuler();

        loadReads(mapper);
        mapReads(mapper);
        clearReads(mapper);
    }

    // Close output file.
    close(mapper.outputStream);

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

template <typename TExecSpace,
          typename TThreading,
          typename TOutputFormat,
          typename TSequencing,
          typename TStrategy,
          typename TAnchoring>
inline void spawnMapper(Options const & options,
                        TExecSpace const & /* tag */,
                        TThreading const & /* tag */,
                        TOutputFormat const & /* tag */,
                        TSequencing const & /* tag */,
                        TStrategy const & /* tag */,
                        TAnchoring const & /* tag */)
{
    typedef ReadMapperConfig<TExecSpace,
                             TThreading,
                             TOutputFormat,
                             TSequencing,
                             TStrategy,
                             TAnchoring> TConfig;

    Mapper<void, TConfig> mapper(options);
    runMapper(mapper);
}

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_H_
