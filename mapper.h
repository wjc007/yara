// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
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

#ifndef APP_YARA_MAPPER_H_
#define APP_YARA_MAPPER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Mapper Options
// ----------------------------------------------------------------------------

struct Options
{
    typedef std::string             TString;
    typedef std::vector<TString>    TList;

    CharString          genomeFile;
    CharString          genomeIndexFile;
    Pair<CharString>    readsFile;
    CharString          outputFile;
    OutputFormat        outputFormat;
    TList               outputFormatList;
    TList               outputFormatExtensions;
    bool                outputSecondary;
    bool                outputHeader;

    MappingMode         mappingMode;
    unsigned            errorRate;
//    unsigned            strataRate;
    bool                quick;

    bool                singleEnd;
    unsigned            libraryLength;
    unsigned            libraryError;
    TList               libraryOrientationList;
    LibraryOrientation  libraryOrientation;
//    bool                anchorOne;

    unsigned            readsCount;
    bool                noCuda;
    unsigned            threadsCount;
    unsigned            hitsThreshold;
    unsigned            verbose;

    CharString          commandLine;
    CharString          version;

    Options() :
        outputFormat(SAM),
        outputSecondary(false),
        outputHeader(true),
        mappingMode(STRATA),
        errorRate(5),
//        strataRate(0),
        quick(false),
        singleEnd(true),
        libraryLength(200),
        libraryError(200),
        libraryOrientation(FWD_REV),
//        anchorOne(false),
        readsCount(100000),
        noCuda(false),
        threadsCount(1),
        hitsThreshold(300),
        verbose(0)
    {
        outputFormatList.push_back("sam");
        outputFormatExtensions.push_back("sam");
#ifdef SEQAN_HAS_ZLIB
        outputFormatList.push_back("bam");
        outputFormatExtensions.push_back("bam");
#endif

        libraryOrientationList.push_back("fwd-rev");
        libraryOrientationList.push_back("fwd-fwd");
        libraryOrientationList.push_back("rev-rev");
    }
};

// ----------------------------------------------------------------------------
// Mapper Configuration
// ----------------------------------------------------------------------------

template <typename TExecSpace_      = ExecHost,
          typename TThreading_      = Parallel,
          typename TOutputFormat_   = Sam,
          typename TSequencing_     = SingleEnd,
          typename TStrategy_       = Strata,
          typename TAnchoring_      = Nothing,
          unsigned BUCKETS_         = 3>
struct ReadMapperConfig : public ContigsConfig<YaraStringSpec>, public ReadsConfig<void>
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

    typedef Index<YaraContigsFM, YaraIndexSpec>                     THostIndex;
    typedef typename Space<THostIndex, TExecSpace>::Type            TIndex;
    typedef typename Size<TIndex>::Type                             TIndexSize;
    typedef typename Fibre<TIndex, FibreSA>::Type                   TSA;

    typedef Reads<TSequencing, TConfig>                             TReads;
    typedef ReadsLoader<TSequencing, TConfig>                       TReadsLoader;
    typedef typename TReads::TReadSeqs                              THostReadSeqs;
    typedef typename Space<THostReadSeqs, TExecSpace>::Type         TReadSeqs;
    typedef typename Value<TReadSeqs>::Type                         TReadSeq;
    typedef typename Size<TReadSeqs>::Type                          TReadSeqsSize;
    typedef String<TReadSeqsSize>                                   TSeedsCount;

    typedef typename TContigs::TContigNames                         TContigNames;
    typedef typename TContigs::TContigNamesCache                    TContigNamesCache;
    typedef Stream<FileStream<File<> > >                            TOutputStream;
    typedef BamIOContext<TContigNames, TContigNamesCache>           TOutputContext;

    typedef ReadsContext<TSpec, TConfig>                            TReadsContext;

    typedef StringSet<TReadSeqs, Segment<TReadSeqs> >               TSeeds;
    typedef Tuple<TSeeds, TConfig::BUCKETS>                         TSeedsBuckets;

    typedef Hit<TIndexSize, HammingDistance>                        THit;
    typedef String<THit>                                            THits;
    typedef Tuple<THits, TConfig::BUCKETS>                          THitsBuckets;
    typedef String<TIndexSize>                                      THitsCounts;

    typedef StringSet<TSeedsCount, Owner<ConcatDirect<> > >         TRanks;
    typedef Tuple<TRanks, TConfig::BUCKETS>                         TRanksBuckets;

    typedef Match<void>                                             TMatch;
    typedef String<TMatch>                                          TMatches;
    typedef StringSet<TMatches, Segment<TMatches> >                 TMatchesSet;

    typedef String<CigarElement<> >                                 TCigar;
    typedef StringSet<TCigar, Segment<TCigar> >                     TCigarSet;
    typedef StringSetLimits<TCigarSet>::Type                        TCigarLimits;

    typedef Multiple<FinderSTree>                                   TAlgorithmExt;
    typedef Multiple<Backtracking<HammingDistance> >                TAlgorithmApx;
    typedef Pattern<TSeeds, TAlgorithmExt>                          TPatternExt;
    typedef Pattern<TSeeds, TAlgorithmApx>                          TPatternApx;
    typedef Finder2<TIndex, TPatternExt, TAlgorithmExt>             TFinderExt;
    typedef Finder2<TIndex, TPatternApx, TAlgorithmApx>             TFinderApx;
};

// ----------------------------------------------------------------------------
// Mapper Stats
// ----------------------------------------------------------------------------

template <typename TValue>
struct Stats
{
    TValue loadGenome;
    TValue loadReads;
    TValue collectSeeds;
    TValue findSeeds;
    TValue classifyReads;
    TValue rankSeeds;
    TValue extendHits;
    TValue sortMatches;
    TValue compactMatches;
    TValue selectPairs;
    TValue alignMatches;
    TValue writeMatches;

    unsigned long loadedReads;
    unsigned long mappedReads;
    unsigned long pairedReads;

    Stats() :
        loadGenome(0),
        loadReads(0),
        collectSeeds(0),
        findSeeds(0),
        classifyReads(0),
        rankSeeds(0),
        extendHits(0),
        sortMatches(0),
        compactMatches(0),
        selectPairs(0),
        alignMatches(0),
        writeMatches(0),
        loadedReads(0),
        mappedReads(0),
        pairedReads(0)
    {}
};

// ----------------------------------------------------------------------------
// Class Mapper
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig = void>
struct Mapper
{
    typedef MapperTraits<TSpec, TConfig>    Traits;

    Options const &                     options;
    Timer<double>                       timer;
    Stats<double>                       stats;

    typename Traits::TContigs           contigs;
    typename Traits::TIndex             index;
    typename Traits::TReads             reads;
    typename Traits::TReadsLoader       readsLoader;
    typename Traits::TOutputStream      outputStream;
    typename Traits::TOutputContext     outputCtx;

    typename Traits::TReadsContext      ctx;
    typename Traits::TSeedsBuckets      seeds;
    typename Traits::THitsBuckets       hits;
    typename Traits::TRanksBuckets      ranks;

    typename Traits::TMatches           matches;
    typename Traits::TMatchesSet        matchesSet;
    typename Traits::TMatches           primaryMatches;

    typename Traits::TCigar             cigars;
    typename Traits::TCigarSet          cigarSet;

    typename Traits::TFinderExt         finderExt;
    typename Traits::TFinderApx         finderApx;

    Mapper(Options const & options) :
        options(options),
        stats(),
        contigs(),
        index(),
        reads(),
        readsLoader(),
        outputStream(),
        outputCtx(contigs.names, contigs.namesCache),
        ctx(),
        seeds(),
        hits(),
        ranks(),
        matches(),
        matchesSet(),
        primaryMatches(),
        cigars(),
        cigarSet(),
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
inline void configureThreads(Mapper<TSpec, TConfig> & me)
{
#ifdef _OPENMP
    omp_set_num_threads(me.options.threadsCount);

    if (me.options.verbose > 0)
        std::cout << "Threads count:\t\t\t" << omp_get_max_threads() << std::endl;
#else
    ignoreUnusedVariableWarning(me);
#endif
}

// ----------------------------------------------------------------------------
// Function loadGenome()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void loadGenome(Mapper<TSpec, TConfig> & me)
{
    start(me.timer);
    if (!open(me.contigs, toCString(me.options.genomeIndexFile)))
        throw RuntimeError("Error while opening genome file.");
    stop(me.timer);
    me.stats.loadGenome += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cout << "Loading genome:\t\t\t" << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function loadGenomeIndex()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void loadGenomeIndex(Mapper<TSpec, TConfig> & me)
{
#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif

    start(me.timer);
    if (!open(me.index, toCString(me.options.genomeIndexFile)))
        throw RuntimeError("Error while opening genome index file.");
    stop(me.timer);
    me.stats.loadGenome += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cout << "Loading genome index:\t\t" << me.timer << std::endl;

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif
}

// ----------------------------------------------------------------------------
// Function openReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void openReads(Mapper<TSpec, TConfig> & me)
{
    _openReadsImpl(me, typename TConfig::TSequencing());
}

template <typename TSpec, typename TConfig, typename TSequencing>
inline void _openReadsImpl(Mapper<TSpec, TConfig> & me, TSequencing const & /*tag */)
{
    open(me.readsLoader, me.options.readsFile);
}

template <typename TSpec, typename TConfig>
inline void _openReadsImpl(Mapper<TSpec, TConfig> & me, SingleEnd const & /* tag */)
{
    open(me.readsLoader, me.options.readsFile.i1);
}

// ----------------------------------------------------------------------------
// Function loadReads()
// ----------------------------------------------------------------------------
// Loads one block of reads.

template <typename TSpec, typename TConfig>
inline void loadReads(Mapper<TSpec, TConfig> & me)
{
    start(me.timer);
    clear(me.reads);
    load(me.reads, me.readsLoader, me.options.readsCount);
    stop(me.timer);
    me.stats.loadReads += getValue(me.timer);
    me.stats.loadedReads += getReadsCount(me.reads.seqs);

    if (me.options.verbose > 1)
    {
        std::cout << "Loading reads:\t\t\t" << me.timer << std::endl;
        std::cout << "Reads count:\t\t\t" << getReadsCount(me.reads.seqs) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function clearReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void clearReads(Mapper<TSpec, TConfig> & me)
{
    clear(me.reads);
}

// ----------------------------------------------------------------------------
// Function initOutput()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void initOutput(Mapper<TSpec, TConfig> & me)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;

    if (!open(me.outputStream, toCString(me.options.outputFile), OPEN_RDWR | OPEN_CREATE))
        throw RuntimeError("Error while opening output file.");

    if (me.options.outputHeader)
    {
        BamHeader header;

        // Fill header.
        fillHeader(header, me.options, me.contigs.seqs, me.contigs.names);

        // Write header to stream.
        write2(me.outputStream, header, me.outputCtx, typename TTraits::TOutputFormat());
    }
}

// ----------------------------------------------------------------------------
// Function initSeeds()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void initSeeds(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs)
{
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
        setHost(me.seeds[bucketId], readSeqs);
}

// ----------------------------------------------------------------------------
// Function clearSeeds()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void clearSeeds(Mapper<TSpec, TConfig> & me)
{
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
    {
        clear(me.seeds[bucketId]);
        clear(me.ranks[bucketId]);
        shrinkToFit(me.seeds[bucketId]);
        shrinkToFit(me.ranks[bucketId]);
    }
}

// ----------------------------------------------------------------------------
// Function initReadsContext()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void initReadsContext(Mapper<TSpec, TConfig> & me, TReadSeqs const & readSeqs)
{
    clear(me.ctx);
    resize(me.ctx, readSeqs);
}

// ----------------------------------------------------------------------------
// Function collectSeeds()
// ----------------------------------------------------------------------------
// Collects seeds from all reads.

template <unsigned ERRORS, typename TSpec, typename TConfig, typename TReadSeqs>
inline void collectSeeds(Mapper<TSpec, TConfig> & me, TReadSeqs const & readSeqs)
{
    typedef MapperTraits<TSpec, TConfig>                TTraits;
    typedef SeedsCollector<Counter, TTraits>            TCounter;
    typedef SeedsCollector<void, TTraits>               TFiller;

    typename TTraits::TSeedsCount seedsCounts;

    start(me.timer);
    TCounter counter(me.ctx, me.seeds[ERRORS], seedsCounts, ERRORS, readSeqs, me.options);
    TFiller filler(me.ctx, me.seeds[ERRORS], seedsCounts, ERRORS, readSeqs, me.options);
    stop(me.timer);
    me.stats.collectSeeds += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cout << "Seeding time:\t\t\t" << me.timer << std::endl;
        std::cout << "Seeds count:\t\t\t" << length(me.seeds[ERRORS]) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function findSeeds()
// ----------------------------------------------------------------------------

template <unsigned ERRORS, typename TSpec, typename TConfig, typename TBucketId>
inline void findSeeds(Mapper<TSpec, TConfig> & me, TBucketId bucketId)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;
    typedef typename TTraits::TPatternExt           TPatternExt;
    typedef typename TTraits::TPatternApx           TPatternApx;

    start(me.timer);
    if (ERRORS > 0)
    {
        setScoreThreshold(me.finderApx, ERRORS);
        // TODO(esiragusa): guess the number of hits.
        reserve(me.hits[bucketId], lengthSum(me.seeds[bucketId]) * Power<ERRORS, 2>::VALUE, Exact());
        _findSeedsImpl(me, me.hits[bucketId], me.seeds[bucketId], me.finderApx, TPatternApx());
    }
    else
    {
        // TODO(esiragusa): guess the number of hits.
        reserve(me.hits[bucketId], length(me.seeds[bucketId]), Exact());
        _findSeedsImpl(me, me.hits[bucketId], me.seeds[bucketId], me.finderExt, TPatternExt());
    }
    stop(me.timer);
    me.stats.findSeeds += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cout << "Filtering time:\t\t\t" << me.timer << std::endl;
        std::cout << "Hits count:\t\t\t" <<
               countHits<unsigned long>(me.hits[bucketId], typename TConfig::TThreading()) << std::endl;
    }
}

template <typename TSpec, typename TConfig, typename THits, typename TSeeds, typename TFinder, typename TPattern>
inline void _findSeedsImpl(Mapper<TSpec, TConfig> & /* me */, THits & hits, TSeeds & seeds, TFinder & finder, TPattern)
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
inline void classifyReads(Mapper<TSpec, TConfig> & me)
{
    typedef MapperTraits<TSpec, TConfig>                TTraits;
    typedef ReadsClassifier<TSpec, TTraits>             TClassifier;

    start(me.timer);
    TClassifier classifier(me.ctx, me.hits[0], me.seeds[0], me.options);
    stop(me.timer);
    me.stats.classifyReads += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cout << "Classification time:\t\t" << me.timer << std::endl;
        std::cout << "Hits count:\t\t\t" <<
               countHits<unsigned long>(me.hits[0], typename TConfig::TThreading()) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function rankSeeds()
// ----------------------------------------------------------------------------
// Rank the seeds in all buckets.

template <typename TSpec, typename TConfig>
inline void rankSeeds(Mapper<TSpec, TConfig> & me)
{
    typedef MapperTraits<TSpec, TConfig>    TTraits;
    typedef SeedsRanker<TSpec, TTraits>     TSeedsRanker;

    typename TTraits::THitsCounts hitsCounts;

    start(me.timer);
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
        TSeedsRanker ranker(hitsCounts, me.ranks[bucketId], me.seeds[bucketId], me.hits[bucketId], me.options);
    stop(me.timer);
    me.stats.rankSeeds += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cout << "Ranking time:\t\t\t" << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function clearHits()
// ----------------------------------------------------------------------------
// Clears the hits in all buckets.

template <typename TSpec, typename TConfig>
inline void clearHits(Mapper<TSpec, TConfig> & me)
{
    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
    {
        clear(me.hits[bucketId]);
        shrinkToFit(me.hits[bucketId]);
    }
}

// ----------------------------------------------------------------------------
// Function countHits()
// ----------------------------------------------------------------------------
// Counts the hits in all buckets.

template <typename TSpec, typename TConfig>
inline unsigned long countHits(Mapper<TSpec, TConfig> const & me)
{
    unsigned long hitsCount = 0;

    for (unsigned bucketId = 0; bucketId < TConfig::BUCKETS; bucketId++)
        hitsCount += countHits<unsigned long>(me.hits[bucketId], typename TConfig::TThreading());

    return hitsCount;
}

// ----------------------------------------------------------------------------
// Function extendHits()
// ----------------------------------------------------------------------------
// Extends the hits in a bucket.

template <unsigned ERRORS, typename TSpec, typename TConfig, typename TBucketId>
inline void extendHits(Mapper<TSpec, TConfig> & me, TBucketId bucketId)
{
    typedef MapperTraits<TSpec, TConfig>    TTraits;
    typedef HitsExtender<TSpec, TTraits>    THitsExtender;

    start(me.timer);
    THitsExtender extender(me.ctx, me.matches, me.contigs.seqs,
                           me.seeds[bucketId], me.hits[bucketId], me.ranks[bucketId], ERRORS,
                           indexSA(me.index), me.options);
    stop(me.timer);
    me.stats.extendHits += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cout << "Extension time:\t\t\t" << me.timer << std::endl;
        std::cout << "Matches count:\t\t\t" << length(me.matches) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function reserveMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void reserveMatches(Mapper<TSpec, TConfig> & me)
{
    // TODO(esiragusa): guess the number of matches.
    reserve(me.matches, countHits(me) / 3);
}

// ----------------------------------------------------------------------------
// Function aggregateMatches()
// ----------------------------------------------------------------------------
// Aggregate matches by readId.

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void aggregateMatches(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs)
{
    typedef MapperTraits<TSpec, TConfig>    TTraits;
    typedef typename TTraits::TMatch        TMatch;

    // Bucket sort matches by readId.
    start(me.timer);
    setHost(me.matchesSet, me.matches);
    sort(me.matches, MatchSorter<TMatch, SortReadId>(), typename TConfig::TThreading());
    bucket(me.matchesSet, Getter<TMatch, SortReadId>(), getReadsCount(readSeqs), typename TConfig::TThreading());
    stop(me.timer);
    me.stats.sortMatches += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cout << "Sorting time:\t\t\t" << me.timer << std::endl;

    start(me.timer);
    removeDuplicates(me.matchesSet, typename TConfig::TThreading());
    stop(me.timer);
    me.stats.compactMatches += getValue(me.timer);

    if (me.options.verbose > 1)
    {
        std::cout << "Compaction time:\t\t" << me.timer << std::endl;
        std::cout << "Matches count:\t\t\t" << lengthSum(me.matchesSet) << std::endl;
    }
}

// ----------------------------------------------------------------------------
// Function verifyMatches()
// ----------------------------------------------------------------------------
// Verifies all mates in within the insert window of their matches.

//template <typename TSpec, typename TConfig, typename TReadSeqs>
//inline void verifyMatches(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs)
//{
//    _verifyMatchesImpl(me, readSeqs, typename TConfig::TAnchoring());
//}
//
//template <typename TSpec, typename TConfig, typename TReadSeqs, typename TAnchoring>
//inline void _verifyMatchesImpl(Mapper<TSpec, TConfig> & /* me */, TReadSeqs & /* readSeqs */, TAnchoring) {}
//
//template <typename TSpec, typename TConfig, typename TReadSeqs>
//inline void _verifyMatchesImpl(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs, AnchorOne)
//{
//    typedef MapperTraits<TSpec, TConfig>    TTraits;
//    typedef AnchorsVerifier<TSpec, TTraits> TMatchesVerifier;
//
//    start(me.timer);
//    TMatchesVerifier verifier(me.ctx, me.pairs,
//                              me.contigs.seqs, readSeqs,
//                              me.matchesSet, me.options);
//    stop(me.timer);
//
//    if (me.options.verbose > 1)
//    {
//        std::cout << "Verification time:\t\t" << me.timer << std::endl;
//        std::cout << "Mates count:\t\t\t" << length(me.pairs) << std::endl;
//        std::cout << "Mapped pairs:\t\t\t" <<
//                countMappedReads(readSeqs, me.pairs, typename TConfig::TThreading()) << std::endl;
//    }
//}

// ----------------------------------------------------------------------------
// Function clearMatches()
// ----------------------------------------------------------------------------
// Clears all matches.

template <typename TSpec, typename TConfig>
inline void clearMatches(Mapper<TSpec, TConfig> & me)
{
    clear(me.matchesSet);
    clear(me.matches);
    shrinkToFit(me.matches);

    clear(me.primaryMatches);
    shrinkToFit(me.primaryMatches);
}

// ----------------------------------------------------------------------------
// Function _sortMatchesImpl()
// ----------------------------------------------------------------------------
// Sort each set of matches by errors.

template <typename TSpec, typename TConfig>
inline void _sortMatchesImpl(Mapper<TSpec, TConfig> & me)
{
    typedef MapperTraits<TSpec, TConfig>                    TTraits;
    typedef typename TTraits::TMatchesSet                   TMatchesSet;
    typedef typename Iterator<TMatchesSet, Standard>::Type  TMatchesIt;
//    typedef typename Value<TMatchesSet>::Type               TMatches;

    start(me.timer);
    iterate(me.matchesSet, sortMatches<TMatchesIt, SortErrors>, Standard(), typename TTraits::TThreading());
//    forEach(me.matchesSet, sortMatches<TMatches, SortErrors>, typename TTraits::TThreading());
    stop(me.timer);
    me.stats.sortMatches += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cout << "Sorting time:\t\t\t" << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function rankMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void rankMatches(Mapper<TSpec, TConfig> & me, TReadSeqs const & readSeqs)
{
    _rankMatchesImpl(me, readSeqs, typename TConfig::TSequencing());
}

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _rankMatchesImpl(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs, SingleEnd)
{
    typedef MapperTraits<TSpec, TConfig>        TTraits;
    typedef typename TTraits::TMatchesSet       TMatchesSet;
    typedef typename Value<TMatchesSet>::Type   TMatches;

    // Sort matches by errors.
    _sortMatchesImpl(me);

    resize(me.primaryMatches, getReadsCount(readSeqs), Exact());
    transform(me.primaryMatches, me.matchesSet, getFirstMatch<TMatches>, typename TTraits::TThreading());

    unsigned long mappedReads = 0;
    if (me.options.verbose > 0)
    {
        mappedReads = countValidMatches(me.primaryMatches, typename TTraits::TThreading());
        me.stats.mappedReads += mappedReads;
    }
    if (me.options.verbose > 1)
        std::cout << "Mapped reads:\t\t\t" << mappedReads << std::endl;
}

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _rankMatchesImpl(Mapper<TSpec, TConfig> & me, TReadSeqs const & readSeqs, PairedEnd)
{
    typedef MapperTraits<TSpec, TConfig>    TTraits;
    typedef PairsSelector<TSpec, TTraits>   TPairsSelector;
    typedef MatchesSelector<TSpec, TTraits> TMatchesSelector;

    unsigned long pairedReads = 0;

    start(me.timer);
    TPairsSelector selector(me.primaryMatches, readSeqs, me.matchesSet, me.options);
    stop(me.timer);
    me.stats.selectPairs += getValue(me.timer);

    if (me.options.verbose > 0)
    {
        pairedReads = countValidMatches(me.primaryMatches, typename TTraits::TThreading());
        me.stats.pairedReads += pairedReads;
    }
    if (me.options.verbose > 1)
    {
        std::cout << "Pairing time:\t\t\t" << me.timer << std::endl;
        std::cout << "Mapped pairs:\t\t\t" << pairedReads / 2 << std::endl;
    }

    // Mark paired reads.
    transform(me.ctx.paired, me.primaryMatches, isValid<void>, typename TTraits::TThreading());

    // Sort matches by errors.
    _sortMatchesImpl(me);

    // Select best match for unpaired reads.
    TMatchesSelector mselector(me.primaryMatches, me.ctx, me.matchesSet);

    unsigned long mappedReads = 0;
    if (me.options.verbose > 0)
    {
        mappedReads = countValidMatches(me.primaryMatches, typename TTraits::TThreading());
        me.stats.mappedReads += mappedReads;
    }
    if (me.options.verbose > 1)
        std::cout << "Mapped reads:\t\t\t" << mappedReads << std::endl;
}

// ----------------------------------------------------------------------------
// Function alignMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void alignMatches(Mapper<TSpec, TConfig> & me)
{
    typedef MapperTraits<TSpec, TConfig>        TTraits;
    typedef MatchesAligner<TSpec, TTraits>      TMatchesAligner;

    start(me.timer);
    setHost(me.cigarSet, me.cigars);
    typename TTraits::TCigarLimits cigarLimits;
    TMatchesAligner aligner(me.cigarSet, cigarLimits, me.primaryMatches, me.contigs.seqs, me.reads.seqs, me.options);
    stop(me.timer);
    me.stats.alignMatches += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cout << "Alignment time:\t\t\t" << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function clearAlignments()
// ----------------------------------------------------------------------------
// Clears all cigars.

template <typename TSpec, typename TConfig>
inline void clearAlignments(Mapper<TSpec, TConfig> & me)
{
    clear(me.cigars);
    clear(me.cigarSet);
    shrinkToFit(me.cigarSet);
}

// ----------------------------------------------------------------------------
// Function writeMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void writeMatches(Mapper<TSpec, TConfig> & me)
{
    typedef MapperTraits<TSpec, TConfig>        TTraits;
    typedef MatchesWriter<TSpec, TTraits>       TMatchesWriter;

    start(me.timer);
    TMatchesWriter writer(me.outputStream, me.outputCtx,
                          me.matchesSet, me.primaryMatches, me.cigarSet,
                          me.ctx, me.contigs, me.reads,
                          me.options);
    stop(me.timer);
    me.stats.writeMatches += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cout << "Output time:\t\t\t" << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void mapReads(Mapper<TSpec, TConfig> & me)
{
    _mapReadsImpl(me, me.reads.seqs, typename TConfig::TStrategy());
}

// ----------------------------------------------------------------------------
// Function _mapReadsImpl(); All
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _mapReadsImpl(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs, All)
{
    initReadsContext(me, readSeqs);
    initSeeds(me, readSeqs);

    collectSeeds<0>(me, readSeqs);
    findSeeds<0>(me, 0);
    classifyReads(me);
    collectSeeds<1>(me, readSeqs);
    collectSeeds<2>(me, readSeqs);
    findSeeds<1>(me, 1);
    if (me.options.quick)
        findSeeds<1>(me, 2);
    else
        findSeeds<2>(me, 2);
    reserveMatches(me);
    extendHits<0>(me, 0);
    extendHits<1>(me, 1);
    extendHits<2>(me, 2);
    clearSeeds(me);
    clearHits(me);
    aggregateMatches(me, readSeqs);
//    verifyMatches(me, readSeqs);
    rankMatches(me, readSeqs);
    alignMatches(me);
    writeMatches(me);
    clearMatches(me);
    clearAlignments(me);
}

// ----------------------------------------------------------------------------
// Function _mapReadsImpl(); Strata
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _mapReadsImpl(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs, Strata)
{
    initReadsContext(me, readSeqs);
    initSeeds(me, readSeqs);

    collectSeeds<0>(me, readSeqs);
    findSeeds<0>(me, 0);
    classifyReads(me);
    collectSeeds<1>(me, readSeqs);
    collectSeeds<2>(me, readSeqs);
    findSeeds<0>(me, 1);
    findSeeds<0>(me, 2);
    rankSeeds(me);
    reserveMatches(me);
    extendHits<0>(me, 0);
    extendHits<0>(me, 1);
    extendHits<0>(me, 2);
    clearSeeds(me);
    clearHits(me);

    initSeeds(me, readSeqs);
    collectSeeds<1>(me, readSeqs);
    findSeeds<1>(me, 1);
    collectSeeds<2>(me, readSeqs);
    findSeeds<1>(me, 2);
    rankSeeds(me);
    // TODO(esiragusa): filter out hits with distance < 1.
    extendHits<1>(me, 1);
    extendHits<1>(me, 2);
    clearSeeds(me);
    clearHits(me);

    if (!me.options.quick)
    {
        initSeeds(me, readSeqs);
        collectSeeds<2>(me, readSeqs);
        findSeeds<2>(me, 2);
        rankSeeds(me);
    // TODO(esiragusa): filter out hits with distance < 2.
        extendHits<2>(me, 2);
        clearHits(me);
        clearSeeds(me);
    }

    aggregateMatches(me, readSeqs);
//    verifyMatches(me, readSeqs);
    rankMatches(me, readSeqs);
    alignMatches(me);
    writeMatches(me);
    clearMatches(me);
    clearAlignments(me);
}

// ----------------------------------------------------------------------------
// Function printStats()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TValue>
inline void printStats(Mapper<TSpec, TConfig> const & me, Timer<TValue> const & timer)
{
    printRuler(std::cout);

    TValue total = getValue(timer) / 100;

    std::cout << "Total time:\t\t\t" << getValue(timer) << " sec" << std::endl;
    std::cout << "Genome loading time:\t\t" << me.stats.loadGenome << " sec" << "\t\t" << me.stats.loadGenome / total << " %" << std::endl;
    std::cout << "Reads loading time:\t\t" << me.stats.loadReads << " sec" << "\t\t" << me.stats.loadReads / total << " %" << std::endl;
    std::cout << "Seeding time:\t\t\t" << me.stats.collectSeeds << " sec" << "\t\t" << me.stats.collectSeeds / total << " %" << std::endl;
    std::cout << "Filtering time:\t\t\t" << me.stats.findSeeds << " sec" << "\t\t" << me.stats.findSeeds / total << " %" << std::endl;
    std::cout << "Classification time:\t\t" << me.stats.classifyReads << " sec" << "\t\t" << me.stats.classifyReads / total << " %" << std::endl;
    std::cout << "Ranking time:\t\t\t" << me.stats.rankSeeds << " sec" << "\t\t" << me.stats.rankSeeds / total << " %" << std::endl;
    std::cout << "Extension time:\t\t\t" << me.stats.extendHits << " sec" << "\t\t" << me.stats.extendHits / total << " %" << std::endl;
    std::cout << "Sorting time:\t\t\t" << me.stats.sortMatches << " sec" << "\t\t" << me.stats.sortMatches / total << " %" << std::endl;
    std::cout << "Compaction time:\t\t" << me.stats.compactMatches << " sec" << "\t\t" << me.stats.compactMatches / total << " %" << std::endl;
    if (IsSameType<typename TConfig::TSequencing, PairedEnd>::VALUE)
        std::cout << "Pairing time:\t\t\t" << me.stats.selectPairs << " sec" << "\t\t" << me.stats.selectPairs / total << " %" << std::endl;
    std::cout << "Alignment time:\t\t\t" << me.stats.alignMatches << " sec" << "\t\t" << me.stats.alignMatches / total << " %" << std::endl;
    std::cout << "Output time:\t\t\t" << me.stats.writeMatches << " sec" << "\t\t" << me.stats.writeMatches / total << " %" << std::endl;

    printRuler(std::cout);

    double totalReads = me.stats.loadedReads / 100;
    std::cout << "Total reads:\t\t\t" << me.stats.loadedReads << std::endl;
    std::cout << "Mapped reads:\t\t\t" << me.stats.mappedReads << "\t\t" << me.stats.mappedReads / totalReads << " %" << std::endl;
    if (IsSameType<typename TConfig::TSequencing, PairedEnd>::VALUE)
        std::cout << "Paired reads:\t\t\t" << me.stats.pairedReads << "\t\t" << me.stats.pairedReads / totalReads << " %" << std::endl;
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void runMapper(Mapper<TSpec, TConfig> & me)
{
    Timer<double> timer;

    start(timer);

    configureThreads(me);

    if (me.options.verbose > 1) printRuler(std::cout);

    loadGenome(me);
    loadGenomeIndex(me);

    // Open reads file.
    openReads(me);

    // Open and init output file.
    initOutput(me);

    // Process reads in blocks.
    while (!atEnd(me.readsLoader))
    {
        if (me.options.verbose > 1) printRuler(std::cout);

        loadReads(me);
        mapReads(me);
        clearReads(me);
    }

    // Close output file.
    close(me.outputStream);

    // Close reads file.
    close(me.readsLoader);

    stop(timer);

    if (me.options.verbose > 0)
        printStats(me, timer);
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

/*
template <typename TSpec, typename TConfig>
inline void runMapper(Mapper<TSpec, TConfig> & me, Parallel)
{
    typedef Mapper<TSpec, TConfig>                      TMapper;
    typedef Reads<void, typename TMapper::TReadsConfig> TReads;
    typedef Logger<std::ostream>                        TLogger;

    Timer<double> timer;
    TLogger cout(std::cout);
    TLogger cerr(std::cerr);

#ifdef _OPENMP
    cout << "Threads count:\t\t\t" << me.options.threadsCount << std::endl;
#endif

#ifdef _OPENMP
    // Disable nested parallelism.
    omp_set_nested(false);
#endif

    loadGenome(me);
    loadGenomeIndex(me);

    // Open reads file.
    open(me.readsLoader);

    // Process reads in parallel.
    SEQAN_OMP_PRAGMA(parallel firstprivate(timer) num_threads(3))
    {
        // Reserve space for reads.
        TReads reads;
        reserve(reads, me.options.readsCount);

        // Process reads.
        while (true)
        {
            // Load a block of reads.
            SEQAN_OMP_PRAGMA(critical(_mapper_readsLoader_load))
            {
                // No more reads.
                if (!atEnd(me.readsLoader))
                {
                    start(me.timer);
                    setReads(me.readsLoader, reads);
                    load(me.readsLoader, me.options.readsCount);
                    stop(me.timer);

                    cout << "Loading reads:\t\t\t" << me.timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;// <<
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
                omp_set_num_threads(me.options.threadsCount);
                #endif

                start(me.timer);
                mapReads(me, me.options, getSeqs(reads));
                stop(me.timer);

                cout << "Mapping reads:\t\t\t" << me.timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;

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
                start(me.timer);
                sleep(reads.readsCount / 1000000.0);
                stop(me.timer);
                
                cout << "Writing results:\t\t" << me.timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;
            }

            // Clear mapped reads.
            clear(reads);
        }
    }

    // Close reads file.
    close(me.readsLoader);
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

#endif  // #ifndef APP_YARA_MAPPER_H_
