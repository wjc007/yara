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
    bool                singleEnd;
    bool                anchorOne;
    unsigned            libraryLength;
    unsigned            libraryError;

    unsigned            mappingBlock;
    bool                noCuda;
    unsigned            threadsCount;
    unsigned            hitsThreshold;

    Options() :
        errorRate(5),
        singleEnd(true),
        anchorOne(false),
        libraryLength(220),
        libraryError(50),
        mappingBlock(200000),
        noCuda(false),
        threadsCount(1),
        hitsThreshold(300)
    {}
};

// ----------------------------------------------------------------------------
// Enum ReadStatus
// ----------------------------------------------------------------------------

enum ReadStatus { STATUS_UNSEEDED, STATUS_SEEDED, STATUS_MAPPED, STATUS_UNMAPPABLE };
//enum ReadAnchor { ANCHOR_FIRST, ANCHOR_SECOND };

// ----------------------------------------------------------------------------
// Class ReadInfo
// ----------------------------------------------------------------------------

template <typename TConfig = void>
struct ReadInfo
{
    unsigned char seedErrors    : 2;
    ReadStatus    status        : 2;
//    ReadAnchor    anchor        : 1;

    ReadInfo() :
        seedErrors(0),
        status(STATUS_UNSEEDED)
    {};
};

// ----------------------------------------------------------------------------
// Mapper Configuration
// ----------------------------------------------------------------------------

template <typename TExecSpace_  = ExecHost,
          typename TSequencing_ = SingleEnd,
          typename TStrategy_   = AnyBest,
          typename TAnchoring_  = Nothing>
struct ReadMapperConfig : public CUDAStoreConfig
{
    typedef TExecSpace_     TExecSpace;
    typedef TSequencing_    TSequencing;
    typedef TStrategy_      TStrategy;
    typedef TAnchoring_     TAnchoring;
};

// ----------------------------------------------------------------------------
// Class Mapper
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig = void>
struct Mapper
{
    typedef Genome<void, TConfig>                                   TGenome;
    typedef GenomeLoader<void, TConfig>                             TGenomeLoader;
    typedef typename Contigs<TGenome>::Type                         TContigs;
    typedef typename Value<TContigs>::Type                          TContig;
    typedef typename StringSetPosition<TContigs>::Type              TContigsPos;

    typedef typename TConfig::TExecSpace                            TExecSpace;
    typedef Index<TFMContigs, TGenomeIndexSpec>                     THostIndex;
    typedef typename Space<THostIndex, TExecSpace>::Type            TIndex;

    typedef typename TConfig::TSequencing                           TSequencing;
    typedef FragmentStore<void, TConfig>                            TStore;
    typedef ReadsConfig<False, False, True, True, TConfig>          TReadsConfig;
    typedef Reads<TSequencing, TReadsConfig>                        TReads;
    typedef ReadsLoader<TSequencing, TReadsConfig>                  TReadsLoader;
    typedef typename TStore::TReadSeqStore                          THostReadSeqs;
    typedef typename Space<THostReadSeqs, TExecSpace>::Type         TReadSeqs;
    typedef typename Value<TReadSeqs>::Type                         TReadSeq;

    typedef String<ReadInfo<TConfig> >                              TReadInfoString;

    typedef StringSet<TReadSeqs, Segment<TReadSeqs> >               TSeedsSet;
    typedef Tuple<TSeedsSet, 3>                                     TSeeds;

    typedef typename Size<TIndex>::Type                             TIndexSize;
    typedef Hit<TIndexSize, HammingDistance>                        THit;
    typedef String<THit>                                            THitsString;
    typedef Tuple<THitsString, 3>                                   THits;

    typedef Match<void>                                             TMatch;
    typedef String<TMatch>                                          TMatches;

    typedef Multiple<FinderSTree>                                   TSeedingExt;
    typedef Multiple<Backtracking<HammingDistance> >                TSeedingApx;
    typedef Pattern<TSeedsSet, TSeedingExt>                         TSeedsExt;
    typedef Pattern<TSeedsSet, TSeedingApx>                         TSeedsApx;
    typedef Finder2<TIndex, TSeedsExt, TSeedingExt>                 TSeederExt;
    typedef Finder2<TIndex, TSeedsApx, TSeedingApx>                 TSeederApx;

    typedef AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_> TMyersSpec;
    typedef Myers<TMyersSpec, True, void>                           TExtenderAlgorithm;
    typedef Extender<TContigs, TReadSeq, TExtenderAlgorithm>        TExtender;

    typedef Verifier<TContigs, TReadSeq, Myers<> >                  TVerifier;
//    typedef Verifier<TContigs, TReadSeq, Filter<MultipleShiftAnd> > TVerifier;

//    typedef WriterConfig<Options, TReadSeqs>                        TWriterConfig;
//    typedef Writer<TSpec, TWriterConfig>                       TWriter;

    Timer<double>       timer;
    Options const &     options;

    TGenome             genome;
    TGenomeLoader       genomeLoader;
    TIndex              index;
    TStore              store;
    TReads              reads;
    TReadsLoader        readsLoader;

    TReadInfoString     info;
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
        setHost(seeds[1], getSeqs(reads));
        setHost(seeds[2], getSeqs(reads));
    };
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function configureThreads()
// ----------------------------------------------------------------------------
// Sets the number of threads that OpenMP can spawn.

template <typename TSpec, typename TConfig>
void configureThreads(Mapper<TSpec, TConfig> & mapper)
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
void loadGenome(Mapper<TSpec, TConfig> & mapper)
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
void loadGenomeIndex(Mapper<TSpec, TConfig> & mapper)
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
void openReads(Mapper<TSpec, TConfig> & mapper)
{
    _openReads(mapper, typename TConfig::TSequencing());
}

template <typename TSpec, typename TConfig, typename TSequencing>
void _openReads(Mapper<TSpec, TConfig> & mapper, TSequencing const & /*tag */)
{
    open(mapper.readsLoader, mapper.options.readsFile);
}

template <typename TSpec, typename TConfig>
void _openReads(Mapper<TSpec, TConfig> & mapper, SingleEnd const & /* tag */)
{
    open(mapper.readsLoader, mapper.options.readsFile.i1);
}

// ----------------------------------------------------------------------------
// Function loadReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void loadReads(Mapper<TSpec, TConfig> & mapper)
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

template <typename TSpec, typename TConfig, typename TReadSeqs, typename TReadSeqId, typename TDelegate>
inline void _selectSeeds(Mapper<TSpec, TConfig> & mapper, TReadSeqs const & readSeqs, TReadSeqId readSeqId, TDelegate & delegate)
{
    typedef typename StringSetPosition<TReadSeqs>::Type     TPos;
    typedef typename Value<TReadSeqs>::Type                 TReadSeq;
    typedef typename Size<TReadSeq>::Type                   TSize;

    TSize readLength = length(readSeqs[readSeqId]);
    TSize errorsPerRead = std::ceil(readLength * (mapper.options.errorRate / 100.0));
    TSize seedsPerRead = errorsPerRead + 1;
    TSize seedsLength = readLength / seedsPerRead;

    for (TSize seedId = 0; seedId < seedsPerRead; ++seedId)
        delegate(TPos(readSeqId, seedId * seedsLength), seedsLength);
}

template <typename TSpec, typename TConfig, typename TReadSeqs, typename TReadSeqId, typename TErrors, typename TDelegate>
inline void _selectSeeds(Mapper<TSpec, TConfig> & mapper, TReadSeqs const & readSeqs, TReadSeqId readSeqId, TErrors seedErrors, TDelegate & delegate)
{
    typedef typename StringSetPosition<TReadSeqs>::Type     TPos;
    typedef typename Value<TReadSeqs>::Type                 TReadSeq;
    typedef typename Size<TReadSeq>::Type                   TSize;

    TSize readLength = length(readSeqs[readSeqId]);
    TSize errorsPerRead = std::ceil(readLength * (mapper.options.errorRate / 100.0));
    TSize seedsPerRead = std::ceil((errorsPerRead + 1) / (seedErrors + 1.0));
    TSize seedsLength = readLength / seedsPerRead;

    for (TSize seedId = 0; seedId < seedsPerRead; ++seedId)
        delegate(TPos(readSeqId, seedId * seedsLength), seedsLength);
}

// ----------------------------------------------------------------------------
// Function selectSeeds()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSeeds, typename TReadSeqs, typename TReadsInfo, typename TErrors>
inline void selectSeeds(Mapper<TSpec, TConfig> & mapper, TSeeds & seeds, TReadSeqs & readSeqs, TReadsInfo & info, Pair<TErrors> seedErrors)
{
    typedef typename Value<TReadSeqs>::Type             TReadSeq;
    typedef typename Id<TReadSeqs>::Type                TId;
    typedef typename Size<TReadSeq>::Type               TSize;
    typedef typename Value<TSeeds>::Type                TSeedsSet;
    typedef SeedsCounter<TSize>                         TCounter;
    typedef SeedsManager<TSeedsSet, String<TSize> >     TManager;
    typedef Tuple<TCounter, 3>                          TCounters;
    typedef Tuple<TManager, 3>                          TManagers;

    TId readsCount = getReadSeqsCount(readSeqs);

    TCounters counters;
    for (TErrors errors = getValueI1(seedErrors); errors <= getValueI2(seedErrors); ++errors)
        resize(counters[errors], readsCount);

    for (TId readSeqId = 0; readSeqId < readsCount; ++readSeqId)
    {
        if (info[readSeqId].status == STATUS_UNSEEDED)
        {
            SEQAN_ASSERT_GEQ(info[readSeqId].seedErrors, getValueI1(seedErrors));
            SEQAN_ASSERT_LEQ(info[readSeqId].seedErrors, getValueI2(seedErrors));
            _selectSeeds(mapper, readSeqs, readSeqId, info[readSeqId].seedErrors, counters[info[readSeqId].seedErrors]);
        }
    }

    TManagers managers;
    for (TErrors errors = getValueI1(seedErrors); errors <= getValueI2(seedErrors); ++errors)
        init(managers[errors], seeds[errors], counters[errors].seedsPerRead);

    for (TId readSeqId = 0; readSeqId < readsCount; ++readSeqId)
    {
        if (info[readSeqId].status == STATUS_UNSEEDED)
        {
            SEQAN_ASSERT_GEQ(info[readSeqId].seedErrors, getValueI1(seedErrors));
            SEQAN_ASSERT_LEQ(info[readSeqId].seedErrors, getValueI2(seedErrors));
            _selectSeeds(mapper, readSeqs, readSeqId, info[readSeqId].seedErrors, managers[info[readSeqId].seedErrors]);
            info[readSeqId].status = STATUS_SEEDED;
        }
    }
}

// ----------------------------------------------------------------------------
// Function findSeeds()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename THitsString, typename TSeeder, typename TPattern>
inline void findSeeds(Mapper<TSpec, TConfig> & /* mapper */, THitsString & hits, TSeeder & seeder, TPattern & pattern)
{
    HitsManager<THitsString> manager(hits);

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif

    // Instantiate a pattern object.
//    TPattern pattern(seeds);

    // Initialize the delegate.
    init(manager, pattern);

    // Find hits.
    find(seeder, pattern, manager);

#ifdef PLATFORM_CUDA
    cudaPrintFreeMemory();
#endif
}

// ----------------------------------------------------------------------------
// Function writeHits()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs, typename THits, typename TSeedsSet, typename TReadsInfo, typename TFilename>
inline void writeHits(Mapper<TSpec, TConfig> & /* mapper */, TReadSeqs & readSeqs, THits & hits, TSeedsSet & seeds, TReadsInfo & info, TFilename const & filename)
{
    typedef Mapper<TSpec, TConfig>                      TMapper;
//    typedef typename TMapper::TSeedsSet                 TSeedsSet;
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

        if (info[readSeqId].status == STATUS_UNSEEDED) continue;
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
// Function reSeed()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs, typename THits, typename TSeeds, typename TReadsInfo>
inline void reSeed(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs, THits & hits, TSeeds & seeds, TReadsInfo & info)
{
    typedef Mapper<TSpec, TConfig>                      TMapper;
    typedef typename TMapper::TSeedsSet                 TSeedsSet;
    typedef typename Id<TSeedsSet>::Type                TSeedId;
    typedef Pair<TSeedId>                               TSeedIds;
    typedef typename TMapper::THit                      THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef Pair<THitId>                                THitIds;
    typedef typename Size<THit>::Type                   THitSize;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    TReadId readSeqsCount = getReadSeqsCount(readSeqs);

    for (TReadId readSeqId = 0; readSeqId < readSeqsCount; ++readSeqId)
    {
        // Count the hits per read.
        TSeedIds readSeedIds = getSeedIds(seeds, readSeqId);
        THitIds readHitIds = getHitIds(hits, readSeedIds);
        THitSize readHits = countHits<THitSize>(hits, readHitIds);

        // Re-seed hard reads.
        if (readHits > mapper.options.hitsThreshold)
        {
            // Guess a good seeding stragegy.
            info[readSeqId].seedErrors = (readHits < 200 * mapper.options.hitsThreshold) ? 1 : 2;
            info[readSeqId].status = STATUS_UNSEEDED;

            // Clear the hits of the read.
            clearHits(hits, readHitIds);
        }
    }
}

// ----------------------------------------------------------------------------
// Function selectAnchors()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs, typename THits, typename TSeeds, typename TReadsInfo>
inline void selectAnchors(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs, THits & hits, TSeeds & seeds, TReadsInfo & info)
{
    typedef Mapper<TSpec, TConfig>                      TMapper;
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

        TReadId anchorSeqId;
        TReadId otherSeqId;
        THitIds anchorHitIds;
        THitIds otherHitIds;

        // Choose the easiest read as the anchor.
        THitSize anchorHits = std::min(readHits, mateHits);

        if (anchorHits == readHits)
        {
            anchorSeqId = readSeqId;
            anchorHitIds = readHitIds;
            otherSeqId = mateSeqId;
            otherHitIds = mateHitIds;
        }
        else
        {
            anchorSeqId = mateSeqId;
            anchorHitIds = mateHitIds;
            otherSeqId = readSeqId;
            otherHitIds = readHitIds;
        }

        // Clear the hits of the other read.
        clearHits(hits, otherHitIds);
        info[otherSeqId].status = STATUS_UNMAPPABLE;

        // Re-seed hard anchors.
        if (anchorHits > mapper.options.hitsThreshold)
        {
            // Guess a good seeding stragegy.
            info[anchorSeqId].seedErrors = (anchorHits < 200 * mapper.options.hitsThreshold) ? 1 : 2;
            info[anchorSeqId].status = STATUS_UNSEEDED;

            // Clear the hits of the anchor.
            clearHits(hits, anchorHitIds);
        }
    }
}

// ----------------------------------------------------------------------------
// Function extendHits()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs, typename THits, typename TSeeds>
inline void extendHits(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs, THits & hits, TSeeds & seeds)
{
    typedef Mapper<TSpec, TConfig>                      TMapper;

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

template <typename TSpec, typename TConfig, typename TContigPos, typename TMatch>
inline void _getMateContigPos(Mapper<TSpec, TConfig> & mapper,
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
inline void _getMateContigPos(Mapper<TSpec, TConfig> & mapper,
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

template <typename TSpec, typename TConfig, typename TReadSeqs, typename TMatches>
inline void findMates(Mapper<TSpec, TConfig> & mapper,
                      TReadSeqs & readSeqs,
                      TMatches const & anchors,
                      TMatches & mates)
{
    typedef Mapper<TSpec, TConfig>                      TMapper;
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

template <typename TSpec, typename TConfig>
void mapReads(Mapper<TSpec, TConfig> & mapper)
{
//SEQAN_OMP_PRAGMA(critical(_mapper_mapReads_filter))
//{
    _mapReads(mapper, getSeqs(mapper.reads), typename TConfig::TAnchoring());
//}
}

// ----------------------------------------------------------------------------
// Function _mapReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
void _mapReads(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs, AnchorBoth)
{
    typedef Mapper<TSpec, TConfig>          TMapper;
    typedef typename TMapper::TSeedsExt     TSeedsExt;
    typedef typename TMapper::TSeedsApx     TSeedsApx;

    clear(mapper.hits[0]);
    clear(mapper.hits[1]);
    clear(mapper.hits[2]);

    clear(mapper.info);
    resize(mapper.info, getReadSeqsCount(readSeqs));

    start(mapper.timer);
    selectSeeds(mapper, mapper.seeds, readSeqs, mapper.info, Pair<unsigned>(0,0));
    std::cout << "Seeds count:\t\t\t" << length(mapper.seeds[0]) << std::endl;
    TSeedsExt seeds0(mapper.seeds[0]);
    findSeeds(mapper, mapper.hits[0], mapper.seederExt, seeds0);
    stop(mapper.timer);
    std::cout << "Seeding time:\t\t\t" << mapper.timer << std::endl;

//#ifdef _OPENMP
//    sortHits(mapper.hits[0]);
//#endif

    reSeed(mapper, readSeqs, mapper.hits[0], mapper.seeds[0], mapper.info);
    std::cout << "Hits count:\t\t\t" << countHits<unsigned long>(mapper.hits[0]) << std::endl;
//    writeHits(mapper, readSeqs, mapper.hits[0], mapper.seeds[0], mapper.info, "hits_0.csv");

    selectSeeds(mapper, mapper.seeds, readSeqs, mapper.info, Pair<unsigned>(1,2));

    start(mapper.timer);
    std::cout << "Seeds count:\t\t\t" << length(mapper.seeds[1]) << std::endl;
    TSeedsApx seeds1(mapper.seeds[1]);
    setScoreThreshold(mapper.seederApx, 1);
    findSeeds(mapper, mapper.hits[1], mapper.seederApx, seeds1);
    stop(mapper.timer);
    std::cout << "Seeding time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Hits count:\t\t\t" << countHits<unsigned long>(mapper.hits[1]) << std::endl;
//    writeHits(mapper, readSeqs, mapper.hits[1], mapper.seeds[1], mapper.info, "hits_1.csv");

    start(mapper.timer);
    std::cout << "Seeds count:\t\t\t" << length(mapper.seeds[2]) << std::endl;
    TSeedsApx seeds2(mapper.seeds[2]);
    setScoreThreshold(mapper.seederApx, 2);
    findSeeds(mapper, mapper.hits[2], mapper.seederApx, seeds2);
    stop(mapper.timer);
    std::cout << "Seeding time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Hits count:\t\t\t" << countHits<unsigned long>(mapper.hits[2]) << std::endl;
//    writeHits(mapper, readSeqs, mapper.hits[2], mapper.seeds[2], mapper.info, "hits_2.csv");

    clear(mapper.anchors);
    reserve(mapper.anchors, countHits<unsigned long>(mapper.hits[0]) + countHits<unsigned long>(mapper.hits[1]) + countHits<unsigned long>(mapper.hits[2]) / 5);

    start(mapper.timer);
    extendHits(mapper, readSeqs, mapper.hits[0], mapper.seeds[0]);
    extendHits(mapper, readSeqs, mapper.hits[1], mapper.seeds[1]);
    extendHits(mapper, readSeqs, mapper.hits[2], mapper.seeds[2]);
    stop(mapper.timer);
    std::cout << "Extension time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << length(mapper.anchors) << std::endl;

    start(mapper.timer);
    removeDuplicateMatches(mapper.anchors);
    stop(mapper.timer);
    std::cout << "Compaction time:\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << length(mapper.anchors) << std::endl;
    std::cout << "Anchored pairs:\t\t\t" << countMatches(readSeqs, mapper.anchors, typename TConfig::TSequencing()) << std::endl;
}

template <typename TSpec, typename TConfig, typename TReadSeqs>
void _mapReads(Mapper<TSpec, TConfig> & mapper, TReadSeqs & readSeqs, AnchorOne)
{
    typedef Mapper<TSpec, TConfig>          TMapper;
    typedef typename TMapper::TSeedsExt     TSeedsExt;
    typedef typename TMapper::TSeedsApx     TSeedsApx;

    clear(mapper.hits[0]);
    clear(mapper.hits[1]);
    clear(mapper.hits[2]);

    clear(mapper.info);
    resize(mapper.info, getReadSeqsCount(readSeqs));

    start(mapper.timer);
    selectSeeds(mapper, mapper.seeds, readSeqs, mapper.info, Pair<unsigned>(0,0));
    std::cout << "Seeds count:\t\t\t" << length(mapper.seeds[0]) << std::endl;
    TSeedsExt seeds0(mapper.seeds[0]);
    findSeeds(mapper, mapper.hits[0], mapper.seederExt, seeds0);
    stop(mapper.timer);
    std::cout << "Seeding time:\t\t\t" << mapper.timer << std::endl;

//#ifdef _OPENMP
//    sortHits(mapper.hits[0]);
//#endif

    selectAnchors(mapper, readSeqs, mapper.hits[0], mapper.seeds[0], mapper.info);
    selectSeeds(mapper, mapper.seeds, readSeqs, mapper.info, Pair<unsigned>(1,2));
    std::cout << "Hits count:\t\t\t" << countHits<unsigned long>(mapper.hits[0]) << std::endl;

    start(mapper.timer);
    std::cout << "Seeds count:\t\t\t" << length(mapper.seeds[1]) << std::endl;
    TSeedsApx seeds1(mapper.seeds[1]);
    setScoreThreshold(mapper.seederApx, 1);
    findSeeds(mapper, mapper.hits[1], mapper.seederApx, seeds1);
    stop(mapper.timer);
    std::cout << "Seeding time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Hits count:\t\t\t" << countHits<unsigned long>(mapper.hits[1]) << std::endl;

    start(mapper.timer);
    std::cout << "Seeds count:\t\t\t" << length(mapper.seeds[2]) << std::endl;
    TSeedsApx seeds2(mapper.seeds[2]);
    setScoreThreshold(mapper.seederApx, 2);
    findSeeds(mapper, mapper.hits[2], mapper.seederApx, seeds2);
    stop(mapper.timer);
    std::cout << "Seeding time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Hits count:\t\t\t" << countHits<unsigned long>(mapper.hits[2]) << std::endl;

    clear(mapper.anchors);
    reserve(mapper.anchors, countHits<unsigned long>(mapper.hits[0]) + countHits<unsigned long>(mapper.hits[1]) + countHits<unsigned long>(mapper.hits[2]) / 5);

    start(mapper.timer);
    extendHits(mapper, readSeqs, mapper.hits[0], mapper.seeds[0]);
    extendHits(mapper, readSeqs, mapper.hits[1], mapper.seeds[1]);
    extendHits(mapper, readSeqs, mapper.hits[2], mapper.seeds[2]);
    stop(mapper.timer);
    std::cout << "Extension time:\t\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << length(mapper.anchors) << std::endl;

    start(mapper.timer);
    removeDuplicateMatches(mapper.anchors);
    stop(mapper.timer);
    std::cout << "Compaction time:\t\t" << mapper.timer << std::endl;
    std::cout << "Anchors count:\t\t\t" << length(mapper.anchors) << std::endl;
    std::cout << "Anchored pairs:\t\t\t" << countMatches(readSeqs, mapper.anchors, typename TConfig::TSequencing()) << std::endl;

    start(mapper.timer);
    clear(mapper.mates);
    reserve(mapper.mates, length(mapper.anchors), Exact());
    findMates(mapper, readSeqs, mapper.anchors, mapper.mates);
    stop(mapper.timer);
    std::cout << "Verification time:\t\t" << mapper.timer << std::endl;
    std::cout << "Mates count:\t\t\t" << length(mapper.mates) << std::endl;
    std::cout << "Mapped pairs:\t\t\t" << countMatches(readSeqs, mapper.mates, typename TConfig::TSequencing()) << std::endl;

//    start(mapper.timer);
//    runWriter(mapper.writer, readSeqs);
//    stop(mapper.timer);
//    std::cout << "Writing time:\t\t\t" << mapper.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void runMapper(Mapper<TSpec, TConfig> & mapper)
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

template <typename TSpec, typename TConfig>
void runMapper(Mapper<TSpec, TConfig> & mapper, Parallel)
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

// ----------------------------------------------------------------------------
// Function spawnMapper()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TSequencing, typename TStrategy, typename TAnchoring>
void spawnMapper(Options const & options,
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
