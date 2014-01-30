// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
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

#ifndef APP_CUDAMAPPER_CLASSIFIER_H_
#define APP_CUDAMAPPER_CLASSIFIER_H_

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

struct Options;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadsClassifier
// ----------------------------------------------------------------------------
// One instance per thread.

template <typename TReadsContext, typename THits, typename TSeeds, typename TConfig>
struct ReadsClassifier
{
    // Shared-memory read-write data.
    TReadsContext &     ctx;
    THits &             hits;

    // Shared-memory read-only data.
    TSeeds const &      seeds;
    // TODO(esiragusa): remove TReadSeqs - it is a global typedef
    TReadSeqs const &   readSeqs;
    Options const &     options;

    ReadsClassifier(TReadsContext & ctx, THits & hits, TSeeds const & seeds, Options const & options) :
        ctx(ctx),
        hits(hits),
        seeds(seeds),
        readSeqs(host(seeds)),
        options(options)
    {
        _classifyReadsImpl(*this, typename TConfig::TAnchoring());
    }

    // NOTE(esiragusa): This is called on firstprivate.
    ReadsClassifier(ReadsClassifier const & other) :
        ctx(other.ctx),
        hits(other.hits),
        seeds(other.seeds),
        readSeqs(other.readSeqs),
        options(other.options)
    {}

    template <typename TReadSeqsIterator>
    void operator() (TReadSeqsIterator const & it)
    {
        _classifyReadImpl(*this, it, typename TConfig::TAnchoring());
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _classifyReadsImpl(); Default
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename THits, typename TSeeds, typename TConfig, typename TAnchoring>
inline void _classifyReadsImpl(ReadsClassifier<TReadsContext, THits, TSeeds, TConfig> & classifier, TAnchoring)
{
    // Iterate over all reads.
    iterate(classifier.readSeqs, classifier, Rooted(), Parallel());
}

// ----------------------------------------------------------------------------
// Function _classifyReadsImpl(); AnchorOne
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename THits, typename TSeeds, typename TConfig>
inline void _classifyReadsImpl(ReadsClassifier<TReadsContext, THits, TSeeds, TConfig> & classifier, AnchorOne)
{
    // TODO(esiragusa): fix this, it returns an empty prefix!
    // Iterate over all pairs.
    iterate(prefix(classifier.readSeqs, getReadsCount(classifier.readSeqs)), classifier, Rooted(), Parallel());
}

// ----------------------------------------------------------------------------
// Function _classifyReadImpl(); Default
// ----------------------------------------------------------------------------
// Raises the seeds errors, mark for reseeding and clears the hits of hard reads.

template <typename TReadsContext, typename THits, typename TSeeds, typename TConfig,
          typename TReadSeqsIterator, typename TAnchoring>
inline void _classifyReadImpl(ReadsClassifier<TReadsContext, THits, TSeeds, TConfig> & classifier,
                              TReadSeqsIterator const & it, TAnchoring)
{
    typedef typename Value<THits>::Type                 THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef Pair<THitId>                                THitIds;
    typedef typename Size<THit>::Type                   THitSize;
    typedef typename Id<TSeeds>::Type                   TSeedId;
    typedef Pair<TSeedId>                               TSeedIds;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    TReadId readSeqId = position(it);

    // Count the hits per read.
    TSeedIds readSeedIds = getSeedIds(classifier.seeds, readSeqId);
    THitIds readHitIds = getHitIds(classifier.hits, readSeedIds);
    THitSize readHits = countHits<THitSize>(classifier.hits, readHitIds);

    // Re-seed hard reads.
    if (readHits > classifier.options.hitsThreshold)
    {
        // Guess a good seeding stragegy.
        setSeedErrors(classifier.ctx, readSeqId, (readHits < 200 * classifier.options.hitsThreshold) ? 1 : 2);
        setStatus(classifier.ctx, readSeqId, STATUS_UNSEEDED);

        // Clear the hits of the read.
        clearHits(classifier.hits, readHitIds);
    }
}

// ----------------------------------------------------------------------------
// Function _classifyReadImpl(); AnchorOne
// ----------------------------------------------------------------------------
// Selects the mate to anchor; raises the seeds errors, mark for reseeding and clears the hits of hard anchors.

template <typename TReadsContext, typename THits, typename TSeeds, typename TConfig, typename TReadSeqsIterator>
inline void _classifyReadImpl(ReadsClassifier<TReadsContext, THits, TSeeds, TConfig> & classifier, TReadSeqsIterator const & it, AnchorOne)
{
    typedef typename Value<THits>::Type                 THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef Pair<THitId>                                THitIds;
    typedef typename Size<THit>::Type                   THitSize;
    typedef typename Id<TSeeds>::Type                   TSeedId;
    typedef Pair<TSeedId>                               TSeedIds;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    // Get readSeqId.
    TReadId readSeqId = position(it);

    // Get mate id.
    TReadId mateSeqId = getMateSeqId(classifier.readSeqs, readSeqId);

    // Get seed ids.
    TSeedIds readSeedIds = getSeedIds(classifier.seeds, readSeqId);
    TSeedIds mateSeedIds = getSeedIds(classifier.seeds, mateSeqId);

    // Get hit ids.
    THitIds readHitIds = getHitIds(classifier.hits, readSeedIds);
    THitIds mateHitIds = getHitIds(classifier.hits, mateSeedIds);

    // Count the hits of each read.
    THitSize readHits = countHits<THitSize>(classifier.hits, readHitIds);
    THitSize mateHits = countHits<THitSize>(classifier.hits, mateHitIds);

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
    clearHits(classifier.hits, otherHitIds);
    setStatus(classifier.ctx, otherSeqId, STATUS_UNMAPPABLE);

    // Re-seed hard anchors.
    if (anchorHits > classifier.options.hitsThreshold)
    {
        // Guess a good seeding stragegy.
        setSeedErrors(classifier.ctx, anchorSeqId, (anchorHits < 200 * classifier.options.hitsThreshold) ? 1 : 2);
        setStatus(classifier.ctx, anchorSeqId, STATUS_UNSEEDED);

        // Clear the hits of the anchor.
        clearHits(classifier.hits, anchorHitIds);
    }
}

#endif  // #ifndef APP_CUDAMAPPER_CLASSIFIER_H_
