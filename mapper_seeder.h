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

#ifndef APP_CUDAMAPPER_MAPPER_SEEDER_H_
#define APP_CUDAMAPPER_MAPPER_SEEDER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadsSeeder
// ----------------------------------------------------------------------------
// One instance per thread.

template <typename TSpec, typename Traits>
struct ReadsSeeder
{
    typedef typename Traits::TReadsContext     TReadsContext;
    typedef typename Traits::TSeeds            TSeeds;
    typedef typename Traits::TReadSeqs         TReadSeqs;

    // Shared-memory read-write data.
    TReadsContext &     ctx;
    TSeeds &            seeds;

    // Shared-memory read-only data.
    TReadSeqs const &   readSeqs;
    Options const &     options;

    ReadsSeeder(TReadsContext & ctx, TSeeds & seeds, TReadSeqs const & readSeqs, Options const & options) :
        ctx(ctx),
        seeds(seeds),
        readSeqs(readSeqs),
        options(options)
    {
        // Iterate over all reads.
        iterate(readSeqs, *this, Rooted(), Parallel());
    }

    template <typename TReadSeqsIterator>
    void operator() (TReadSeqsIterator const & it)
    {
        _seedReadImpl(*this, it);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _seedReadImpl()
// ----------------------------------------------------------------------------
// Samples seeds for one read.

template <typename TSpec, typename Traits, typename TErrors>
inline void _seedReadImpl(ReadsSeeder<TSpec, Traits> & me, Pair<TErrors> seedErrors)
{
    typedef typename Traits::TReadSeqs                  TReadSeqs;
    typedef typename Value<TReadSeqs>::Type             TReadSeq;
    typedef typename Id<TReadSeqs>::Type                TId;
    typedef typename Size<TReadSeq>::Type               TSize;
    typedef typename Traits::TSeeds                     TSeeds;
    typedef SeedsCounter<TSize>                         TCounter;
    typedef SeedsManager<TSeeds, String<TSize> >        TManager;
    typedef Tuple<TCounter, Traits::BUCKETS>            TCounters;
    typedef Tuple<TManager, Traits::BUCKETS>            TManagers;

    TId readsCount = getReadSeqsCount(me.readSeqs);

    // Initialize counters.
    TCounters counters;
    for (TErrors errors = getValueI1(seedErrors); errors <= getValueI2(seedErrors); ++errors)
        resize(counters[errors], readsCount);

    // Count seeds.
    // Counters(ctx)
//    iterate(readSeqs, counters, Rooted(), Parallel());

    for (TId readSeqId = 0; readSeqId < readsCount; ++readSeqId)
    {
        if (getStatus(me.ctx, readSeqId) == STATUS_UNSEEDED)
        {
            unsigned char errors = getSeedErrors(me.ctx, readSeqId);
            SEQAN_ASSERT_GEQ(errors, getValueI1(seedErrors));
            SEQAN_ASSERT_LEQ(errors, getValueI2(seedErrors));
            _getSeedsPerRead(me, readSeqId, counters[errors]);
        }
    }

    // Initialize managers.
    TManagers managers;
    for (TErrors errors = getValueI1(seedErrors); errors <= getValueI2(seedErrors); ++errors)
        init(managers[errors], me.seeds[errors], counters[errors].seedsPerRead);

    // Select seeds.
//    iterate(readSeqs, managers, Rooted(), Parallel());

    for (TId readSeqId = 0; readSeqId < readsCount; ++readSeqId)
    {
        if (getStatus(me.ctx, readSeqId) == STATUS_UNSEEDED)
        {
            unsigned char errors = getSeedErrors(me.ctx, readSeqId);
            SEQAN_ASSERT_GEQ(errors, getValueI1(seedErrors));
            SEQAN_ASSERT_LEQ(errors, getValueI2(seedErrors));
            _getSeedsPerRead(me, readSeqId, managers[errors]);
            setStatus(me.ctx, readSeqId, STATUS_SEEDED);
        }
    }
}

// ----------------------------------------------------------------------------
// Function _getSeedsPerRead()
// ----------------------------------------------------------------------------
// Enumerates the seeds for a given read sequence.

template <typename TSpec, typename Traits, typename TReadSeqId, typename TDelegate>
inline void _getSeedsPerRead(ReadsSeeder<TSpec, Traits> & me, TReadSeqId readSeqId, TDelegate & delegate)
{
    typedef typename Traits::TReadSeqs                      TReadSeqs;
    typedef typename StringSetPosition<TReadSeqs>::Type     TPos;
    typedef typename Value<TReadSeqs>::Type                 TReadSeq;
    typedef typename Size<TReadSeq>::Type                   TSize;

    TSize readLength = length(me.readSeqs[readSeqId]);
    TSize readErrors = getReadErrors(me.options, readLength);
    TSize seedErrors = getSeedErrors(me.ctx, readSeqId);
    TSize seedsCount = std::ceil((readErrors + 1) / (seedErrors + 1.0));
    TSize seedsLength = readLength / seedsCount;

    for (TSize seedId = 0; seedId < seedsCount; ++seedId)
        delegate(TPos(readSeqId, seedId * seedsLength), seedsLength);
}

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_SEEDER_H_
