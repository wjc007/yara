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

#ifndef APP_CUDAMAPPER_MAPPER_WRITER_H_
#define APP_CUDAMAPPER_MAPPER_WRITER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class MatchesWriter
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
struct MatchesWriter
{
    typedef typename Traits::TStore            TStore;
    typedef typename Traits::TMatchesSet       TMatchesSet;
    typedef typename Traits::TOutputStream     TOutputStream;
    typedef typename Traits::TOutputContext    TOutputContext;
    typedef typename Traits::TReadsContext     TReadsContext;

    // Thread-private data.
    BamAlignmentRecord      record;

    // Shared-memory read-write data.
    TOutputStream &         outputStream;
    TOutputContext &        outputCtx;

    // Shared-memory read-only data.
    TReadsContext const &   ctx;
    TStore const &          store;
    TMatchesSet const &     matchesSet;
    Options const &         options;

    MatchesWriter(TOutputStream & outputStream,
                  TOutputContext & outputCtx,
                  TReadsContext const & ctx,
                  TStore const & store,
                  TMatchesSet const & matchesSet,
                  Options const & options) :
        outputStream(outputStream),
        outputCtx(outputCtx),
        ctx(ctx),
        store(store),
        matchesSet(matchesSet),
        options(options)
    {
        // Process all matches.
        // TODO(esiragusa): insure that forEach() does not copy the functor.
        forEach(matchesSet, *this, Serial());
    }

    template <typename TMatches>
    void operator() (TMatches const & matches)
    {
        _writeMatchesImpl(*this, matches);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _writeMatchesImpl()
// ----------------------------------------------------------------------------
// Writes one block of matches.

template <typename TSpec, typename Traits, typename TMatches>
inline void _writeMatchesImpl(MatchesWriter<TSpec, Traits> & me, TMatches const & matches)
{
    typedef typename Value<TMatches>::Type          TMatch;

    // The first match is supposed to be the best one.
    TMatch const & primary = front(matches);

    clear(me.record.tags);
    me.record.flag = 0;

    // Add primary alignment information.
    setName(me.record, me.store, primary);
//    setSeqAndQual(me.record, me.store, primary);
    setOrientation(me.record, me.store, primary);
//    setPosition(me.record, me.store, primary);
//    me.record.rID = getContigId(primary);
    me.record.beginPos = getContigBegin(primary);

//    setAlignment(record, store, primary, primary, alignFunctor);
    setScore(me.record, me.store, primary);

    // Clear mate information.
//    clearMateInfo(me.record, me.store, primary);
    clearMatePosition(me.record, me.store);

    // Add secondary match information.
//    addSecondaryMatch(me.record, me.store, itBegin + 1, itEnd);

    // Write record to output stream.
    write2(me.outputStream, me.record, me.outputCtx, typename Traits::TOutputFormat());
}

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_WRITER_H_
