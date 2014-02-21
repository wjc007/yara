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
    typedef typename Traits::TContigs          TContigs;
    typedef typename Traits::TReads            TReads;
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
    TContigs const &        contigs;
    TReads const &          reads;
    TMatchesSet const &     matchesSet;
    Options const &         options;

    MatchesWriter(TOutputStream & outputStream,
                  TOutputContext & outputCtx,
                  TReadsContext const & ctx,
                  TContigs const & contigs,
                  TReads const & reads,
                  TMatchesSet const & matchesSet,
                  Options const & options) :
        outputStream(outputStream),
        outputCtx(outputCtx),
        ctx(ctx),
        contigs(contigs),
        reads(reads),
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

// ----------------------------------------------------------------------------
// Class QualityExtractor
// ----------------------------------------------------------------------------
// TODO(esiragusa): remove this when new tokenization gets into develop.

template <typename TValue>
struct QualityExtractor : public std::unary_function<TValue, char>
{
    inline char operator()(TValue const & x) const
    {
        return '!' + static_cast<char>(getQualityValue(x));
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setQual()
// ----------------------------------------------------------------------------

template <typename TString>
inline void setQual(BamAlignmentRecord & record, TString const & string)
{
    typedef typename Value<TString>::Type                               TAlphabet;
    typedef QualityExtractor<TAlphabet>                                 TQualityExtractor;
    typedef ModifiedString<TString const, ModView<TQualityExtractor> >  TQualities;

    TQualities qual(string);
    record.qual = qual;
}

// ----------------------------------------------------------------------------
// Function _writeMatchesImpl()
// ----------------------------------------------------------------------------
// Writes one block of matches.

template <typename TSpec, typename Traits, typename TMatches>
inline void _writeMatchesImpl(MatchesWriter<TSpec, Traits> & me, TMatches const & matches)
{
    typedef typename Value<TMatches>::Type          TMatch;
    typedef typename Size<TMatches>::Type           TSize;

    // The first match is supposed to be the best one.
    TMatch const & primary = front(matches);
    TSize bestCount = countBestMatches(matches);

    // Clear the record.
    clear(me.record);

    // Set primary alignment information.
    me.record.qName = me.reads.names[getReadId(primary)];
    me.record.seq = me.reads.seqs[getReadSeqId(primary, me.reads.seqs)];
    setQual(me.record, me.reads.seqs[getReadSeqId(primary, me.reads.seqs)]);

    // Set orientation.
    if (onReverseStrand(primary))
        me.record.flag |= BAM_FLAG_RC;

    // Set contig position.
    me.record.rID = getContigId(primary);
    me.record.beginPos = getContigBegin(primary);

    // Set mapping quality.
    me.record.mapQ = getScore(primary);

    // Set alignment.
//    setAlignment(record, me.contigs, primary, primary, alignFunctor);

    // Set number of errors.
    appendTagValue(me.record.tags, "NM", getErrors(primary));

    // Set number of cooptimal and suboptimal hits.
    appendTagValue(me.record.tags, "X0", bestCount, 'i');
    appendTagValue(me.record.tags, "X1", length(matches) - bestCount, 'i');

    // Set type as unique or repeat.
    char type = (bestCount == 1) ? 'U' : 'R';
    appendTagValue(me.record.tags, "XT", type, 'A');

//    appendTagValue(me.record.tags, "XA", "chr,pos,CIGAR,NM;...", 'Z');

    // Write record to output stream.
    write2(me.outputStream, me.record, me.outputCtx, typename Traits::TOutputFormat());
}

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_WRITER_H_
