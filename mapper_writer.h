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
    typedef typename Traits::TPairsSet         TPairsSet;
    typedef typename Traits::TOutputStream     TOutputStream;
    typedef typename Traits::TOutputContext    TOutputContext;
    typedef typename Traits::TReadsContext     TReadsContext;

    // Thread-private data.
    BamAlignmentRecord      record;
    CharString              xa;

    // Shared-memory read-write data.
    TOutputStream &         outputStream;
    TOutputContext &        outputCtx;
    TMatchesSet &           matchesSet;
    TPairsSet &             pairsSet;

    // Shared-memory read-only data.
    TReadsContext const &   ctx;
    TContigs const &        contigs;
    TReads const &          reads;
    Options const &         options;

    MatchesWriter(TOutputStream & outputStream,
                  TOutputContext & outputCtx,
                  TMatchesSet & matchesSet,
                  TPairsSet & pairsSet,
                  TReadsContext const & ctx,
                  TContigs const & contigs,
                  TReads const & reads,
                  Options const & options) :
        outputStream(outputStream),
        outputCtx(outputCtx),
        matchesSet(matchesSet),
        pairsSet(pairsSet),
        ctx(ctx),
        contigs(contigs),
        reads(reads),
        options(options)
    {
        _writeAllMatchesImpl(*this, typename Traits::TSequencing());
    }

    template <typename TIterator>
    void operator() (TIterator const & it)
    {
        _writeMatchesImpl(*this, it, typename Traits::TSequencing());
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
// Function fillHeader()
// ----------------------------------------------------------------------------

template <typename TContigSeqs, typename TContigNames>
inline void fillHeader(BamHeader & header, TContigSeqs const & seqs, TContigNames const & names)
{
    typedef typename Iterator<TContigSeqs const, Standard>::Type    TContigSeqsIter;
    typedef typename Iterator<TContigNames const, Standard>::Type   TContigNamesIter;

    typedef BamHeader::TSequenceInfo                                TSequenceInfo;
    typedef BamHeaderRecord::TTag                                   TTag;

    // Fill first header line.
    BamHeaderRecord firstRecord;
    firstRecord.type = BAM_HEADER_FIRST;
    appendValue(firstRecord.tags, TTag("VN", "1.4"));
    appendValue(firstRecord.tags, TTag("SO", "queryname"));
    appendValue(header.records, firstRecord);

    // Fill sequence info header line.
    TContigSeqsIter sit = begin(seqs, Standard());
    TContigSeqsIter sitEnd = end(seqs, Standard());
    TContigNamesIter nit = begin(names, Standard());
    TContigNamesIter nitEnd = end(names, Standard());

    for (; sit != sitEnd && nit != nitEnd; ++sit, ++nit)
        appendValue(header.sequenceInfos, TSequenceInfo(*nit, length(*sit)));

    // Fill program header line.
    BamHeaderRecord pgRecord;
    pgRecord.type = BAM_HEADER_PROGRAM;
    appendValue(pgRecord.tags, TTag("ID", "SeqAn"));
    appendValue(header.records, pgRecord);
}

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
// Function append*()
// ----------------------------------------------------------------------------

template <typename TErrors>
inline void appendErrors(BamAlignmentRecord & record, TErrors errors)
{
    appendTagValue(record.tags, "NM", errors, 'i');
}

template <typename TCount>
inline void appendCooptimalCount(BamAlignmentRecord & record, TCount count)
{
    appendTagValue(record.tags, "X0", count, 'i');
}

template <typename TCount>
inline void appendSuboptimalCount(BamAlignmentRecord & record, TCount count)
{
    appendTagValue(record.tags, "X1", count, 'i');
}

inline void appendType(BamAlignmentRecord & record, bool unique)
{
    appendTagValue(record.tags, "XT", unique ? 'U' : 'R', 'A');
}

template <typename TString>
inline void appendAlignments(BamAlignmentRecord & record, TString const & xa)
{
    // XA:Z:(chr,pos,CIGAR,NM;)+
    appendTagValue(record.tags, "XA", xa, 'Z');
}

template <typename TSpec, typename Traits, typename TMatches>
inline void _fillXa(MatchesWriter<TSpec, Traits> & me, TMatches const & matches)
{
    typedef typename Iterator<TMatches const, Standard>::Type   TIter;

    clear(me.xa);
    TIter itEnd = end(matches, Standard());
    for (TIter it = begin(matches, Standard()) + 1; it != itEnd; ++it)
    {
        append(me.xa, nameStore(me.outputCtx)[getContigId(value(it))]);
        appendValue(me.xa, ',');
//        append(me.xa, getContigBegin(value(it)) + 1);
        appendValue(me.xa, '1');
        appendValue(me.xa, ',');
        appendValue(me.xa, '*');
        appendValue(me.xa, ',');
        appendValue(me.xa, '0' + getErrors(value(it)));
        appendValue(me.xa, ';');
    }
}

// ----------------------------------------------------------------------------
// Function _writeAllMatchesImpl(SingleEnd)
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
inline void _writeAllMatchesImpl(MatchesWriter<TSpec, Traits> & me, SingleEnd)
{
    typedef typename Traits::TMatchesSet                    TMatchesSet;
    typedef typename Iterator<TMatchesSet, Standard>::Type  TMatchesIt;
//    typedef typename Value<TMatchesSet>::Type             TMatches;

    // Sort each set of matches by errors.
    iterate(me.matchesSet, sortMatches<TMatchesIt, SortErrors>, Standard(), typename Traits::TThreading());
//    forEach(me.matchesSet, sortMatches<TMatches, SortErrors>, typename TConfig::TThreading());

    // Write all matches.
    iterate(me.matchesSet, me, Standard(), Serial());
}

// ----------------------------------------------------------------------------
// Function _writeAllMatchesImpl(PairedEnd)
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
inline void _writeAllMatchesImpl(MatchesWriter<TSpec, Traits> & me, PairedEnd)
{
    typedef typename Traits::TMatchesSet                    TMatchesSet;
    typedef typename Iterator<TMatchesSet, Standard>::Type  TMatchesIt;
    typedef typename Traits::TPairsSet                      TPairsSet;
    typedef typename Iterator<TPairsSet, Standard>::Type    TPairsIt;

    // Sort each set of pairs by errors and insert size deviation.
//    iterate(me.pairsSet, sortPairs<TPairsIt, SortErrors>, Standard(), typename Traits::TThreading());

    // Sort each set of unpaired matches by errors.
    iterate(me.matchesSet, sortMatches<TMatchesIt, SortErrors>, Standard(), typename Traits::TThreading());

    // Write all pairs.
    iterate(me.pairsSet, me, Standard(), Serial());
}

// ----------------------------------------------------------------------------
// Function _writeMatchesImpl(SingleEnd)
// ----------------------------------------------------------------------------
// Writes one block of matches.

template <typename TSpec, typename Traits, typename TMatchesIt>
inline void _writeMatchesImpl(MatchesWriter<TSpec, Traits> & me, TMatchesIt const & it, SingleEnd)
{
    typedef typename Value<TMatchesIt const>::Type  TMatches;
    typedef typename Value<TMatches>::Type          TMatch;

    TMatches const & matches = value(it);

    if (empty(matches))
        _writeUnmappedReadImpl(me, position(it, me.matchesSet));
    else
        _writeMappedReadImpl(me, matches);
}

// ----------------------------------------------------------------------------
// Function _writeMatchesImpl(PairedEnd)
// ----------------------------------------------------------------------------
// Writes one block of pairs.

template <typename TSpec, typename Traits, typename TPairsIt>
inline void _writeMatchesImpl(MatchesWriter<TSpec, Traits> & me, TPairsIt const & it, PairedEnd)
{
    typedef typename Traits::TMatchesSet                    TMatchesSet;
    typedef typename Iterator<TMatchesSet, Standard>::Type  TMatchesIt;
    typedef typename Value<TPairsIt const>::Type            TPairs;
    typedef typename Value<TPairs>::Type                    TPair;
    typedef typename Traits::TReadSeqs                      TReadSeqs;
    typedef typename Size<TReadSeqs>::Type                  TReadId;

    TPairs const & pairs = value(it);

    if (empty(pairs))
    {
        TReadId firstId = getFirstMateFwdSeqId(me.reads.seqs, position(it, me.pairsSet));
        TReadId secondId = getSecondMateFwdSeqId(me.reads.seqs, position(it, me.pairsSet));

        _writeMatchesImpl(me, begin(me.matchesSet, Standard()) + firstId, SingleEnd());
        _writeMatchesImpl(me, begin(me.matchesSet, Standard()) + secondId, SingleEnd());
    }
    else
    {
        // Write paired reads.
//        _writePairedReadsImpl(me, pairs);
    }
}

// ----------------------------------------------------------------------------
// Function _writeUnmappedReadImpl()
// ----------------------------------------------------------------------------
// Writes one unmapped read.

template <typename TSpec, typename Traits, typename TReadId>
inline void _writeUnmappedReadImpl(MatchesWriter<TSpec, Traits> & me, TReadId readId)
{
    // Clear the record.
    clear(me.record);

    // Set primary alignment information.
    me.record.qName = me.reads.names[readId];
    me.record.seq = me.reads.seqs[readId];
    setQual(me.record, me.reads.seqs[readId]);

    // Set read as unmapped.
    me.record.flag |= BAM_FLAG_UNMAPPED;

    // Write record to output stream.
    write2(me.outputStream, me.record, me.outputCtx, typename Traits::TOutputFormat());
}

// ----------------------------------------------------------------------------
// Function _writeMappedReadImpl()
// ----------------------------------------------------------------------------
// Writes one block of matches.

template <typename TSpec, typename Traits, typename TMatches>
inline void _writeMappedReadImpl(MatchesWriter<TSpec, Traits> & me, TMatches const & matches)
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
    appendErrors(me.record, getErrors(primary));

    // Set number of secondary alignments.
//    appendTagValue(me.record.tags, "NH", 1, 'i');

    // Set hit index.
//    appendTagValue(me.record.tags, "HI", 1, 'i');

    // Set number of cooptimal and suboptimal hits.
    appendCooptimalCount(me.record, bestCount);
    appendSuboptimalCount(me.record, length(matches) - bestCount);

    // Set type as unique or repeat.
    appendType(me.record, bestCount == 1);

    // Append secondary matches.
    if (length(matches) > 1)
    {
        _fillXa(me, matches);
        appendAlignments(me.record, me.xa);
    }

    // Write record to output stream.
    write2(me.outputStream, me.record, me.outputCtx, typename Traits::TOutputFormat());
}

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_WRITER_H_
