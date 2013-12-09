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
// This file contains the Genome class.
// ==========================================================================

#ifndef APP_CUDAMAPPER_GENOME_H_
#define APP_CUDAMAPPER_GENOME_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/random.h>

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

template <typename TObject>
struct Contigs;

template <typename TObject>
struct GenomeHost;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class GenomeConfig
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct GenomeConfig
{
    typedef String<Dna5>            TContigSeq;
    typedef Owner<ConcatDirect<> >  TContigSpec;
};

// ----------------------------------------------------------------------------
// Class Genome
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = GenomeConfig<TSpec> >
struct Genome
{
    Holder<typename Contigs<Genome>::Type>  _contigs;

    Genome() {}

    template <typename TContigs>
    Genome(TContigs & contigs) :
        _contigs(contigs)
    {}
};

// ----------------------------------------------------------------------------
// Class GenomeLoader
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = GenomeConfig<> >
struct GenomeLoader
{
    typedef std::fstream                            TStream;
    typedef RecordReader<TStream, SinglePass<> >    TRecordReader;
    typedef typename GenomeHost<GenomeLoader>::Type TGenome;

    TStream                         _file;
    unsigned long                   _fileSize;
    AutoSeqStreamFormat             _fileFormat;
    std::auto_ptr<TRecordReader>    _reader;
    Rng<MersenneTwister>            _rng;
    Holder<TGenome>                 genome;

    GenomeLoader(TGenome & genome) :
        _fileSize(0),
        _rng(0xDEADBEEF),
        genome(genome)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GenomeHost<T>::Type                                   [TObject]
// ----------------------------------------------------------------------------

template <typename TObject>
struct GenomeHost {};

template <typename TObject>
struct GenomeHost<TObject const>
{
    typedef typename GenomeHost<TObject>::Type const    Type;
};

// ----------------------------------------------------------------------------
// Metafunction GenomeHost<T>::Type                              [GenomeLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
struct GenomeHost<GenomeLoader<TSpec, TConfig> >
{
    typedef Genome<TSpec, TConfig>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Contigs                                               [TObject]
// ----------------------------------------------------------------------------

template <typename TObject>
struct Contigs {};

template <typename TObject>
struct Contigs<TObject const>
{
    typedef typename Contigs<TObject>::Type const   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Contig                                                [TObject]
// ----------------------------------------------------------------------------

template <typename TObject>
struct Contig {};

template <typename TObject>
struct Contig<TObject const>
{
    typedef typename Contig<TObject>::Type const   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Contigs                                                [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
struct Contigs<Genome<TSpec, TConfig> >
{
    typedef StringSet<typename TConfig::TContigSeq, typename TConfig::TContigSpec>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Contig                                                 [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
struct Contig<Genome<TSpec, TConfig> > :
    Value<Contigs<Genome<TSpec, TConfig> > > {};

// ----------------------------------------------------------------------------
// Metafunction Size                                                   [Genome]
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSpec, typename TConfig>
struct Size<Genome<TSpec, TConfig> > :
    Size<typename Contigs<Genome<TSpec, TConfig> >::Type> {};
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setGenome()                                               [TObject]
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec, typename TConfig>
inline void
setGenome(TObject & object, Genome<TSpec, TConfig> const & genome)
{
    setValue(object.genome, genome);
}

// ----------------------------------------------------------------------------
// Function getGenome()                                               [TObject]
// ----------------------------------------------------------------------------

template <typename TObject>
inline typename GenomeHost<TObject>::Type &
getGenome(TObject & object)
{
    return value(object.genome);
}

// ----------------------------------------------------------------------------
// Function contigs()                                                  [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline typename Contigs<Genome<TSpec, TConfig> >::Type &
contigs(Genome<TSpec, TConfig> & genome)
{
    return value(genome._contigs);
}

template <typename TSpec, typename TConfig>
inline typename Contigs<Genome<TSpec, TConfig> const>::Type &
contigs(Genome<TSpec, TConfig> const & genome)
{
    return value(genome._contigs);
}

// ----------------------------------------------------------------------------
// Function contigs()                                                  [Genome]
// ----------------------------------------------------------------------------
//
//template <typename TSpec, typename TConfig, typename TContigId>
//inline typename Contigs<Genome<TSpec, TConfig> >::Type &
//contig(Genome<TSpec, TConfig> & genome, TContigId contigId)
//{
//    return contigs(genome)[contigId];
//}
//
//template <typename TSpec, typename TConfig, typename TContigId>
//inline typename Contigs<Genome<TSpec, TConfig> const>::Type &
//contig(Genome<TSpec, TConfig> const & genome, TContigId contigId)
//{
//    return contigs(genome)[contigId];
//}

// ----------------------------------------------------------------------------
// Function contigLength()                                             [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TContigId>
inline typename Size<Genome<TSpec, TConfig> const>::Type
contigLength(Genome<TSpec, TConfig> const & genome, TContigId contigId)
{
    return length(contigs(genome)[contigId]);
}

// ----------------------------------------------------------------------------
// Function clear()                                                    [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void clear(Genome<TSpec, TConfig> & genome)
{
    clear(contigs(genome));
}

// ----------------------------------------------------------------------------
// Function reserve()                                                  [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize>
void reserve(Genome<TSpec, TConfig> & genome, TSize newCapacity)
{
    reserve(contigs(genome), newCapacity, Exact());
}

// ----------------------------------------------------------------------------
// Function reverse()                                                  [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void reverse(Genome<TSpec, TConfig> & genome)
{
    for (unsigned contigId = 0; contigId < length(contigs(genome)); ++contigId)
        reverse(contigs(genome)[contigId]);
}

// ----------------------------------------------------------------------------
// Function open()                                               [GenomeLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TString>
void open(GenomeLoader<TSpec, TConfig> & loader, TString const & genomeFile)
{
    typedef GenomeLoader<TSpec, TConfig>            TGenomeLoader;
    typedef typename TGenomeLoader::TRecordReader   TRecordReader;

    // Open file.
    loader._file.open(toCString(genomeFile), std::ios::binary | std::ios::in);

    if (!loader._file.is_open())
        throw std::runtime_error("Error while opening genome file.");

    // Compute file size.
    loader._file.seekg(0, std::ios::end);
    loader._fileSize = loader._file.tellg();
    loader._file.seekg(0, std::ios::beg);

    // Initialize record reader.
    loader._reader.reset(new TRecordReader(loader._file));

    // Autodetect file format.
    if (!guessStreamFormat(*(loader._reader), loader._fileFormat))
        throw std::runtime_error("Error while guessing genome file format.");
}

// ----------------------------------------------------------------------------
// Function convertContig()                                      [GenomeLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void convertContig(GenomeLoader<TSpec, TConfig> & /* loader */, Dna5String & /* contig */, Dna5) {}

template <typename TSpec, typename TConfig>
void convertContig(GenomeLoader<TSpec, TConfig> & loader, Dna5String & contig, Dna)
{
    typedef typename Iterator<Dna5String, Standard>::Type   TContigIt;

    TContigIt cIt = begin(contig, Standard());
    TContigIt cEnd = end(contig, Standard());

    while (cIt != cEnd)
    {
        for (; cIt != cEnd && value(cIt) != Dna5('N'); ++cIt) ;

        if (cIt == cEnd) break;

        for (; cIt != cEnd && value(cIt) == Dna5('N'); ++cIt)
            value(cIt) = pickRandomNumber(loader._rng) % ValueSize<Dna>::VALUE;
    }
}

// ----------------------------------------------------------------------------
// Function load()                                               [GenomeLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void load(GenomeLoader<TSpec, TConfig> & loader)
{
    typedef typename TConfig::TContigSeq        TContigSeq;
    typedef typename Value<TContigSeq>::Type    TContigAlphabet;

    load(loader, TContigAlphabet());
}

template <typename TSpec, typename TConfig, typename TAlphabet>
void load(GenomeLoader<TSpec, TConfig> & loader, TAlphabet const & /* tag */)
{
    // Reserve space for contigs.
    reserve(getGenome(loader), loader._fileSize);

    CharString contigName;
    Dna5String contigSeq;

    // Read records.
    while (!atEnd(*(loader._reader)))
    {
        if (readRecord(contigName, contigSeq, *(loader._reader), loader._fileFormat) != 0)
            throw std::runtime_error("Error while reading genome contig.");

        convertContig(loader, contigSeq, TAlphabet());
        appendValue(contigs(getGenome(loader)), contigSeq);
    }
}

#endif  // #ifndef APP_CUDAMAPPER_GENOME_H_
