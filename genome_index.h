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
// This file contains the GenomeIndex class.
// ==========================================================================

#ifndef APP_CUDAMAPPER_GENOME_INDEX_H_
#define APP_CUDAMAPPER_GENOME_INDEX_H_

#include <stdexcept>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class GenomeIndex
// ----------------------------------------------------------------------------

template <typename TGenome, typename TIndexSpec, typename TSpec = void>
struct GenomeIndex
{
    typedef Index<typename Contigs<TGenome>::Type, TIndexSpec>  TIndex;

    Holder<TGenome>     genome;
    TIndex              index;

    GenomeIndex() {}

    GenomeIndex(TGenome & genome) :
        genome(genome)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GenomeHost<T>::Type                               [GenomeIndex]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TIndexSpec, typename TSpec>
struct GenomeHost<GenomeIndex<TGenome, TIndexSpec, TSpec> >
{
    typedef TGenome Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function load()                                                [GenomeIndex]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TIndexSpec, typename TSpec, typename TString>
void load(GenomeIndex<TGenome, TIndexSpec, TSpec> & genomeIndex, TString const & genomeIndexFile)
{
    setValue(genomeIndex.index.text, contigs(getGenome(genomeIndex)));

    if (!open(genomeIndex.index, toCString(genomeIndexFile)))
        throw std::runtime_error("Error while opening genome index file.");
}

// ----------------------------------------------------------------------------
// Function build()                                               [GenomeIndex]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TIndexSpec, typename TSpec>
void build(GenomeIndex<TGenome, TIndexSpec, TSpec> & genomeIndex)
{
    typedef typename GenomeIndex<TGenome, TIndexSpec, TSpec>::TIndex    TIndex;

    setValue(genomeIndex.index.text, contigs(getGenome(genomeIndex)));

    // Iterator instantiation calls automatic index construction.
    typename Iterator<TIndex, TopDown<> >::Type it(genomeIndex.index);
}

template <typename TGenome, typename TIndexSpec, typename TIndexConfig, typename TSpec>
void build(GenomeIndex<TGenome, FMIndex<TIndexSpec, TIndexConfig>, TSpec> & genomeIndex)
{
    typedef typename GenomeIndex<TGenome, FMIndex<TIndexSpec, TIndexConfig>, TSpec>::TIndex   TIndex;

    // IndexFM is built on the reversed genome.
    reverse(getGenome(genomeIndex));

    setValue(genomeIndex.index.text, contigs(getGenome(genomeIndex)));

    // Iterator instantiation calls automatic index construction.
    typename Iterator<TIndex, TopDown<> >::Type it(genomeIndex.index);

    ignoreUnusedVariableWarning(it);

    reverse(getGenome(genomeIndex));
}

// ----------------------------------------------------------------------------
// Function dump()                                                [GenomeIndex]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TIndexSpec, typename TSpec, typename TString>
void dump(GenomeIndex<TGenome, TIndexSpec, TSpec> & genomeIndex, TString const & genomeIndexFile)
{
    if (!save(genomeIndex.index, toCString(genomeIndexFile)))
        throw std::runtime_error("Error while saving genome index file.");
}

// ----------------------------------------------------------------------------
// Function clear()                                               [GenomeIndex]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TIndexSpec, typename TSpec>
void clear(GenomeIndex<TGenome, TIndexSpec, TSpec> & genomeIndex)
{
    clear(genomeIndex.index);
}

#endif  // #ifndef APP_CUDAMAPPER_GENOME_INDEX_H_
