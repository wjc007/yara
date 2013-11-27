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
// This file contains the cuda_indexer application.
// ==========================================================================

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/index.h>

// ----------------------------------------------------------------------------
// I/O and options
// ----------------------------------------------------------------------------

#include "genome.h"

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "types.h"
#include "misc.h"
#include "options.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString genomeFile;
    CharString genomeIndexFile;
};

// ----------------------------------------------------------------------------
// Class Indexer
// ----------------------------------------------------------------------------

template <typename TIndexSpec, typename TSpec = void>
struct Indexer
{
    typedef Genome<void, CUDAStoreConfig>                           TGenome;
    typedef GenomeLoader<void, CUDAStoreConfig>                     TGenomeLoader;
    typedef Index<typename Contigs<TGenome>::Type, TIndexSpec>      TIndex;

    TGenome             genome;
    TGenomeLoader       genomeLoader;
    TIndex              index;
    Timer<double>       timer;

    Indexer() :
        genome(),
        genomeLoader(genome),
        index(contigs(genome))
    {};
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser, Options const & /* options */)
{
    setAppName(parser, "cuda_indexer");
    setShortDescription(parser, "CUDA Indexer");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta fa");

    addSection(parser, "Genome Index Options");

    setIndexPrefix(parser);
    
    addSection(parser, "Output Options");

    setTmpFolder(parser);
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Parse genome input file.
    getArgumentValue(options.genomeFile, parser, 0);

    // Parse genome index prefix.
    getIndexPrefix(options, parser);

    // Parse tmp folder.
    getTmpFolder(options, parser);

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function loadGenome()
// ----------------------------------------------------------------------------

template <typename TIndexSpec, typename TSpec>
void loadGenome(Indexer<TIndexSpec, TSpec> & indexer, Options const & options)
{
    std::cout << "Loading genome:\t\t\t" << std::flush;
    start(indexer.timer);

    open(indexer.genomeLoader, options.genomeFile);
    load(indexer.genomeLoader);

    stop(indexer.timer);
    std::cout << indexer.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function buildIndex()
// ----------------------------------------------------------------------------

template <typename TIndexSpec, typename TSpec>
void buildIndex(Indexer<TIndexSpec, TSpec> & indexer)
{
    typedef typename Indexer<TIndexSpec, TSpec>::TIndex TIndex;

    std::cout << "Building genome index:\t\t" << std::flush;
    start(indexer.timer);

    // IndexFM is built on the reversed genome.
    reverse(contigs(indexer.genome));

    // Set the Index text.
//    setValue(indexer.index.text, contigs(indexer.genome));

    // Iterator instantiation calls automatic index construction.
    typename Iterator<TIndex, TopDown<> >::Type it(indexer.index);
    ignoreUnusedVariableWarning(it);

    reverse(contigs(indexer.genome));

    stop(indexer.timer);
    std::cout << indexer.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function saveIndex()
// ----------------------------------------------------------------------------

template <typename TIndexSpec, typename TSpec>
void saveIndex(Indexer<TIndexSpec, TSpec> & indexer, Options const & options)
{
    std::cout << "Dumping genome index:\t\t" << std::flush;
    start(indexer.timer);

    if (!save(indexer.index, toCString(options.genomeIndexFile)))
        throw RuntimeError("Error while dumping genome index file.");

    stop(indexer.timer);
    std::cout << indexer.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function runIndexer()
// ----------------------------------------------------------------------------

template <typename TIndexSpec, typename TSpec>
int runIndexer(Indexer<TIndexSpec, TSpec> & indexer, Options const & options)
{
    loadGenome(indexer, options);
    buildIndex(indexer);
    saveIndex(indexer, options);

    return 0;
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    Indexer<TGenomeIndexSpec, void> indexer;
    return runIndexer(indexer, options);
}
