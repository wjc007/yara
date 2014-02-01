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

// ============================================================================
// Forwards
// ============================================================================

struct Options;

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/parallel.h>

// ----------------------------------------------------------------------------
// I/O and options
// ----------------------------------------------------------------------------

#include "misc_tags.h"
#include "store_reads.h"
#include "store_genome.h"

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "misc_timer.h"
#include "misc_options.h"
#include "misc_types.h"
#include "index_fm.h"
#include "bits_hits.h"
#include "bits_context.h"
#include "bits_matches.h"
#include "bits_seeds.h"
#include "find_verifier.h"
#include "find_extender.h"
#include "mapper_seeder.h"
#include "mapper_classifier.h"
#include "mapper_extender.h"
#include "mapper_verifier.h"
#include "mapper_writer.h"
#include "mapper.h"
#ifndef CUDA_DISABLED
#include "mapper.cuh"
#endif

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "cuda_mapper");
    setShortDescription(parser, "CUDA Mapper");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    // Setup mandatory arguments.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE FILE\\fP> <\\fISE-READS FILE\\fP>");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE FILE\\fP> <\\fIPE-READS FILE 1\\fP> <\\fIPE-READS FILE 2\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta fa");
    setHelpText(parser, 0, "A reference genome file.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "READS", true));
    setValidValues(parser, 1, "fastq fasta fa");
    setHelpText(parser, 1, "Either one single-end or two paired-end read files.");

    addOption(parser, ArgParseOption("v", "verbose", "Displays verbose output."));

    // Setup mapping options.
    addSection(parser, "Mapping Options");

    addOption(parser, ArgParseOption("mm", "mapping-mode", "Selects a mapping strategy.", ArgParseOption::STRING));
    setValidValues(parser, "mapping-mode", options.mappingModeList);
    setDefaultValue(parser, "mapping-mode", options.mappingModeList[options.mappingMode]);

    addOption(parser, ArgParseOption("e", "error-rate", "Considers mapping locations up to this specified error rate.", ArgParseOption::INTEGER));
    setMinValue(parser, "error-rate", "0");
    setMaxValue(parser, "error-rate", "10");
    setDefaultValue(parser, "error-rate", options.errorRate);

    addOption(parser, ArgParseOption("ll", "library-length", "Paired-end library length.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-length", "1");
    setDefaultValue(parser, "library-length", options.libraryLength);

    addOption(parser, ArgParseOption("le", "library-error", "Paired-end library length tolerance.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-error", "0");
    setDefaultValue(parser, "library-error", options.libraryError);

    addOption(parser, ArgParseOption("a", "anchor", "Anchor one read and verify its mate."));

    // Setup index options.
    addSection(parser, "Indexing Options");

    setIndexPrefix(parser);

    // Setup performance options.
    addSection(parser, "Performance Options");

#ifdef _OPENMP
    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", options.threadsCount);
#endif

#ifndef CUDA_DISABLED
    addOption(parser, ArgParseOption("nc", "no-cuda", "Do not use CUDA accelerated code."));
#endif

    addOption(parser, ArgParseOption("mb", "mapping-block", "Maximum number of reads to be mapped at once.", ArgParseOption::INTEGER));
    setMinValue(parser, "mapping-block", "1000");
    setDefaultValue(parser, "mapping-block", options.mappingBlock);
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse genome input file.
    getArgumentValue(options.genomeFile, parser, 0);

    // Parse read input files.
    switch (getArgumentValueCount(parser, 1))
    {
    case 1:
        getArgumentValue(options.readsFile.i1, parser, 1, 0);
        options.singleEnd = true;
        break;
    case 2:
        getArgumentValue(options.readsFile.i1, parser, 1, 0);
        getArgumentValue(options.readsFile.i2, parser, 1, 1);
        options.singleEnd = false;
        break;
    default:
        return ArgumentParser::PARSE_ERROR;
    }

    // Parse mapping options.
    getOptionValue(options.mappingMode, parser, "mapping-mode", options.mappingModeList);
    getOptionValue(options.errorRate, parser, "error-rate");
    getOptionValue(options.libraryLength, parser, "library-length");
    getOptionValue(options.libraryError, parser, "library-error");
    getOptionValue(options.anchorOne, parser, "anchor");

    // Parse genome index prefix.
    getIndexPrefix(options, parser);

#ifdef _OPENMP
    // Parse the number of threads.
    getOptionValue(options.threadsCount, parser, "threads");
#endif

    // Parse CUDA options.
#ifndef CUDA_DISABLED
    getOptionValue(options.noCuda, parser, "no-cuda");
#endif

    // Parse mapping block option.
    getOptionValue(options.mappingBlock, parser, "mapping-block");

    // Parse verbose output option.
    getOptionValue(options.verbose, parser, "verbose");

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function configureAnchoring()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TSequencing, typename TStrategy>
void configureAnchoring(Options const & options, TExecSpace const & execSpace, TSequencing const & sequencing, TStrategy const & strategy)
{
    if (options.anchorOne)
        spawnMapper(options, execSpace, sequencing, strategy, AnchorOne());
    else
        spawnMapper(options, execSpace, sequencing, strategy, AnchorBoth());
}

// ----------------------------------------------------------------------------
// Function configureStrategy()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TSequencing>
void configureStrategy(Options const & options, TExecSpace const & execSpace, TSequencing const & sequencing)
{
    switch (options.mappingMode)
    {
    case Options::ANY_BEST:
        return spawnMapper(options, execSpace, sequencing, AnyBest(), Nothing());

    case Options::ALL_BEST:
        return spawnMapper(options, execSpace, sequencing, AllBest(), Nothing());

    case Options::ALL:
        return spawnMapper(options, execSpace, sequencing, All(), Nothing());

    default:
        return;
    }
}

// ----------------------------------------------------------------------------
// Function configureSequencing()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void configureSequencing(Options const & options, TExecSpace const & execSpace)
{
    if (options.singleEnd)
        configureStrategy(options, execSpace, SingleEnd());
    else
        configureAnchoring(options, execSpace, PairedEnd(), All());
}

// ----------------------------------------------------------------------------
// Function configureMapper()
// ----------------------------------------------------------------------------

void configureMapper(Options const & options)
{
#ifndef CUDA_DISABLED
    if (options.noCuda)
#endif
        configureSequencing(options, ExecHost());
#ifndef CUDA_DISABLED
    else
        configureSequencing(options, ExecDevice());
#endif
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

    try
    {
        configureMapper(options);
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
