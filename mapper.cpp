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

#include "tags.h"
#include "reads.h"
#include "genome.h"

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "types.h"
#include "misc.h"
#include "options.h"
#include "index.h"
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

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIREADS FILE\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta fa");
    setValidValues(parser, 1, "fastq fasta fa");

    addSection(parser, "Global Options");

#ifndef CUDA_DISABLED
    addOption(parser, ArgParseOption("nc", "no-cuda", "Do not use CUDA accelerated code."));
#endif

#ifdef _OPENMP
    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", options.threadsCount);
#endif

    addSection(parser, "Mapping Options");

    addOption(parser, ArgParseOption("mb", "mapping-block", "Maximum number of reads to be mapped at once.", ArgParseOption::INTEGER));
    setMinValue(parser, "mapping-block", "1000");
    setDefaultValue(parser, "mapping-block", options.mappingBlock);

    addOption(parser, ArgParseOption("sl", "seed-length", "Minimum seed length.", ArgParseOption::INTEGER));
    setMinValue(parser, "seed-length", "10");
    setMaxValue(parser, "seed-length", "100");
    setDefaultValue(parser, "seed-length", options.seedLength);

    addOption(parser, ArgParseOption("se", "seed-errors", "Maximum number of errors per seed.", ArgParseOption::INTEGER));
    setMinValue(parser, "seed-errors", "0");
    setMaxValue(parser, "seed-errors", "2");
    setDefaultValue(parser, "seed-errors", options.errorsPerSeed);

    addSection(parser, "Genome Index Options");

    setIndexPrefix(parser);
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

    // Parse reads input file.
    getArgumentValue(options.readsFile, parser, 1);

    // Parse CUDA options.
#ifndef CUDA_DISABLED
    getOptionValue(options.noCuda, parser, "no-cuda");
#endif

#ifdef _OPENMP
    // Parse the number of threads.
    getOptionValue(options.threadsCount, parser, "threads");
#endif

    // Parse mapping block.
    getOptionValue(options.mappingBlock, parser, "mapping-block");

    // Parse mapping options.
    getOptionValue(options.seedLength, parser, "seed-length");
    getOptionValue(options.errorsPerSeed, parser, "seed-errors");

    // Parse genome index prefix.
    getIndexPrefix(options, parser);

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function configureMapper()
// ----------------------------------------------------------------------------

int configureMapper(Options const & options)
{
#ifndef CUDA_DISABLED
    if (options.noCuda)
#endif
        spawnMapper(options, ExecHost());
#ifndef CUDA_DISABLED
    else
        spawnMapper(options, ExecDevice());
#endif

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

    return configureMapper(options);
}
