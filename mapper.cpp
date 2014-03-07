// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
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
#include "misc_options.h"
#include "store_reads.h"
#include "store_genome.h"

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "misc_timer.h"
#include "misc_types.h"
#include "index_fm.h"
#include "bits_hits.h"
#include "bits_context.h"
#include "bits_matches.h"
#include "bits_seeds.h"
#include "find_verifier.h"
#include "find_extender.h"
#include "mapper_collector.h"
#include "mapper_classifier.h"
#include "mapper_ranker.h"
#include "mapper_filter.h"
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
    setAppName(parser, "yara_mapper");
    setShortDescription(parser, "Yara Mapper");
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
    setHelpText(parser, 1, "Either one single-end or two paired-end / mate-pairs read files.");

    addOption(parser, ArgParseOption("v", "verbose", "Displays verbose output."));

    // Setup index options.
    addSection(parser, "Input Options");

    setIndexPrefix(parser);

    // Setup output options.
    addSection(parser, "Output Options");

    setOutputFile(parser, options);

    // Setup mapping options.
    addSection(parser, "Mapping Options");

    addOption(parser, ArgParseOption("e", "error-rate", "Consider mapping locations within this error rate.", ArgParseOption::INTEGER));
    setMinValue(parser, "error-rate", "0");
    setMaxValue(parser, "error-rate", "10");
    setDefaultValue(parser, "error-rate", options.errorRate);

//    addOption(parser, ArgParseOption("s", "strata-rate", "Report found suboptimal mapping locations within this error rate from the optimal one. Note that strata-rate << error-rate.", ArgParseOption::STRING));
//    setMinValue(parser, "strata-rate", "0");
//    setMaxValue(parser, "strata-rate", "10");
//    setDefaultValue(parser, "strata-rate", options.strataRate);

//    addOption(parser, ArgParseOption("a", "all", "Report all suboptimal mapping locations."));// Shortcut for strata-rate = error-rate."));
    addOption(parser, ArgParseOption("s", "strata", "Report only cooptimal mapping locations."));

    // Setup paired-end mapping options.
    addSection(parser, "Paired-End / Mate-Pairs Options");

    addOption(parser, ArgParseOption("ll", "library-length", "Mean template length.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-length", "1");
    setDefaultValue(parser, "library-length", options.libraryLength);

    addOption(parser, ArgParseOption("le", "library-error", "Deviation from the mean template length.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-error", "0");
    setDefaultValue(parser, "library-error", options.libraryError);

    addOption(parser, ArgParseOption("lo", "library-orientation", "Expected orientation of segments in the template.", ArgParseOption::STRING));
    setValidValues(parser, "library-orientation", options.libraryOrientationList);
    setDefaultValue(parser, "library-orientation", options.libraryOrientationList[options.libraryOrientation]);

//    addOption(parser, ArgParseOption("la", "anchor", "Anchor one read and verify its mate."));

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

    addOption(parser, ArgParseOption("r", "reads-count", "Maximum number of reads to process at once.", ArgParseOption::INTEGER));
    setMinValue(parser, "reads-count", "1000");
    setDefaultValue(parser, "reads-count", options.readsCount);
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

    // Parse output file.
    getOutputFile(options.outputFile, options, parser, options.readsFile.i1, "");

    // Parse output format.
    getOutputFormat(options, options.outputFile);

    // Parse genome index prefix.
    getIndexPrefix(options, parser);

    // Parse mapping options.
    getOptionValue(options.errorRate, parser, "error-rate");
//    getOptionValue(options.strataRate, parser, "strata-rate");

    bool strata = false;
    getOptionValue(strata, parser, "strata");
    if (strata) options.mappingMode = STRATA;

    // Parse paired-end mapping options.
    getOptionValue(options.libraryLength, parser, "library-length");
    getOptionValue(options.libraryError, parser, "library-error");
    getOptionValue(options.libraryOrientation, parser, "library-orientation", options.libraryOrientationList);
//    getOptionValue(options.anchorOne, parser, "anchor");

#ifdef _OPENMP
    // Parse the number of threads.
    getOptionValue(options.threadsCount, parser, "threads");
#endif

    // Parse CUDA options.
#ifndef CUDA_DISABLED
    getOptionValue(options.noCuda, parser, "no-cuda");
#endif

    // Parse mapping block option.
    getOptionValue(options.readsCount, parser, "reads-count");

    // Parse verbose output option.
    getOptionValue(options.verbose, parser, "verbose");

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function configureAnchoring()
// ----------------------------------------------------------------------------

//template <typename TExecSpace, typename TThreading, typename TOutputFormat, typename TSequencing, typename TStrategy>
//void configureAnchoring(Options const & options, TExecSpace const & execSpace, TThreading const & threading,
//                        TOutputFormat const & format, TSequencing const & sequencing, TStrategy const & strategy)
//{
//    if (options.anchorOne)
//        spawnMapper(options, execSpace, threading, format, sequencing, strategy, AnchorOne());
//    else
//        spawnMapper(options, execSpace, threading, format, sequencing, strategy, AnchorBoth());
//}

// ----------------------------------------------------------------------------
// Function configureStrategy()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TThreading, typename TOutputFormat, typename TSequencing>
void configureStrategy(Options const & options, TExecSpace const & execSpace, TThreading const & threading,
                       TOutputFormat const & format, TSequencing const & sequencing)
{
    switch (options.mappingMode)
    {
    case STRATA:
        return spawnMapper(options, execSpace, threading, format, sequencing, Strata(), Nothing());

    case ALL:
        return spawnMapper(options, execSpace, threading, format, sequencing, All(), Nothing());

    default:
        return;
    }
}

// ----------------------------------------------------------------------------
// Function configureSequencing()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TThreading, typename TOutputFormat>
void configureSequencing(Options const & options, TExecSpace const & execSpace, TThreading const & threading,
                         TOutputFormat const & format)
{
    if (options.singleEnd)
        configureStrategy(options, execSpace, threading, format, SingleEnd());
    else
        configureStrategy(options, execSpace, threading, format, PairedEnd());

//        configureAnchoring(options, execSpace, threading, format, PairedEnd(), All());
}

// ----------------------------------------------------------------------------
// Function configureOutputFormat()
// ----------------------------------------------------------------------------

template <typename TExecSpace, typename TThreading>
void configureOutputFormat(Options const & options, TExecSpace const & execSpace, TThreading const & threading)
{
    switch (options.outputFormat)
    {
    case SAM:
        return configureSequencing(options, execSpace, threading, Sam());

#ifdef SEQAN_HAS_ZLIB
    case BAM:
        return configureSequencing(options, execSpace, threading, Bam());
#endif

    default:
        return;
    }
}

// ----------------------------------------------------------------------------
// Function configureThreading()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void configureThreading(Options const & options, TExecSpace const & execSpace)
{
#ifdef _OPENMP
    if (options.threadsCount > 1)
        configureOutputFormat(options, execSpace, Parallel());
    else
#endif
        configureOutputFormat(options, execSpace, Serial());
}

// ----------------------------------------------------------------------------
// Function configureMapper()
// ----------------------------------------------------------------------------

void configureMapper(Options const & options)
{
#ifndef CUDA_DISABLED
    if (options.noCuda)
#endif
        configureThreading(options, ExecHost());
#ifndef CUDA_DISABLED
    else
        configureThreading(options, ExecDevice());
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
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}
