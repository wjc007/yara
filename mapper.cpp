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
// Masai headers; I/O and options
// ----------------------------------------------------------------------------

#include "tags.h"
#include "options.h"
#include "reads.h"
#include "genome.h"
#include "genome_index.h"

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "types.h"
#include "mapper.h"
#ifndef CUDA_DISABLED
#include "mapper.cuh"
#endif

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString  genomeFile;
    CharString  genomeIndexFile;
    CharString  readsFile;

    bool        noCuda;
    unsigned    threadsCount;
    int         mappingBlock;
    unsigned    seedLength;
    unsigned    errorsPerSeed;

    Options() :
        noCuda(false),
        threadsCount(1),
        mappingBlock(200000),
        seedLength(33),
        errorsPerSeed(0)
    {}
};

// ----------------------------------------------------------------------------
// Class App
// ----------------------------------------------------------------------------

template <typename TExecSpace>
struct App
{
    typedef Genome<void, CUDAStoreConfig>                           TGenome;
    typedef GenomeLoader<void, CUDAStoreConfig>                     TGenomeLoader;
    typedef GenomeIndex<TGenome, TGenomeIndexSpec, void>            TGenomeIndex;
    typedef FragmentStore<void, CUDAStoreConfig>                    TStore;
    typedef ReadsConfig<False, False, True, True, CUDAStoreConfig>  TReadsConfig;
    typedef ReadsLoader<void, TReadsConfig>                         TReadsLoader;

    TGenome             genome;
#ifdef ENABLE_GENOME_LOADING
    TGenomeLoader       genomeLoader;
#endif
    TGenomeIndex        genomeIndex;
    TStore              store;
    TReadsLoader        readsLoader;

    App() :
        genome(),
#ifdef ENABLE_GENOME_LOADING
        genomeLoader(genome),
#endif
        genomeIndex(genome),
        store()
    {};
};

// ----------------------------------------------------------------------------
// Class Timer
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec = void>
struct Timer
{
    TValue _begin, _end;

    Timer() : _begin(0), _end(0) {};
};

template <typename TValue, typename TSpec>
inline void start(Timer<TValue, TSpec> & timer)
{
    timer._begin = sysTime();
}

template <typename TValue, typename TSpec>
inline void stop(Timer<TValue, TSpec> & timer)
{
    timer._end = sysTime();
}

template <typename TValue, typename TSpec>
inline void clear(Timer<TValue, TSpec> & timer)
{
    timer._begin = 0;
    timer._end = 0;
}

template <typename TValue, typename TSpec>
inline TValue getValue(Timer<TValue, TSpec> & timer)
{
    return timer._end - timer._begin;
}

template <typename TValue, typename TSpec>
std::ostream & operator<<(std::ostream & os, Timer<TValue, TSpec> & timer)
{
    os << getValue(timer) << " sec";
    return os;
}

// ----------------------------------------------------------------------------
// Class Logger
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec = void>
struct Logger
{
    TStream &   stream;
    bool        quiet;

    Logger(TStream & stream) :
        stream(stream),
        quiet(false)
    {};
};

template <typename TStream, typename TSpec, typename TStream2>
void write(Logger<TStream, TSpec> & logger, TStream2 const & stream)
{
    if (logger.quiet) return;

    // TODO(esiragusa): use a scoped lock per logger instance.
    SEQAN_OMP_PRAGMA(critical(_logger))
    logger.stream << stream;
}

template <typename TStream, typename TSpec, typename TStream2>
void write(Logger<TStream, TSpec> & logger, TStream2 & stream)
{
    write(logger, reinterpret_cast<TStream2 const &>(stream));
}

template <typename TStream, typename TSpec, typename TStream2>
TStream & operator<<(Logger<TStream, TSpec> & logger, TStream2 & stream)
{
    write(logger, stream);
    return logger.stream;
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()                              [ArgumentParser]
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
// Function parseCommandLine()                                        [Options]
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
// Function runApp()
// ----------------------------------------------------------------------------

template <typename TExecSpace>
void runApp(App<TExecSpace> & app, Options const & options)
{
    Timer<double> timer;

#ifdef _OPENMP
    // Set the number of threads that OpenMP can spawn.
    omp_set_num_threads(options.threadsCount);
    std::cout << "Threads count:\t\t\t" << omp_get_max_threads() << std::endl;
#endif

#ifdef ENABLE_GENOME_LOADING
    // Load genome.
    open(app.genomeLoader, options.genomeFile);

    start(timer, "Loading genome:\t\t\t");
    load(app.genomeLoader);
    stop(timer, " sec");
#endif

    // Load genome index.
    start(timer, "Loading genome index:\t\t");
    load(app.genomeIndex, options.genomeIndexFile);
    stop(timer, " sec");

    // Open reads file.
    open(app.readsLoader, options.readsFile);

    // Reserve space for reads.
    if (options.mappingBlock < MaxValue<int>::VALUE)
        reserve(app.reads, options.mappingBlock);
    else
        reserve(app.reads);

    // Process reads in blocks.
    while (!atEnd(app.readsLoader))
    {
        // Load reads.
        start(timer, "Loading reads:\t\t\t");
        load(app.readsLoader, options.mappingBlock);
        stop(timer, " sec");

        std::cout << "Reads count:\t\t\t" << app.reads.readsCount << std::endl;

        // Map reads.
        mapReads(app.genomeIndex.index, getSeqs(app.reads), options.seedLength, options.errorsPerSeed, TExecSpace());

        // Clear mapped reads.
        clear(app.reads);
    }

    // Close reads file.
    close(app.readsLoader);
}

// ----------------------------------------------------------------------------
// Function runApp()
// ----------------------------------------------------------------------------

#if false
template <typename TExecSpace>
void runApp(App<TExecSpace> & app, Options const & options)
{
    typedef Timer<double>                               TTimer;
    typedef Logger<std::ostream>                        TLogger;
    typedef App<TExecSpace>                             TApp;
    typedef Reads<void, typename TApp::TReadsConfig>    TReads;

    TTimer timer;
    TLogger cout(std::cout);
    TLogger cerr(std::cerr);

#ifdef _OPENMP
    cout << "Threads count:\t\t\t" << options.threadsCount << std::endl;
#endif

#ifdef _OPENMP
    // Disable nested parallelism.
    omp_set_nested(false);
#endif

#ifdef ENABLE_GENOME_LOADING
    // Load genome.
    open(app.genomeLoader, options.genomeFile);

    start(timer);
    load(app.genomeLoader);
    stop(timer);
    cout << "Loading genome:\t\t\t" << timer << std::endl;
#endif

    // Load genome index.
    start(timer);
    load(app.genomeIndex, options.genomeIndexFile);
    stop(timer);
    cout << "Loading genome index:\t\t" << timer << std::endl;

    // Open reads file.
    open(app.readsLoader, options.readsFile);

    // Process reads in parallel.
    SEQAN_OMP_PRAGMA(parallel firstprivate(timer) num_threads(3))
    {
        // Reserve space for reads.
        TReads reads;
        reserve(reads, options.mappingBlock);
//        if (options.mappingBlock < MaxValue<int>::VALUE)
//            reserve(reads, options.mappingBlock);
//        else
//            reserve(reads);

        // Process reads.
        while (true)
        {
            // Load a block of reads.
            SEQAN_OMP_PRAGMA(critical(_mapper_readsLoader_load))
            {
                // No more reads.
                if (!atEnd(app.readsLoader))
                {
                    start(timer);
                    setReads(app.readsLoader, reads);
                    load(app.readsLoader, options.mappingBlock);
                    stop(timer);

                    cout << "Loading reads:\t\t\t" << timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;// <<
//                            "Reads count:\t\t\t" << reads.readsCount << "\t\t\t[" << omp_get_thread_num() << "]" << std::endl;
                }
            }

            // No more reads.
            if (!reads.readsCount) break;

            // Map reads.
            SEQAN_OMP_PRAGMA(critical(_mapper_mapReads))
            {
                #ifdef _OPENMP
                // Enable nested parallelism.
                omp_set_nested(true);
                #endif

                #ifdef _OPENMP
                omp_set_num_threads(options.threadsCount);
                #endif

                start(timer);
                mapReads(app.genomeIndex.index, getSeqs(reads), options.seedLength, options.errorsPerSeed, TExecSpace());
                stop(timer);

                cout << "Mapping reads:\t\t\t" << timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;

                #ifdef _OPENMP
                omp_set_num_threads(1);
                #endif

                #ifdef _OPENMP
                // Disable nested parallelism.
                omp_set_nested(false);
                #endif
            }

            // Writer results.
            SEQAN_OMP_PRAGMA(critical(_mapper_samWriter_write))
            {
                start(timer);
                sleep(reads.readsCount / 1000000.0);
                stop(timer);
                
                cout << "Writing results:\t\t" << timer << "\t\t[" << omp_get_thread_num() << "]" << std::endl;
            }

            // Clear mapped reads.
            clear(reads);
        }
    }

    // Close reads file.
    close(app.readsLoader);
}
#endif

// ----------------------------------------------------------------------------
// Function configureApp()
// ----------------------------------------------------------------------------

template <typename TOptions>
int configureApp(TOptions & options)
{
#ifndef CUDA_DISABLED
    if (options.noCuda)
    {
        App<ExecHost> app;
        runApp(app, options);
    }
    else
    {
        App<ExecDevice> app;
        runApp(app, options);
    }
#else
    App<ExecHost> app;
    runApp(app, options);
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

    return configureApp(options);
}
