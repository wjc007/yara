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

#ifndef APP_CUDAMAPPER_MISC_H_
#define APP_CUDAMAPPER_MISC_H_

using namespace seqan;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Space
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec = void>
struct Space
{
    typedef TObject Type;
};

template <typename TObject>
struct Space<TObject, ExecDevice>
{
    typedef typename Device<TObject>::Type  Type;
};

// ============================================================================
// Classes
// ============================================================================

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

#endif // APP_CUDAMAPPER_MISC_H_
