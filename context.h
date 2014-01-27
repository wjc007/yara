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

#ifndef APP_CUDAMAPPER_CONTEXT_H_
#define APP_CUDAMAPPER_CONTEXT_H_

using namespace seqan;

// ============================================================================
// Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Enum ReadStatus
// ----------------------------------------------------------------------------

enum ReadStatus { STATUS_UNSEEDED, STATUS_SEEDED, STATUS_MAPPED, STATUS_UNMAPPABLE };

// ----------------------------------------------------------------------------
// Enum ReadAnchor
// ----------------------------------------------------------------------------

enum ReadAnchor { ANCHOR_FIRST, ANCHOR_SECOND };

// ----------------------------------------------------------------------------
// Class ReadContext
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = void>
struct ReadContext
{
    unsigned char stratum       : 4;
    unsigned char seedErrors    : 2;
    ReadStatus    status        : 2;
//    ReadAnchor    anchor        : 1;

    ReadContext() :
        stratum(0),
        seedErrors(0),
        status(STATUS_UNSEEDED)
    {};
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function ctxIsMapped()
// ----------------------------------------------------------------------------

template <typename TReadContext>
inline bool ctxIsMapped(TReadContext const & readCtx)
{
    return readCtx.status == STATUS_MAPPED;
}

// ----------------------------------------------------------------------------
// Function getStratum()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadSeqId>
inline unsigned char getStratum(TReadsContext const & ctx, TReadSeqId readSeqId)
{
    return ctx[readSeqId].stratum;
}

// ----------------------------------------------------------------------------
// Function incStratum()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadSeqId>
inline void incStratum(TReadsContext & ctx, TReadSeqId readSeqId)
{
    ctx[readSeqId].stratum++;
}

// ----------------------------------------------------------------------------
// Function getSeedErrors()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadSeqId>
inline unsigned char getSeedErrors(TReadsContext const & ctx, TReadSeqId readSeqId)
{
    return ctx[readSeqId].seedErrors;
}

// ----------------------------------------------------------------------------
// Function setSeedErrors()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadSeqId, typename TErrors>
inline void setSeedErrors(TReadsContext & ctx, TReadSeqId readSeqId, TErrors errors)
{
    ctx[readSeqId].seedErrors = errors;
}

// ----------------------------------------------------------------------------
// Function getStatus()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadSeqId>
inline ReadStatus getStatus(TReadsContext const & ctx, TReadSeqId readSeqId)
{
    return ctx[readSeqId].status;
}

// ----------------------------------------------------------------------------
// Function setStatus()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadSeqId>
inline void setStatus(TReadsContext & ctx, TReadSeqId readSeqId, ReadStatus status)
{
    ctx[readSeqId].status = status;
}

// ----------------------------------------------------------------------------
// Function isMapped()
// ----------------------------------------------------------------------------

template <typename TReadsContext, typename TReadSeqId>
inline bool isMapped(TReadsContext const & ctx, TReadSeqId readSeqId)
{
    return ctx[readSeqId].status == STATUS_MAPPED;
}

// ----------------------------------------------------------------------------
// Function countMapped()
// ----------------------------------------------------------------------------

template <typename TReadsContext>
inline unsigned long countMapped(TReadsContext const & ctx)
{
    typedef typename Value<TReadsContext>::Type     TReadContext;

    return std::count_if(begin(ctx, Standard()), end(ctx, Standard()), ctxIsMapped<TReadContext>) / 2;
}

#endif  // #ifndef APP_CUDAMAPPER_CONTEXT_H_
