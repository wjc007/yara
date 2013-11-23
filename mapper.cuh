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

#ifndef SEQAN_EXTRAS_CUDAMAPPER_MAPPER_CUH_
#define SEQAN_EXTRAS_CUDAMAPPER_MAPPER_CUH_

// ============================================================================
// Forwards
// ============================================================================

// --------------------------------------------------------------------------
// Function mapReads()                                                  [GPU]
// --------------------------------------------------------------------------

void mapReads(TGenomeIndex & index, TReadSeqs & readSeqs, unsigned seedLength, unsigned errorsPerSeed, ExecDevice const & /* tag */);

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function assign()                                                  [FMIndex]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): We do not assign the text to the device index!

#ifdef PLATFORM_CUDA
namespace seqan {
template <typename TValue, typename TAlloc, typename TSSetSpec, typename TOccSpec, typename TSpec,
          typename TText2, typename TOccSpec2, typename TSpec2>
inline void
assign(Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, FMIndex<TOccSpec, TSpec> > & index,
       Index<TText2, FMIndex<TOccSpec2, TSpec2> > & source)
{
    cudaPrintFreeMemory();

    assign(indexSA(index), indexSA(source));
    assign(indexLF(index), indexLF(source));

    cudaPrintFreeMemory();
}
}
#endif

#endif  // #ifndef SEQAN_EXTRAS_CUDAMAPPER_MAPPER_CUH_
