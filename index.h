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
// This file contains Index overloads.
// ==========================================================================

#ifndef APP_CUDAMAPPER_INDEX_H_
#define APP_CUDAMAPPER_INDEX_H_

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------
// This function is overloaded to avoid loading the text.

namespace seqan {
template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig>
inline bool open(Index<StringSet<TText, TSSetSpec>, FMIndex<TSpec, TConfig> > & index,
                 const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;    append(name, ".sa");
    if (!open(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!open(getFibre(index, FibreLF()), toCString(name), openMode)) return false;

    setFibre(getFibre(index, FibreSA()), getFibre(index, FibreLF()), FibreLF());

    return true;
}
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------
// This function is overloaded to load a CUDA index.
// TODO(esiragusa): remove this as soon as open() works on thrust::device_vector

#ifdef PLATFORM_CUDA
namespace seqan {
template <typename TValue, typename TTraits, typename TSSetSpec, typename TSpec, typename TConfig>
inline bool open(Index<StringSet<thrust::device_vector<TValue, TTraits>, TSSetSpec>, FMIndex<TSpec, TConfig> > & index,
                 const char * fileName, int openMode)
{
    typedef StringSet<String<TValue>, TSSetSpec>    TText;
    typedef Index<TText, FMIndex<TSpec, TConfig> >  THostIndex;

    THostIndex hindex;
    if (!open(hindex, fileName, openMode)) return false;

    assign(index, hindex);

    return true;
}
}
#endif

// ----------------------------------------------------------------------------
// Function assign()                                                  [FMIndex]
// ----------------------------------------------------------------------------
// This function is overloaded to avoid the text.

namespace seqan {
template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig, typename TText2, typename TSSetSpec2,
          typename TOccSpec2, typename TSpec2>
inline void
assign(Index<StringSet<TText, TSSetSpec>, FMIndex<TSpec, TConfig> > & index,
       Index<StringSet<TText2, TSSetSpec2>, FMIndex<TOccSpec2, TSpec2> > const & source)
{
    assign(indexLF(index), indexLF(source));
    assign(indexSA(index), indexSA(source));

    // Set the CSA pointer to LF.
    setFibre(indexSA(index), indexLF(index), FibreLF());
}
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------
// This function is overloaded to avoid the text.

namespace seqan {
template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig>
typename View<Index<StringSet<TText, TSSetSpec>, FMIndex<TSpec, TConfig> > >::Type
view(Index<StringSet<TText, TSSetSpec>, FMIndex<TSpec, TConfig> > & index)
{
    typename View<Index<StringSet<TText, TSSetSpec>, FMIndex<TSpec, TConfig> > >::Type indexView;

    indexLF(indexView) = view(indexLF(index));
    indexSA(indexView) = view(indexSA(index));

    return indexView;
}
}

#endif  // #ifndef APP_CUDAMAPPER_INDEX_H_
