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
// This file contains type definitions.
// ==========================================================================

#ifndef APP_CUDAMAPPER_TYPES_H_
#define APP_CUDAMAPPER_TYPES_H_

using namespace seqan;

// ============================================================================
// Global Types
// ============================================================================

// ----------------------------------------------------------------------------
// Fragment Store Configuration
// ----------------------------------------------------------------------------

struct CUDAStoreConfig
{
    typedef String<Dna>             TReadSeq;
    typedef String<Dna>             TContigSeq;
    typedef Owner<ConcatDirect<> >  TContigSpec;

    typedef double                  TMean;
    typedef double                  TStd;
    typedef signed char             TMappingQuality;

    typedef void                    TReadStoreElementSpec;
    typedef Owner<ConcatDirect<> >  TReadSeqStoreSpec;
    typedef Alloc<>                 TReadNameSpec;
    typedef Owner<ConcatDirect<> >  TReadNameStoreSpec;
    typedef void                    TMatePairStoreElementSpec;
    typedef void                    TLibraryStoreElementSpec;
    typedef void                    TContigStoreElementSpec;
    typedef void                    TContigFileSpec;
    typedef void                    TAlignedReadStoreElementSpec;
    typedef Owner<ConcatDirect<> >  TAlignedReadTagStoreSpec;
    typedef void                    TAnnotationStoreElementSpec;
};

// ============================================================================
// Other Store Types
// ============================================================================

typedef StringSet<CUDAStoreConfig::TContigSeq, CUDAStoreConfig::TContigSpec>        TContigs;
typedef StringSet<CUDAStoreConfig::TReadSeq, CUDAStoreConfig::TReadSeqStoreSpec>    TReadSeqs;

// ----------------------------------------------------------------------------
// ReadSeqs Size
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Register usage does not decrease.

//namespace seqan {
//template <>
//struct Size<TReadSeqs>
//{
//    typedef __uint32 Type;
//};
//
//#ifdef PLATFORM_CUDA
//template <>
//struct Size<typename Device<TReadSeqs>::Type>
//{
//    typedef __uint32 Type;
//};
//
//template <>
//struct Size<typename View<typename Device<TReadSeqs>::Type>::Type>
//{
//    typedef __uint32 Type;
//};
//#endif
//}

// ----------------------------------------------------------------------------
// ContainerView Size
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Register usage does not decrease.

//#ifdef PLATFORM_CUDA
//namespace seqan {
//template <typename TAlloc>
//struct Size<thrust::device_vector<Dna, TAlloc> >
//{
//    typedef __uint32 Type;
//};
//
//template <typename TAlloc>
//struct Size<thrust::device_vector<bool, TAlloc> >
//{
//    typedef __uint32 Type;
//};
//
//template <typename TValue, typename TAlloc, typename TViewSpec>
//struct Size<ContainerView<thrust::device_vector<TValue, TAlloc>, TViewSpec> >
//{
//    typedef __uint32 Type;
//};
//
//template <typename TValue, typename TAlloc, typename TViewSpec, typename TSSetSpec>
//struct Size<StringSet<ContainerView<thrust::device_vector<TValue, TAlloc>, TViewSpec>, TSSetSpec> >
//{
//    typedef __uint32 Type;
//};
//}
//#endif

// ============================================================================
// Index Types
// ============================================================================

// ----------------------------------------------------------------------------
// Suffix Array Value Type
// ----------------------------------------------------------------------------

namespace seqan {
template <>
struct SAValue<TContigs>
{
    typedef Pair<__uint8, __uint32, Pack> Type;
};
}

// ----------------------------------------------------------------------------
// FM Index Fibres
// ----------------------------------------------------------------------------

struct CUDAFMIndexConfig
{
    typedef TwoLevels<void>    TValuesSpec;
    typedef Naive<void>        TSentinelsSpec;

    static const unsigned SAMPLING = 10;
};

typedef FMIndex<void, CUDAFMIndexConfig>        TGenomeIndexSpec;
typedef Index<TContigs, TGenomeIndexSpec>       TGenomeIndex;

// ----------------------------------------------------------------------------
// FM Index Size
// ----------------------------------------------------------------------------

namespace seqan {
template <>
struct Size<TGenomeIndex>
{
    typedef __uint32 Type;
};

template <>
struct Size<typename View<TGenomeIndex>::Type>
{
    typedef __uint32 Type;
};

#ifdef PLATFORM_CUDA
template <>
struct Size<typename Device<TGenomeIndex>::Type>
{
    typedef __uint32 Type;
};

template <>
struct Size<typename View<typename Device<TGenomeIndex>::Type>::Type>
{
    typedef __uint32 Type;
};
#endif
}

// ----------------------------------------------------------------------------
// FibreLF Size
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSpec, typename TConfig>
struct Size<LF<TContigs, TSpec, TConfig> >
{
    typedef __uint32 Type;
};
}

// ----------------------------------------------------------------------------
// Rank Dictionary Size
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSpec>
struct Size<RankDictionary<Dna, TwoLevels<TSpec> > >
{
    typedef __uint32 Type;
};

template <typename TSpec>
struct Size<RankDictionary<bool, TwoLevels<TSpec> > >
{
    typedef __uint32 Type;
};

template <typename TSpec>
struct Size<RankDictionary<bool, Naive<TSpec> > >
{
    typedef __uint32 Type;
};
}

// ----------------------------------------------------------------------------
// CSA Size
// ----------------------------------------------------------------------------

//namespace seqan {
//// TODO(esiragusa): Overload Size<CSA> instead of Size<SparseString>
//template <typename TValueString>
//struct Size<SparseString<TValueString, void> >
//{
//    typedef __uint32    Type;
//};
//}

// ----------------------------------------------------------------------------
// Shape Size
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TValue, unsigned q>
struct Value<Shape<TValue, UngappedShape<q> > >
{
    typedef __uint32    Type;
};

// ----------------------------------------------------------------------------
// StringSet MinLength
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct MinLength<StringSet<TString, TSpec> >
{
    static const unsigned VALUE = 10;
};

#ifdef PLATFORM_CUDA
template <typename TValue, typename TSpec>
struct MinLength<thrust::device_vector<TValue, TSpec> >
{
    static const unsigned VALUE = 10;
};
#endif

template <typename TContainer, typename TSpec>
struct MinLength<ContainerView<TContainer, TSpec> >
{
    static const unsigned VALUE = 10;
};
}

#endif  // #ifndef APP_CUDAMAPPER_TYPES_H_
