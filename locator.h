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

#ifndef APP_CUDAMAPPER_LOCATOR_H_
#define APP_CUDAMAPPER_LOCATOR_H_

using namespace seqan;

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Ranges_
// ----------------------------------------------------------------------------

struct Ranges_;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Hits
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec = void>
struct Hits
{
    typename Member<Hits, Ranges_>::Type    ranges;

    template <typename TFinder>
    inline SEQAN_HOST_DEVICE void
    operator() (TFinder const & finder)
    {
        ranges[finder._patternIt] = range(textIterator(finder));
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Member Ranges_
// ----------------------------------------------------------------------------

struct Ranges_;

namespace seqan {
template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, TSpec>, Ranges_>
{
    typedef Pair<typename Size<TIndex>::Type>   TRange_;
    typedef String<TRange_>                     Type;
};

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, Device<TSpec> >, Ranges_>
{
    typedef Pair<typename Size<TIndex>::Type>   TRange_;
    typedef thrust::device_vector<TRange_>      Type;
};
#endif

template <typename TIndex, typename TSpec>
struct Member<Hits<TIndex, View<TSpec> >, Ranges_>
{
    typedef typename Member<Hits<TIndex, TSpec>, Ranges_>::Type  TRanges_;
    typedef typename View<TRanges_>::Type   Type;
};
}

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TIndex, typename TSpec>
struct View<Hits<TIndex, TSpec> >
{
    typedef Hits<TIndex, View<TSpec> >  Type;
};
}

// ----------------------------------------------------------------------------
// Metafunction Device
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TIndex, typename TSpec>
struct Device<Hits<TIndex, TSpec> >
{
    typedef Hits<TIndex, Device<TSpec> >  Type;
};
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TPattern>
inline void init(Hits<TIndex, TSpec> & hits, TPattern const & pattern)
{
    resize(hits.ranges, length(needle(pattern)), Exact());
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
inline typename View<Hits<TIndex, TSpec> >::Type
view(Hits<TIndex, TSpec> & hits)
{
    typename View<Hits<TIndex, TSpec> >::Type hitsView;

    hitsView.ranges = view(hits.ranges);

    return hitsView;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
inline void clear(Hits<TIndex, TSpec> & hits)
{
    clear(hits.ranges);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSize>
inline bool isValid(Pair<TSize> range)
{
    return range.i1 < range.i2;
}

template <typename TIndex, typename TSpec>
inline typename Size<TIndex>::Type
getCount(Hits<TIndex, TSpec> const & hits)
{
    return std::count_if(begin(hits.ranges, Standard()), end(hits.ranges, Standard()), isValid<typename Size<TIndex>::Type>);
}

#endif  // #ifndef APP_CUDAMAPPER_LOCATOR_H_
