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

template <typename TSize, typename TSpec = void>
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

template <typename TSize, typename TSpec = void>
struct HitsCounter
{
    TSize count;

    HitsCounter() :
        count(0)
    {}

    void operator() (Pair<TSize> const & range)
    {
        count += getValueI2(range) - getValueI1(range);
    }
};

template <typename TSize, typename TSpec = void>
struct IsHardToMap
{
    typedef Hits<TSize, TSpec>         THits;

    THits const &   hits;
    TSize           threshold;

    template <typename TTreshold>
    IsHardToMap(THits const & hits, TTreshold threshold) :
        hits(hits),
        threshold(threshold)
    {}

    template <typename TReadId>
    inline SEQAN_HOST_DEVICE bool
    operator() (TReadId readId)
    {
        return getCount(hits, readId) > threshold;
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
template <typename TSize, typename TSpec>
struct Member<Hits<TSize, TSpec>, Ranges_>
{
    typedef Pair<TSize>                         TRange_;
    typedef String<TRange_>                     Type;
};

#ifdef PLATFORM_CUDA
template <typename TSize, typename TSpec>
struct Member<Hits<TSize, Device<TSpec> >, Ranges_>
{
    typedef Pair<TSize>                         TRange_;
    typedef thrust::device_vector<TRange_>      Type;
};
#endif

template <typename TSize, typename TSpec>
struct Member<Hits<TSize, View<TSpec> >, Ranges_>
{
    typedef typename Member<Hits<TSize, TSpec>, Ranges_>::Type  TRanges_;
    typedef typename View<TRanges_>::Type   Type;
};
}

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSize, typename TSpec>
struct View<Hits<TSize, TSpec> >
{
    typedef Hits<TSize, View<TSpec> >  Type;
};
}

// ----------------------------------------------------------------------------
// Metafunction Device
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSize, typename TSpec>
struct Device<Hits<TSize, TSpec> >
{
    typedef Hits<TSize, Device<TSpec> >  Type;
};
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec, typename TPattern>
inline void init(Hits<TSize, TSpec> & hits, TPattern const & pattern)
{
    resize(hits.ranges, length(needle(pattern)), Exact());
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
inline typename View<Hits<TSize, TSpec> >::Type
view(Hits<TSize, TSpec> & hits)
{
    typename View<Hits<TSize, TSpec> >::Type hitsView;

    hitsView.ranges = view(hits.ranges);

    return hitsView;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
inline void clear(Hits<TSize, TSpec> & hits)
{
    clear(hits.ranges);
}

//template <typename TSize, typename TSpec>
//inline TSize getCount(Hits<TSize, TSpec> const & hits, TReadId readId)
//{
//    return std::count_if(begin(hits.ranges, Standard()), end(hits.ranges, Standard()), isValid<TSize>);
//}

// ----------------------------------------------------------------------------
// Function getCount()
// ----------------------------------------------------------------------------

template <typename TSize>
inline bool isValid(Pair<TSize> const & range)
{
    return range.i1 < range.i2;
}

template <typename TSize, typename TSpec>
inline TSize countRanges(Hits<TSize, TSpec> const & hits)
{
    return std::count_if(begin(hits.ranges, Standard()), end(hits.ranges, Standard()), isValid<TSize>);
}

template <typename TSize, typename TSpec>
inline unsigned long countHits(Hits<TSize, TSpec> const & hits)
{
    typedef HitsCounter<unsigned long, TSpec>   TCounter;

    return std::for_each(begin(hits.ranges, Standard()), end(hits.ranges, Standard()), TCounter()).count;
}

template <typename TSize, typename TSpec, typename TReadId>
inline TSize countHits(Hits<TSize, TSpec> const & hits, TReadId readId)
{
    typedef Hits<TSize, TSpec> const                    THits;
    typedef HitsCounter<TSize, TSpec>                   THitsCounter;
    typedef typename Member<THits, Ranges_>::Type       TRanges;
    typedef typename Iterator<TRanges const, Standard>::Type    TRangesIt;

    TRangesIt rangesBegin = begin(hits.ranges, Standard()) + (readId * (5u + 1));
    TRangesIt rangesEnd = begin(hits.ranges, Standard()) + ((readId + 1) * (5u + 1));

    return std::for_each(rangesBegin, rangesEnd, THitsCounter()).count;
}


#endif  // #ifndef APP_CUDAMAPPER_LOCATOR_H_
