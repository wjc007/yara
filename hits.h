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

#ifndef APP_CUDAMAPPER_HITS_H_
#define APP_CUDAMAPPER_HITS_H_

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
    typedef typename Member<Hits, Ranges_>::Type    TRanges;
    typedef typename Size<TRanges>::Type            THitId;
    typedef Pair<THitId>                            THitIds;
    typedef Pair<TSize>                             THitRange;
    typedef unsigned char                           THitErrors;

    TRanges ranges;

    template <typename TFinder>
    inline SEQAN_HOST_DEVICE void
    operator() (TFinder const & finder)
    {
        ranges[finder._patternIt] = range(textIterator(finder));
    }
};

//template <typename TSize>
//struct Hits<TSize, HammingDistance>
//{
//    typename Member<Hits, Ids_>::Type       seedIds;
//    typename Member<Hits, Ranges_>::Type    ranges;
//    typename Member<Hits, Errors_>::Type    errors;
//};

// ----------------------------------------------------------------------------
// Class HitsCounter
// ----------------------------------------------------------------------------

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
// Function getHitErrors()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec, typename THitId>
inline typename Hits<TSize, TSpec>::THitErrors
getHitErrors(Hits<TSize, TSpec> const & /* hits */, THitId /* hitId */)
{
    return 0;
}

// ----------------------------------------------------------------------------
// Function getHitRange()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec, typename THitId>
inline typename Hits<TSize, TSpec>::THitRange
getHitRange(Hits<TSize, TSpec> const & hits, THitId hitId)
{
    return hits.ranges[hitId];
}

// ----------------------------------------------------------------------------
// Function getHitIds()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec, typename TSeedId>
inline typename Hits<TSize, TSpec>::THitIds
getHitIds(Hits<TSize, TSpec> const & /* hits */, TSeedId seedId)
{
    typedef Hits<TSize, TSpec> const    THits;
    typedef typename THits::THitIds     THitIds;

    return THitIds(seedId, seedId + 1);
}

// ----------------------------------------------------------------------------
// Function countHits()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
inline unsigned long countHits(Hits<TSize, TSpec> const & hits)
{
    return std::for_each(begin(hits.ranges, Standard()),
                         end(hits.ranges, Standard()),
                         HitsCounter<unsigned long, TSpec>()).count;
}

// ----------------------------------------------------------------------------
// Function countHits()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec, typename TSeedId>
inline TSize countHits(Hits<TSize, TSpec> const & hits, Pair<TSeedId> seedIds)
{
    return std::for_each(begin(hits.ranges, Standard()) + getValueI1(seedIds),
                         begin(hits.ranges, Standard()) + getValueI2(seedIds),
                         HitsCounter<TSize, TSpec>()).count;
}

// ----------------------------------------------------------------------------
// Function clearHits()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
inline void clearHits(Hits<TSize, TSpec> & hits)
{
    clear(hits.ranges);
}

// ----------------------------------------------------------------------------
// Function clearHits()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec, typename TSeedId>
inline void clearHits(Hits<TSize, TSpec> & hits, Pair<TSeedId> seedIds)
{
    Pair<TSize> emptyRange;
    setValueI1(emptyRange, 0);
    setValueI2(emptyRange, 0);

    std::fill(begin(hits.ranges, Standard()) + getValueI1(seedIds),
              begin(hits.ranges, Standard()) + getValueI2(seedIds),
              emptyRange);
}

#endif  // #ifndef APP_CUDAMAPPER_HITS_H_
