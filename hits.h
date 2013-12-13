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
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Hit<Exact>
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec = Exact>
struct Hit
{
    typename Position<Hit>::Type    range;
};

// ----------------------------------------------------------------------------
// Class Hit<HammingDistance>
// ----------------------------------------------------------------------------

template <typename TSize>
struct Hit<TSize, HammingDistance>
{
    typename Position<Hit>::Type    range;
    typename Id<Hit>::Type          seedId;
    unsigned char                   errors;
};

namespace seqan
{

// ----------------------------------------------------------------------------
// Metafunction Size<Hit>
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
struct Size<Hit<TSize, TSpec> >
{
    typedef TSize Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position<Hit>
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
struct Position<Hit<TSize, TSpec> >
{
    typedef Pair<TSize> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Id<Hit>
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
struct Id<Hit<TSize, TSpec> >
{
    typedef unsigned int Type;
};

// ----------------------------------------------------------------------------
// Metafunction Spec<Hit>
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
struct Spec<Hit<TSize, TSpec> >
{
    typedef TSpec Type;
};
}

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

    template <typename THit>
    void operator() (THit const & hit)
    {
        count += getValueI2(hit.range) - getValueI1(hit.range);
    }
};

// ----------------------------------------------------------------------------
// Class HitsManager
// ----------------------------------------------------------------------------

template <typename THits, typename TSpec = void>
struct HitsManager
{
    typedef typename Value<THits>::Type  THit;
    typedef typename Spec<THit>::Type    THitSpec;

    THits & hits;

    HitsManager(THits & hits) :
        hits(hits)
    {}

    template <typename TFinder>
    SEQAN_HOST_DEVICE void
    operator() (TFinder const & finder)
    {
        _addHit(*this, finder, THitSpec());
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

namespace seqan {
template <typename THits, typename TSpec>
struct View<HitsManager<THits, TSpec> >
{
    typedef HitsManager<typename View<THits>::Type, TSpec>  Type;
};
}

// ============================================================================
// Functions
// ============================================================================

template <typename TSize, typename TSpec>
inline bool operator< (Hit<TSize, TSpec> const & a, Hit<TSize, TSpec> const & b)
{
    return a.seedId < b.seedId;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSize>
inline void clear(Hit<TSize, Exact> & hit)
{
    setValueI1(hit.range, 0);
    setValueI2(hit.range, 0);
}

template <typename TSize>
inline void clear(Hit<TSize, HammingDistance> & hit)
{
    setValueI1(hit.range, 0);
    setValueI2(hit.range, 0);
    hit.seedId = 0;
    hit.errors = 0;
}

// ----------------------------------------------------------------------------
// Function getErrors()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
inline unsigned char
getErrors(Hit<TSize, TSpec> const & /* hit */)
{
    return 0;
}

template <typename TSize>
inline unsigned char
getErrors(Hit<TSize, HammingDistance> const & hit)
{
    return hit.errors;
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename THits, typename TSpec>
inline typename View<HitsManager<THits, TSpec> >::Type
view(HitsManager<THits, TSpec> & manager)
{
    return typename View<HitsManager<THits, TSpec> >::Type(view(manager.hits));
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec, typename TPattern>
inline void init(HitsManager<TSize, TSpec> & manager, TPattern const & pattern)
{
    typedef HitsManager<TSize, TSpec>   TManager;
    typedef typename TManager::THitSpec THitSpec;

    _init(manager, pattern, THitSpec());
}

// ----------------------------------------------------------------------------
// Function _init()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec, typename TPattern>
inline void _init(HitsManager<TSize, TSpec> & manager, TPattern const & pattern, Exact)
{
    resize(manager.hits, length(needle(pattern)), Exact());
}

template <typename TSize, typename TSpec, typename TPattern>
inline void _init(HitsManager<TSize, TSpec> & manager, TPattern const & pattern, HammingDistance)
{
    // TODO(esiragusa): reserve more than this.
    reserve(manager.hits, length(needle(pattern)), Exact());
}

// ----------------------------------------------------------------------------
// Function _addHit()
// ----------------------------------------------------------------------------

template <typename THits, typename TSpec, typename TFinder>
inline SEQAN_HOST_DEVICE void
_addHit(HitsManager<THits, TSpec> & manager, TFinder const & finder, Exact)
{
    manager.hits[finder._patternIt].range = range(textIterator(finder));
}

template <typename THits, typename TSpec, typename TFinder>
inline SEQAN_HOST_DEVICE void
_addHit(HitsManager<THits, TSpec> & manager, TFinder const & finder, HammingDistance)
{
    typedef typename Value<THits>::Type THit;

    // TODO(esiragusa): implement getScore(finder) for multiple finder.
    THit hit = { range(textIterator(finder)), finder._patternIt, 1 };

    // TODO(esiragusa): atomic append.
    appendValue(manager.hits, hit);
}

// ----------------------------------------------------------------------------
// Function getHitErrors()
// ----------------------------------------------------------------------------

template <typename THits, typename THitId>
inline unsigned char
getHitErrors(THits const & hits, THitId hitId)
{
    return getErrors(hits[hitId]);
}

// ----------------------------------------------------------------------------
// Function getHitRange()
// ----------------------------------------------------------------------------

template <typename THits, typename THitId>
inline typename Position<typename Value<THits>::Type>::Type
getHitRange(THits const & hits, THitId hitId)
{
    return hits[hitId].range;
}

// ----------------------------------------------------------------------------
// Function getHitIds()
// ----------------------------------------------------------------------------

template <typename THits, typename TSeedId>
inline Pair<typename Id<typename Value<THits>::Type>::Type>
getHitIds(THits const & hits, TSeedId seedId)
{
    typedef typename Value<THits>::Type THit;
    typedef typename Spec<THit>::Type   THitSpec;

    return _getHitIds(hits, seedId, THitSpec());
}

template <typename THits, typename TSeedId>
inline Pair<typename Id<typename Value<THits>::Type>::Type>
_getHitIds(THits const & /* hits */, TSeedId seedId, Exact)
{
    typedef typename Value<THits>::Type THit;
    typedef typename Id<THit>::Type     THitId;
    typedef Pair<THitId>                THitIds;

    return THitIds(seedId, seedId + 1);
}

template <typename THits, typename TSeedId>
inline Pair<typename Id<typename Value<THits>::Type>::Type>
_getHitIds(THits const & hits, TSeedId seedId, HammingDistance)
{
    typedef typename Value<THits>::Type                     THit;
    typedef typename Id<THit>::Type                         THitId;
    typedef Pair<THitId>                                    THitIds;
    typedef typename Iterator<THits const, Standard>::Type  THitsIterator;

    THitsIterator hitsBegin = begin(hits, Standard());
    THitsIterator hitsEnd = end(hits, Standard());

    THit key;
    key.seedId = seedId;

    THitsIterator firstHit = std::lower_bound(hitsBegin, hitsEnd, key);
    THitsIterator lastHit = std::upper_bound(hitsBegin, hitsEnd, key);

    return THitIds(position(firstHit, hits), position(lastHit, hits));
}

// ----------------------------------------------------------------------------
// Function sortHits()
// ----------------------------------------------------------------------------

template <typename THits>
inline void sortHits(THits & hits)
{
    typedef typename Value<THits>::Type THit;
    typedef typename Spec<THit>::Type   THitSpec;

    _sortHits(hits, THitSpec());
}

template <typename THits>
inline void _sortHits(THits & /* hits */, Exact) {}

template <typename THits>
inline void _sortHits(THits & hits, HammingDistance)
{
    return std::stable_sort(begin(hits, Standard()), end(hits, Standard()));
}

// ----------------------------------------------------------------------------
// Function countHits()
// ----------------------------------------------------------------------------

template <typename TSize, typename THits>
inline TSize countHits(THits const & hits)
{
    return std::for_each(begin(hits, Standard()),
                         end(hits, Standard()),
                         HitsCounter<TSize>()).count;
}

// ----------------------------------------------------------------------------
// Function countHits()
// ----------------------------------------------------------------------------

template <typename TSize, typename THits, typename THitId>
inline TSize countHits(THits const & hits, Pair<THitId> hitIds)
{
    return std::for_each(begin(hits, Standard()) + getValueI1(hitIds),
                         begin(hits, Standard()) + getValueI2(hitIds),
                         HitsCounter<TSize>()).count;
}

// ----------------------------------------------------------------------------
// Function clearHits()
// ----------------------------------------------------------------------------

template <typename THits, typename THitId>
inline void clearHits(THits & hits, Pair<THitId> hitIds)
{
    typedef typename Value<THits>::Type THit;

    THit emptyHit;
    clear(emptyHit);

    std::fill(begin(hits, Standard()) + getValueI1(hitIds),
              begin(hits, Standard()) + getValueI2(hitIds),
              emptyHit);
}

#endif  // #ifndef APP_CUDAMAPPER_HITS_H_
