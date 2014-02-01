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

#ifndef APP_CUDAMAPPER_MAPPER_FILTER_H_
#define APP_CUDAMAPPER_MAPPER_FILTER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class FilterDelegate
// ----------------------------------------------------------------------------

template <typename THits, typename TSpec = void>
struct FilterDelegate
{
    typedef typename Value<THits>::Type  THit;
    typedef typename Spec<THit>::Type    THitSpec;

    THits & hits;

    FilterDelegate(THits & hits) :
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
struct View<FilterDelegate<THits, TSpec> >
{
    typedef FilterDelegate<typename View<THits>::Type, TSpec>  Type;
};
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename THits, typename TSpec>
inline typename View<FilterDelegate<THits, TSpec> >::Type
view(FilterDelegate<THits, TSpec> & me)
{
    return typename View<FilterDelegate<THits, TSpec> >::Type(view(me.hits));
}

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename THits, typename TSpec, typename TPattern>
inline void init(FilterDelegate<THits, TSpec> & me, TPattern const & pattern)
{
    typedef FilterDelegate<THits, TSpec>   TManager;
    typedef typename TManager::THitSpec THitSpec;

    _init(me, pattern, THitSpec());
}

template <typename THits, typename TSpec, typename TPattern>
inline void _init(FilterDelegate<THits, TSpec> & me, TPattern const & pattern, Exact)
{
    resize(me.hits, length(needle(pattern)), Exact());
}

template <typename THits, typename TSpec, typename TPattern>
inline void _init(FilterDelegate<THits, TSpec> & me, TPattern const & pattern, HammingDistance)
{
    // TODO(esiragusa): reserve more than this.
    reserve(me.hits, length(needle(pattern)) * 5, Exact());
}

// ----------------------------------------------------------------------------
// Function _addHit()
// ----------------------------------------------------------------------------

template <typename THits, typename TSpec, typename TFinder>
inline SEQAN_HOST_DEVICE void
_addHit(FilterDelegate<THits, TSpec> & me, TFinder const & finder, Exact)
{
    me.hits[finder._patternIt].range = range(textIterator(finder));
}

template <typename THits, typename TSpec, typename TFinder>
inline SEQAN_HOST_DEVICE void
_addHit(FilterDelegate<THits, TSpec> & me, TFinder const & finder, HammingDistance)
{
    typedef typename Value<THits>::Type THit;

    THit hit = { range(textIterator(finder)), finder._patternIt, getScore(finder) };

    appendValue(me.hits, hit, Insist(), Parallel());
}

#endif  // #ifndef APP_CUDAMAPPER_MAPPER_FILTER_H_
