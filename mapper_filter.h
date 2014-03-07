// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
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

#ifndef APP_YARA_MAPPER_FILTER_H_
#define APP_YARA_MAPPER_FILTER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class FilterDelegate
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
struct FilterDelegate
{
    typedef typename Traits::THits      THits;
    typedef typename Value<THits>::Type THit;
    typedef typename Spec<THit>::Type   THitSpec;

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
template <typename TSpec, typename Traits>
struct View<FilterDelegate<TSpec, Traits> >
{
    typedef FilterDelegate<typename View<typename Traits::THits>::Type, TSpec>   Type;
};
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
inline typename View<FilterDelegate<TSpec, Traits> >::Type
view(FilterDelegate<TSpec, Traits> & me)
{
    return typename View<FilterDelegate<TSpec, Traits> >::Type(view(me.hits));
}

// ----------------------------------------------------------------------------
// Function _addHit()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TFinder>
inline SEQAN_HOST_DEVICE void
_addHit(FilterDelegate<TSpec, Traits> & me, TFinder const & finder, Exact)
{
    // NOTE(esiragusa): resize(hits, length(needle(pattern)), Exact())
    me.hits[finder._patternIt].range = range(textIterator(finder));
}

template <typename TSpec, typename Traits, typename TFinder>
inline SEQAN_HOST_DEVICE void
_addHit(FilterDelegate<TSpec, Traits> & me, TFinder const & finder, HammingDistance)
{
    typedef typename Value<typename Traits::THits>::Type    THit;

    THit hit = { range(textIterator(finder)), finder._patternIt, getScore(finder) };

    appendValue(me.hits, hit, Insist(), typename Traits::TThreading());
}

#endif  // #ifndef APP_YARA_MAPPER_FILTER_H_
