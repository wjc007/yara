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
// Class Locator
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec = void>
struct Locator
{
    typename Member<Locator, Ranges_>::Type    ranges;

    template <typename TFinder>
    inline SEQAN_HOST_DEVICE void
    operator() (TFinder const & finder)
    {
        appendRange(*this, finder);
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
struct Member<Locator<TIndex, TSpec>, Ranges_>
{
    typedef Pair<typename Size<TIndex>::Type>   TRange_;
    typedef String<TRange_>                     Type;
};

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec>
struct Member<Locator<TIndex, Device<TSpec> >, Ranges_>
{
    typedef Pair<typename Size<TIndex>::Type>   TRange_;
    typedef thrust::device_vector<TRange_>      Type;
};
#endif

template <typename TIndex, typename TSpec>
struct Member<Locator<TIndex, View<TSpec> >, Ranges_>
{
    typedef typename Member<Locator<TIndex, TSpec>, Ranges_>::Type  TRanges_;
    typedef ContainerView<TRanges_, Resizable<TSpec> >              Type;
};
}

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TIndex, typename TSpec>
struct View<Locator<TIndex, TSpec> >
{
    typedef Locator<TIndex, View<TSpec> >  Type;
};
}

// ----------------------------------------------------------------------------
// Metafunction Device
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TIndex, typename TSpec>
struct Device<Locator<TIndex, TSpec> >
{
    typedef Locator<TIndex, Device<TSpec> >  Type;
};
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TPattern>
inline void
init(Locator<TIndex, TSpec> & locator, TPattern const & pattern)
{
    reserve(locator.ranges, length(needle(pattern)), Exact());
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
typename View<Locator<TIndex, TSpec> >::Type
view(Locator<TIndex, TSpec> & locator)
{
    typename View<Locator<TIndex, TSpec> >::Type locatorView;

    locatorView.ranges = view(locator.ranges);

    return locatorView;
}

// ----------------------------------------------------------------------------
// Function appendRange()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TFinder>
inline void
appendRange(Locator<TIndex, TSpec> & locator, TFinder const & finder)
{
    SEQAN_OMP_PRAGMA(critical(Locator_appendRange))
    appendValue(locator.ranges, range(textIterator(finder)));
}

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec, typename TFinder>
inline SEQAN_HOST_DEVICE void
appendRange(Locator<TIndex, View<Device<TSpec> > > & /* locator */, TFinder const & /* finder */)
{
    // TODO(esiragusa): Global lock.
//    appendValue(locator.ranges, range(textIterator(finder)));
}
#endif


#endif  // #ifndef APP_CUDAMAPPER_LOCATOR_H_
