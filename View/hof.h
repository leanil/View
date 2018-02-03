#pragma once

#include "View.h"

namespace hof {
    template<typename F, typename PA, typename A, typename DA, typename PB, typename B, typename DB>
    __forceinline void map(View<PB, B, DB> result, F f, View<PA, A, DA> v) {
        static_assert(DA::head::dim == DB::head::dim, "input and output sizes must be equal");
        for (size_t i = 0; i < DA::head::dim; ++i) {
            f(result[i], v[i]);
        }
    }

    template<typename F, typename PA, typename A, typename DA, typename PB, typename B, typename DB>
    __forceinline void reduce(View<PB, B, DB> result, F f, View<PA, A, DA> v) {
        result = v[0];
        for (size_t i = 1; i < DA::head::dim; ++i) {
            f(result, v[i], result);
        }
    }

    template<typename R, typename F, typename P, typename A, typename D, typename... V>
    __forceinline void zip(R result, F f, View<P, A, D> v1, V... v) {
        for (size_t i = 0; i < D::head::dim; ++i) {
            f(result[i], v1[i], v[i]...);
        }
    }

    template<typename R, typename T, typename F, typename G, typename P, typename A, typename D, typename... V>
    __forceinline void rnz(R result, T tmp, F f, G g, View<P, A, D> v1, V... v) {
        g(result, v1[0], v[0]...);
        for (size_t i = 1; i < D::head::dim; ++i) {
            g(tmp, v1[i], v[i]...);
            f(result, result, tmp);
        }
    }

    template<typename F, typename Ptr, typename A, size_t D, size_t S1, size_t S2, size_t S3>
    __forceinline void zip(View<Ptr, A, List<P<D, S3>, EmptyList>> result, F f, View<Ptr, A, List<P<D, S1>, EmptyList>> v1, View<Ptr, A, List<P<D, S2>, EmptyList>> v2) {
        auto p1 = v1.data, p2 = v2.data, pr = result.data;
        auto end = p1 + D*S1;
        for (; p1 < end; p1 += S1, p2 += S2, pr += S3) {
            *pr = f(*p1, *p2);
        }
    }

    template<typename R, typename T, typename F, typename G, typename Ptr, typename A, size_t D, size_t S1, size_t S2>
    __forceinline void rnz(R result, T, F f, G g, View<Ptr, A, List<P<D, S1>, EmptyList>> v1, View<Ptr, A, List<P<D, S2>, EmptyList>> v2) {
        auto p1 = v1.data, p2 = v2.data;
        auto end = p1 + D*S1;
        auto r = g(*p1, *p2);
        for (p1 += S1; p1 < end; p1 += S1) {
            p2 += S2;
            auto tmp = g(*p1, *p2);
            r = f(r, tmp);
        }
        result = r;
    }

    const auto add = [](auto a, auto b) {return a + b; };
    const auto mul = [](auto a, auto b) {return a * b; };
}