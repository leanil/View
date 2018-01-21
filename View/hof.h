#pragma once

#include "View.h"

namespace hof {
    //TODO: how to pass the lambdas
    template<typename F, typename PA, typename A, typename DA, typename PB, typename B, typename DB>
    INLINE void map(View<PB, B, DB> HOF_REF result, F HOF_REF f, View<PA, A, DA> HOF_REF v) {
        static_assert(DA::head::dim == DB::head::dim, "input and output sizes must be equal");
        for (size_t i = 0; i < DA::head::dim; ++i) {
            f(result[i], v[i]);
        }
    }

    template<typename F, typename PA, typename A, typename DA, typename PB, typename B, typename DB>
    INLINE void reduce(View<PB, B, DB> HOF_REF result, F HOF_REF f, View<PA, A, DA> HOF_REF v) {
        result = v[0];
        for (size_t i = 1; i < DA::head::dim; ++i) {
            f(result, v[i], result);
        }
    }

    template<typename R, typename F, typename P, typename A, typename D, typename... V>
    INLINE void zip(R HOF_REF result, F HOF_REF f, View<P, A, D> HOF_REF v1, V HOF_REF... v) {
        for (size_t i = 0; i < D::head::dim; ++i) {
            f(result[i], v1[i], v[i]...);
        }
    }

    template<typename R, typename T, typename F, typename G, typename P, typename A, typename D, typename... V>
    INLINE void rnz(R HOF_REF result, T HOF_REF tmp, F HOF_REF f, G HOF_REF g, View<P, A, D> HOF_REF v1, V HOF_REF... v) {
        g(result, v1[0], v[0]...);
        for (size_t i = 1; i < D::head::dim; ++i) {
            g(tmp, v1[i], v[i]...);
            f(result, result, tmp);
        }
    }

    auto add = [](auto LAMBDA_REF result, auto LAMBDA_REF a, auto LAMBDA_REF b) {result = a + b; };
    auto mul = [](auto LAMBDA_REF result, auto LAMBDA_REF a, auto LAMBDA_REF b) {result = a * b; };
}