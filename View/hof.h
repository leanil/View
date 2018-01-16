#pragma once

#include "View.h"

namespace hof {
    template<typename F, typename PA, typename A, typename DA, typename PB, typename B, typename DB>
    void map(View<PB, B, DB> result, F f, View<PA, A, DA> v) {
        static_assert(DA::head::dim == DB::head::dim, "input and output sizes must be equal");
        for (size_t i = 0; i < DA::head::dim; ++i) {
            f(result[i], v[i]);
        }
    }

    template<typename F, typename PA, typename A, typename DA, typename PB, typename B, typename DB>
    void reduce(View<PB, B, DB> result, F f, View<PA, A, DA> v) {
        result = v[0];
        for (size_t i = 1; i < DA::head::dim; ++i) {
            f(result, v[i], result);
        }
    }

    template<typename R, typename F, typename P, typename A, typename D, typename... V>
    void zip(R result, F f, View<P, A, D> v1, V... v) {
        for (size_t i = 0; i < D::head::dim; ++i) {
            f(result[i], v1[i], v[i]...);
        }
    }

    template<typename R, typename T, typename F, typename G, typename P, typename A, typename D, typename... V>
    void rnz(R result, T tmp, F f, G g, View<P, A, D> v1, V... v) {
        g(result, v1[0], v[0]...);
        for (size_t i = 1; i < D::head::dim; ++i) {
            g(tmp, v1[i], v[i]...);
            f(result, result, tmp);
        }
    }

    auto add = [](auto result, auto a, auto b) {result = a + b; };
    auto mul = [](auto result, auto a, auto b) {result = a * b; };
}