#pragma once

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

#ifdef HOF_SPEC
    // TODO: specialize the other hof-s as well
    // TODO: implement n-ary versions
    template<typename F, typename Ptr, typename A, size_t D, size_t S1, size_t S2, size_t S3>
    INLINE void zip(View<Ptr, A, List<P<D, S3>, EmptyList>> HOF_REF result, F HOF_REF f, View<Ptr, A, List<P<D, S1>, EmptyList>> HOF_REF v1, View<Ptr, A, List<P<D, S2>, EmptyList>> HOF_REF v2) {
        auto p1 = v1.data;
        auto p2 = v2.data;
        auto p3 = result.data;
#ifndef VIEW_PTR_ONLY
        p1 += v1.base;
        p2 += v2.base;
        p3 += result.base;
#endif
        auto end = p1 + D*S1;
        for (; p1 < end; p1 += S1, p2 += S2, p3 += S3) {
            f(*p3, *p1, *p2);
        }
    }

    template<typename R, typename T, typename F, typename G, typename Ptr, typename A, size_t D, size_t S1, size_t S2>
    INLINE void rnz(R HOF_REF result, T HOF_REF tmp, F HOF_REF f, G HOF_REF g, View<Ptr, A, List<P<D, S1>, EmptyList>> HOF_REF v1, View<Ptr, A, List<P<D, S2>, EmptyList>> HOF_REF v2) {
        auto p1 = v1.data;
        auto p2 = v2.data;
#ifndef VIEW_PTR_ONLY
        p1 += v1.base;
        p2 += v2.base;
#endif
        auto end = p1 + D*S1;
        g(result, *p1, *p2);
        for (p1 += S1, p2 += S2; p1 < end; p1 += S1, p2 += S2) {
            g(tmp, *p1, *p2);
            f(result, result, tmp);
        }
    }

    auto add = [](auto& result, auto LAMBDA_REF a, auto LAMBDA_REF b) {result = a + b; };
    auto mul = [](auto& result, auto LAMBDA_REF a, auto LAMBDA_REF b) {result = a * b; };
#else
    auto add = [](auto LAMBDA_REF result, auto LAMBDA_REF a, auto LAMBDA_REF b) {result = a + b; };
    auto mul = [](auto LAMBDA_REF result, auto LAMBDA_REF a, auto LAMBDA_REF b) {result = a * b; };
#endif
}