#include "benchmark_util.h"
#include "hof.h"
#include "View.h"
#include <chrono>
#include <utility>
#include <vector>

using namespace hof;

template<typename F>
auto lift(F f) {
    return [&](auto res, auto A, auto B) {
        zip(res, f, A, B);
    };
}

namespace t1 {
    template<typename T0, typename T1, typename T2, typename T3>
    void kernel(T0&& res, T1&& A2, T2&& B2, T3&& t00)
    {
        map(res, [&](auto res, auto A1)
        {
            map(res, [&](auto res, auto B1)
            {
                rnz(res, t00, add, mul, A1, B1);
            }, B2);
        }, A2);
    }

    template<int n, typename T>
    auto test(std::vector<T> const& A, std::vector<T> const& B)
    {
        std::vector<T> C(n*n);
        View<T const*, T, to_list_t<P<n, n>, P<n, 1>>> vA(A.data());
        View<T const*, T, to_list_t<P<n, 1>, P<n, n>>> vB(B.data());
        View<T*, T, to_list_t<P<n, n>, P<n, 1>>> vC(C.data());
        View<T*, T, EmptyList> tmp(new T);
        auto t0 = std::chrono::high_resolution_clock::now();
        kernel(vC, vA, vB, tmp);
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::make_pair(C, ms(t0, t1));;
    }
}

namespace t2 {
    // Subdividing the rnz :
    template<typename T0, typename T1, typename T2, typename T3, typename T4>
    void kernel(T0&& res, T1&& A3, T2&& B3, T3&& t00, T4&& t10)
    {
        map(res, [&](auto res, auto A2)
        {
            map(res, [&](auto res, auto B2)
            {
                rnz(res, t10, add3, [&](auto res, auto A1, auto B1)
                {
                    rnz(res, t00, add, mul, A1, B1);
                }, A2, B2);
            }, B3);
        }, A3);
    }

    template<int n, int b, typename T>
    auto test(std::vector<T> const& A, std::vector<T> const& B)
    {
        std::vector<T> C(n*n);
        View<T const*, T, subdiv_t<1, b, to_list_t<P<n, n>, P<n, 1>>>> vA(A.data());
        View<T const*, T, subdiv_t<1, b, to_list_t<P<n, 1>, P<n, n>>>> vB(B.data());
        View<T*, T, to_list_t<P<n, n>, P<n, 1>>> vC(C.data());
        View<T*, T, EmptyList> tmp1(new T), tmp2(new T);
        auto t0 = std::chrono::high_resolution_clock::now();
        kernel(vC, vA, vB, tmp1, tmp2);
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::make_pair(C, ms(t0, t1));;
    }
}

namespace t3 {
    // Raising the upper rnz :
    template<typename T0, typename T1, typename T2, typename T3, typename T4>
    void kernel(T0&& res, T1&& A3, T2&& B3, T3&& t00, T4&& t11)
    {
        map(res, [&](auto res, auto A2)
        {
            rnz(res, t11, lift(add), [&](auto res, auto A1, auto B2)
            {
                map(res, [&](auto res, auto B1)
                {
                    rnz(res, t00, add, mul, A1, B1);
                }, B2);
            }, A2, B3);
        }, A3);
    }

    template<int n, int b, typename T>
    auto test(std::vector<T> const& A, std::vector<T> const& B)
    {
        std::vector<T> C(n*n), tmp(n);
        View<T const*, T, subdiv_t<1, b, to_list_t<P<n, n>, P<n, 1>>>> vA(A.data());
        View<T const*, T, flip_t<0, subdiv_t<1, b, to_list_t<P<n, 1>, P<n, n>>>>> vB(B.data());
        View<T*, T, to_list_t<P<n, n>, P<n, 1>>> vC(C.data());
        View<T*, T, EmptyList> tmp1(new T);
        View<T*, T, to_list_t<P<n, 1>>> tmp2(tmp.data());
        auto t0 = std::chrono::high_resolution_clock::now();
        kernel(vC, vA, vB, tmp1, tmp2);
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::make_pair(C, ms(t0, t1));;
    }
}

namespace t4 {
    // Raising the upper rnz :
    template<int b, typename T0, typename T1, typename T2, typename T3, typename T4>
    void kernel(T0&& res, T1&& _A3, T2&& _B3, T3&& t00, T4&& t12)
    {
        auto A3 = flip<0>(subdiv<1, b>(_A3));
        auto B3 = flip<0>(subdiv<1, b>(_B3));
        rnz(res, t12, lift(lift(add)), [&](auto res, auto A2, auto B2)
        {
            map(res, [&](auto res, auto A1)
            {
                map(res, [&](auto res, auto B1)
                {
                    rnz(res, t00, add, mul, A1, B1);
                }, B2);
            }, A2);
        }, A3, B3);
    }

    template<int n, int b, typename T>
    auto test(std::vector<T> const& A, std::vector<T> const& B)
    {
        std::vector<T> C(n*n), tmp(n*n);
        View<T const*, T, to_list_t<P<n, n>, P<n, 1>>> vA(A.data());
        View<T const*, T, to_list_t<P<n, 1>, P<n, n>>> vB(B.data());
        View<T*, T, to_list_t<P<n, n>, P<n, 1>>> vC(C.data());
        View<T*, T, EmptyList> tmp1(new T);
        View<T*, T, to_list_t<P<n, n>, P<n, 1>>> tmp2(tmp.data());
        auto t0 = std::chrono::high_resolution_clock::now();
        kernel<b>(vC, vA, vB, tmp1, tmp2);
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::make_pair(C, ms(t0, t1));;
    }
}