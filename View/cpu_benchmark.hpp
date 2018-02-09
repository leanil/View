#pragma once

#include "View.h"
#include "hof.h"
#include "benchmark_util.h"
#include "gen_test.hpp"
#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>

namespace cpu {

    template<typename T>
    auto naive(std::vector<T> const& A, std::vector<T> const& B)
    {
        auto n = (int)sqrt(A.size());
        std::vector<T> C(n*n);
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                T sum = (T)0;
                for (int k = 0; k < n; ++k)
                {
                    sum += A[i*n + k] * B[k*n + j];
                }
                C[i*n + j] = sum;
            }
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::make_pair(C, ms(t0, t1));;
    }

    template<int b, typename T>
    auto blocked(std::vector<T> const& A, std::vector<T> const& B)
    {
        auto n = (int)sqrt(A.size());
        std::vector<T> C(n*n);
        auto Bs = n / b;
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int bi = 0; bi < Bs; ++bi)
        { //block index 1
            for (int bj = 0; bj < Bs; ++bj)
            { //block index 2
                for (int bk = 0; bk < Bs; ++bk)
                { //block index 3
                    auto i0 = bi * b; auto j0 = bj * b; auto k0 = bk * b;
                    for (int i = 0; i < b; ++i)
                    {
                        auto ii = i0 + i;
                        for (int j = 0; j < b; ++j)
                        {
                            auto jj = j0 + j;
                            double sum = 0.0;
                            for (int k = 0; k < b; ++k)
                            {
                                sum += A[ii*n + k0 + k] * B[(k0 + k)*n + jj];
                            }
                            C[ii*n + jj] += sum;
                        }
                    }
                }
            }
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::make_pair(C, ms(t0, t1));;
    }

    template<int n, int b, typename T>
    auto view_blocked(std::vector<T> const& A, std::vector<T> const& B)
    {
        std::vector<T> C(n*n);
        constexpr auto Bs = n / b;
        auto Adata = A.data(), Bdata = B.data();
        auto Cdata = C.data();
        View<T const*, T, to_list_t<P<Bs, b*n>, P<Bs, b>, P<b, n>, P<b, 1>>> vA(Adata);
        View<T const*, T, to_list_t<P<Bs, b>, P<Bs, b*n>, P<b, 1>, P<b, n>>> vB(Bdata);
        View<T*, T, to_list_t<P<Bs, b*n>, P<Bs, b>, P<b, n>, P<b, 1>>> vC(Cdata);
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int bi = 0; bi < Bs; ++bi)
        { //block index 1
            for (int bj = 0; bj < Bs; ++bj)
            { //block index 2
                for (int bk = 0; bk < Bs; ++bk)
                { //block index 3
                    auto bA = vA[bi][bk];
                    auto bB = vB[bj][bk];
                    auto bC = vC[bi][bj];
                    for (int i = 0; i < b; ++i)
                    {
                        for (int j = 0; j < b; ++j)
                        {
                            T sum = T();
                            for (int k = 0; k < b; ++k)
                            {
                                sum += bA[i][k] * bB[j][k];
                            }
                            bC[i][j] = bC[i][j] + sum;
                        }
                    }
                }
            }
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::make_pair(C, ms(t0, t1));;
    }

    using namespace hof;

    template<int n, typename T>
    auto functional_naive(std::vector<T> const& A, std::vector<T> const& B)
    {
        std::vector<T> C(n*n);
        auto Adata = A.data(), Bdata = B.data();
        auto Cdata = C.data();
        View<T const*, T, to_list_t<P<n, n>, P<n, 1>>> vA(Adata);
        View<T const*, T, to_list_t<P<n, 1>, P<n, n>>> vB(Bdata);
        View<T*, T, to_list_t<P<n, n>, P<n, 1>>> vC(Cdata);
        T* t = new T;
        View<T*, T, EmptyList> tmp(t);
        auto t0 = std::chrono::high_resolution_clock::now();
        map(vC,
            [&](auto result, auto r) {
            map(result,
                [&](auto result, auto c) {
                rnz(result, tmp, add, mul, r, c); },
                vB); },
            vA);
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::make_pair(C, ms(t0, t1));;
    }

    template<int n, typename T>
    auto functional_naive_exchange(std::vector<T> const& A, std::vector<T> const& B)
    {
        std::vector<T> C(n*n);
        auto Adata = A.data(), Bdata = B.data();
        auto Cdata = C.data();
        View<T const*, T, to_list_t<P<n, n>, P<n, 1>>> vA(Adata);
        View<T const*, T, to_list_t<P<n, n>, P<n, 1>>> vB(Bdata);
        View<T*, T, to_list_t<P<n, n>, P<n, 1>>> vC(Cdata);
        std::vector<T> t(n);
        View<T*, T, to_list_t<P<n, 1>>> tmp(t.data());
        auto t0 = std::chrono::high_resolution_clock::now();
        map(vC,
            [&](auto result, auto r) {
            rnz(result, tmp,
                [&](auto result, auto v1, auto v2) {
                zip(result, add, v1, v2); },
                [&](auto result, auto q, auto r2) {
                    map(result,
                        [&](auto result, auto x) { result = mul(q, x); },
                        r2); },
                    r, vB); },
            vA);
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::make_pair(C, ms(t0, t1));;
    }

    template<int n, int b, typename T>
    auto functional_view_blocked(std::vector<T> const& A, std::vector<T> const& B)
    {
        constexpr auto Bs = n / b;
        std::vector<T> C(n*n), D(b*b);
        auto Adata = A.data(), Bdata = B.data();
        auto Cdata = C.data(), Ddata = D.data();
        View<T const*, T, to_list_t<P<Bs, b*n>, P<Bs, b>, P<b, n>, P<b, 1>>> vA(Adata);
        View<T const*, T, to_list_t<P<Bs, b>, P<Bs, b*n>, P<b, 1>, P<b, n>>> vB(Bdata);
        View<T*, T, to_list_t<P<Bs, b*n>, P<Bs, b>, P<b, n>, P<b, 1>>> vC(Cdata);
        View<T*, T, to_list_t<P<b, b>, P<b, 1>>> vD(Ddata);
        T q;
        T* t = &q;
        View<T*, T, EmptyList> tmp(t);
        auto t0 = std::chrono::high_resolution_clock::now();
        map(vC, [&](auto result, auto br) {
            map(result, [&](auto result, auto bc) {
                rnz(result, vD,
                    [&](auto result, auto b1, auto b2) {
                    zip(result, [&](auto result, auto r1, auto r2) {
                        zip(result, add, r1, r2); },
                        b1, b2); },
                    [&](auto result, auto bA, auto bB) {
                            map(result, [&](auto result, auto r) {
                                map(result, [&](auto result, auto c) {
                                    rnz(result, tmp, add, mul, r, c); },
                                    bB); },
                                bA); },
                            br, bc); },
                vB); },
            vA);
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::make_pair(C, ms(t0, t1));;
    }

    template<int n>
    void invoke() {
        using T = double;
        using R = std::pair<std::vector<T>, double>;

        std::cout << "\nn = " << n << "\n\n";
        // Host vectors
        std::vector<T> A(n*n);
        std::vector<T> B(n*n);

        // Initialize vectors on host
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                A[i*n + j] = (T)((i + j + 1) / (1.*n*n));//(T)uniform_rnd(state);//(i+j+1) / (1.*n*n);
                B[i*n + j] = (T)((i - j + 2) / (1.*n*n));//(T)uniform_rnd(state);//(i-j+2) / (1.*n*n);
            }
        }

        auto ref = functional_view_blocked<n, 16>(A, B);

        auto summary = [&](std::string const& title, std::vector<R> const& v)
        {
            std::cout << title << ": ";
            for (auto const& r : v) { std::cout << r.second << " ms " "(" << (is_same(r.first, ref.first) ? '+' : '-') << ") "; }
            std::cout << "\n";
        };

        summary("CPU Blocked 4", { functional_view_blocked<n, 4>(A, B), t1::test<n>(A,B), t2::test<n,4>(A,B), t3::test<n,4>(A,B), t4::test<n,4>(A,B) });
        summary("CPU Blocked 8", { functional_view_blocked<n, 8>(A, B), t1::test<n>(A,B), t2::test<n,8>(A,B), t3::test<n,8>(A,B), t4::test<n,8>(A,B) });
        summary("CPU Blocked 16", { functional_view_blocked<n, 16>(A, B), t1::test<n>(A,B), t2::test<n,16>(A,B), t3::test<n,16>(A,B), t4::test<n,16>(A,B) });
        summary("CPU Blocked 32", { functional_view_blocked<n, 32>(A, B), t1::test<n>(A,B), t2::test<n,32>(A,B), t3::test<n,32>(A,B), t4::test<n,32>(A,B) });
    }

    void benchmark() {
        std::cout << "           func_view_blocked     t1          t2          t3              t4\n";
        invoke<64>();
        invoke<128>();
        invoke<256>();
        invoke<512>();
        invoke<1024>();
    }

}
