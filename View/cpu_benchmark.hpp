#include "View.h"
#include "hof.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>

using Time = decltype(std::chrono::high_resolution_clock::now());
auto ms(Time const& t0, Time const& t1) { return std::chrono::duration_cast<std::chrono::microseconds>(std::max(t1, t0) - std::min(t1, t0)).count() / 1000.0; }

template<typename T>
bool is_same(std::vector<T> const& u, std::vector<T> const& v)
{
    if (u.size() != v.size()) { return false; }
    size_t i = 0;
    auto err = (T)0;
    for (; i < u.size(); ++i)
    {
        err = abs((u[i] - v[i]) / (u[i] + v[i]));
        if (err > 5e-6)
        {
            std::cout << "mismatch " << i << "  " << err << "\n";
            return false;
        }
    }
    return true;
}

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

double uniform_rnd(long& state)
{
    const long A = 48271;      /* multiplier*/
    const long M = 2147483647; /* modulus */
    const long Q = M / A;      /* quotient */
    const long R = M % A;      /* remainder */
    long t = A * (state % Q) - R * (state / Q);
    if (t > 0) { state = t; }
    else { state = t + M; }
    return ((double)state / M);
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

    auto ref = naive(A, B), func_ref = functional_naive<n>(A, B);

    auto summary = [&](std::string const& title, std::vector<R> const& v)
    {
        std::cout << title << ": ";
        for (auto const& r : v) { std::cout << r.second << " ms " "(" << (is_same(r.first, ref.first) ? '+' : '-') << ") "; }
        std::cout << "\n";
    };

    summary("CPU Blocked 4", { ref, func_ref, blocked<4>(A, B), view_blocked<n, 4>(A, B), functional_view_blocked<n, 4>(A, B) });
    summary("CPU Blocked 8", { ref, func_ref, blocked<8>(A, B), view_blocked<n, 8>(A, B), functional_view_blocked<n, 8>(A, B) });
    summary("CPU Blocked 16", { ref, func_ref, blocked<16>(A, B), view_blocked<n, 16>(A, B), functional_view_blocked<n, 16>(A, B) });
    summary("CPU Blocked 32", { ref, func_ref, blocked<32>(A, B), view_blocked<n, 32>(A, B), functional_view_blocked<n, 32>(A, B) });
}

void cpu_benchmark() {
    std::cout << "               naive       func naive   blocked   view blocked   func view blocked\n";
    invoke<64>();
    invoke<128>();
    invoke<256>();
    invoke<512>();
    invoke<1024>();
}