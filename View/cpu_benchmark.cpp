#define LOOP_REF const&
#include "View.h"
#include "hof.h"
#include <random>

// Standard C++ includes
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
    View<T const*, T, to_list_t<P<Bs, b*n>, P<Bs, b>, P<b, n>, P<b, 1>>> vA(A.data()), vB(B.data());
    View<T*, T, to_list_t<P<Bs, b*n>, P<Bs, b>, P<b, n>, P<b, 1>>> vC(C.data());
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int bi = 0; bi < Bs; ++bi)
    { //block index 1
        for (int bj = 0; bj < Bs; ++bj)
        { //block index 2
            for (int bk = 0; bk < Bs; ++bk)
            { //block index 3
                auto bA = vA[bi][bk], bB = vB[bk][bj];
                auto bC = vC[bi][bj];
                for (int i = 0; i < b; ++i)
                {
                    for (int j = 0; j < b; ++j)
                    {
                        T sum = T();
                        for (int k = 0; k < b; ++k)
                        {
                            sum += bA[i][k] * bB[k][j];
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
    View<T const*, T, to_list_t<P<n, n>, P<n, 1>>> vA(A.data());
    View<T const*, T, to_list_t<P<n, 1>, P<n, n>>> vB(B.data());
    View<T*, T, to_list_t<P<n, n>, P<n, 1>>> vC(C.data());
    T t;
    View<T*, T, EmptyList> tmp(&t);
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
    std::vector<T> C(n*n);
    constexpr auto Bs = n / b;
    View<T const*, T, to_list_t<P<Bs, b*n>, P<Bs, b>, P<b, n>, P<b, 1>>> vA(A.data());
    View<T const*, T, to_list_t<P<Bs, b>, P<Bs, b*n>, P<b, 1>, P<b, n>>> vB(B.data());
    View<T*, T, to_list_t<P<Bs, b*n>, P<Bs, b>, P<b, n>, P<b, 1>>> vC(C.data());
    std::vector<T> D(b*b);
    View<T*, T, to_list_t<P<b, b>, P<b, 1>>> vD(D.data());
    T t;
    View<T*, T, EmptyList> tmp(&t);
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

template<int n, int b, typename T>
auto functional_loop_convert(std::vector<T> const& A, std::vector<T> const& B)
{
    std::vector<T> C(n*n);
    constexpr auto Bs = n / b;
    View<T const*, T, to_list_t<P<Bs, b*n>, P<Bs, b>, P<b, n>, P<b, 1>>> vA(A.data());
    View<T const*, T, to_list_t<P<Bs, b>, P<Bs, b*n>, P<b, 1>, P<b, n>>> vB(B.data());
    View<T*, T, to_list_t<P<Bs, b*n>, P<Bs, b>, P<b, n>, P<b, 1>>> vC(C.data());
    std::vector<T> D(b*b);
    View<T*, T, to_list_t<P<b, b>, P<b, 1>>> vD(D.data());
    T t;
    View<T*, T, EmptyList> tmp(&t);
    auto t0 = std::chrono::high_resolution_clock::now();

    auto LOOP_REF result0 = vC;
    auto LOOP_REF v0 = vA;
    // map(vC, [&](auto result, auto br), vA)
    for (size_t i = 0; i < v0.size; ++i) {
        auto LOOP_REF result1 = result0[i];
        auto LOOP_REF br1 = v0[i];
        // [&](auto result, auto br)
        auto LOOP_REF result2 = result1;
        auto LOOP_REF v2 = vB;
        // map(result, [&](auto result, auto bc), vB)
        for (size_t i = 0; i < v2.size; ++i) {
            auto LOOP_REF result3 = result2[i];
            auto LOOP_REF bc3 = v2[i];
            // [&](auto result, auto bc)
            auto LOOP_REF result4 = result3;
            auto LOOP_REF tmp4 = vD;
            auto LOOP_REF v14 = br1;
            auto LOOP_REF v24 = bc3;
            // rnz(result, vD, [&](auto result, auto b1, auto b2), [&](auto result, auto bA, auto bB), br, bc)
            auto LOOP_REF result5 = result4;
            auto LOOP_REF bA5 = v14[0];
            auto LOOP_REF bB5 = v24[0];
            // [&](auto result, auto bA, auto bB)
            auto LOOP_REF result6 = result5;
            auto LOOP_REF v6 = bA5;
            // map(result, [&](auto result, auto r), bA)
            for (size_t i = 0; i < v6.size; ++i) {
                auto LOOP_REF result7 = result6[i];
                auto LOOP_REF r7 = v6[i];
                // [&](auto result, auto r)
                auto LOOP_REF result8 = result7;
                auto LOOP_REF v8 = bB5;
                // map(result, [&](auto result, auto c), bB)
                for (size_t i = 0; i < v8.size; ++i) {
                    auto LOOP_REF result9 = result8[i];
                    auto LOOP_REF c9 = v8[i];
                    // [&](auto result, auto c)
                    auto LOOP_REF result10 = result9;
                    auto LOOP_REF tmp10 = tmp;
                    auto LOOP_REF v110 = r7;
                    auto LOOP_REF v210 = c9;
                    // rnz(result, tmp, add, mul, r, c)
                    auto LOOP_REF result11 = result10;
                    auto LOOP_REF a11 = v110[0];
                    auto LOOP_REF b11 = v210[0];
                    // mul
                    result11 = a11 * b11;
                    for (size_t i = 1; i < v110.size; ++i) {
                        auto LOOP_REF result11 = tmp10;
                        auto LOOP_REF a11 = v110[i];
                        auto LOOP_REF b11 = v210[i];
                        // mul
                        result11 = a11 * b11;
                        auto LOOP_REF result12 = result10;
                        auto LOOP_REF a12 = result10;
                        auto LOOP_REF b12 = tmp10;
                        // add
                        result12 = a12 + b12;
                    }
                }
            }
            for (size_t i = 1; i < v14.size; ++i) {
                auto LOOP_REF result5 = tmp4;
                auto LOOP_REF bA5 = v14[i];
                auto LOOP_REF bB5 = v24[i];
                // [&](auto result, auto bA, auto bB)
                auto LOOP_REF result6 = result5;
                auto LOOP_REF v6 = bA5;
                // map(result, [&](auto result, auto r), bA)
                for (size_t i = 0; i < v6.size; ++i) {
                    auto LOOP_REF result7 = result6[i];
                    auto LOOP_REF r7 = v6[i];
                    // [&](auto result, auto r)
                    auto LOOP_REF result8 = result7;
                    auto LOOP_REF v8 = bB5;
                    // map(result, [&](auto result, auto c), bB)
                    for (size_t i = 0; i < v8.size; ++i) {
                        auto LOOP_REF result9 = result8[i];
                        auto LOOP_REF c9 = v8[i];
                        // [&](auto result, auto c)
                        auto LOOP_REF result10 = result9;
                        auto LOOP_REF tmp10 = tmp;
                        auto LOOP_REF v110 = r7;
                        auto LOOP_REF v210 = c9;
                        // rnz(result, tmp, add, mul, r, c)
                        auto LOOP_REF result11 = result10;
                        auto LOOP_REF a11 = v110[0];
                        auto LOOP_REF b11 = v210[0];
                        // mul
                        result11 = a11 * b11;
                        for (size_t i = 1; i < v110.size; ++i) {
                            auto LOOP_REF result11 = tmp10;
                            auto LOOP_REF a11 = v110[i];
                            auto LOOP_REF b11 = v210[i];
                            // mul
                            result11 = a11 * b11;
                            auto LOOP_REF result12 = result10;
                            auto LOOP_REF a12 = result10;
                            auto LOOP_REF b12 = tmp10;
                            // add
                            result12 = a12 + b12;
                        }
                    }
                }
                auto LOOP_REF result13 = result4;
                auto LOOP_REF b113 = result4;
                auto LOOP_REF b213 = tmp4;
                // [&](auto result, auto b1, auto b2)
                auto LOOP_REF result14 = result13;
                auto LOOP_REF v114 = b113;
                auto LOOP_REF v214 = b213;
                // zip(result, [&](auto result, auto r1, auto r2), b1, b2)
                for (size_t i = 0; i < v114.size; ++i) {
                    auto LOOP_REF result15 = result14[i];
                    auto LOOP_REF r115 = v114[i];
                    auto LOOP_REF r215 = v214[i];
                    // [&](auto result, auto r1, auto r2)
                    auto LOOP_REF result16 = result15;
                    auto LOOP_REF v116 = r115;
                    auto LOOP_REF v216 = r215;
                    // zip(result, add, r1, r2)
                    for (size_t i = 0; i < v116.size; ++i) {
                        auto LOOP_REF result17 = result16[i];
                        auto LOOP_REF a17 = v116[i];
                        auto LOOP_REF b17 = v216[i];
                        // add
                        result17 = a17 + b17;
                    }
                }
            }
        }
    }
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
        for (auto const& r : v) { std::cout << r.second << " ms (" << (is_same(r.first, ref.first) ? '+' : '-') << ") "; }
        std::cout << "\n";
    };

    summary("CPU Blocked 4", { ref, func_ref, blocked<4>(A, B), view_blocked<n, 4>(A, B), functional_view_blocked<n, 4>(A, B), functional_loop_convert<n, 4>(A, B) });
    summary("CPU Blocked 8", { ref, func_ref, blocked<8>(A, B), view_blocked<n, 8>(A, B), functional_view_blocked<n, 8>(A, B), functional_loop_convert<n, 8>(A, B) });
    summary("CPU Blocked 16", { ref, func_ref, blocked<16>(A, B), view_blocked<n, 16>(A, B), functional_view_blocked<n, 16>(A, B), functional_loop_convert<n, 16>(A, B) });
    summary("CPU Blocked 32", { ref, func_ref, blocked<32>(A, B), view_blocked<n, 32>(A, B), functional_view_blocked<n, 32>(A, B), functional_loop_convert<n, 32>(A, B) });
}

void test() {
    std::vector<int> d(36);
    iota(d.begin(), d.end(), 0);
    View<int*, int, to_list_t<P<2, 18>, P<2, 3>, P<3, 6>, P<3, 1>>> v(d.data());
    std::cout << v << std::endl;
    std::cout << v[1] << std::endl;
    auto vv = v(1, 1);
    std::cout << vv << std::endl;
    std::cout << v(1, 1, 2, 2) << std::endl;
    std::cout << v[1](1, 2)[2] << std::endl;
    using namespace std::placeholders;
    std::cout << v(_2, _1, _4, _3) << std::endl;
    auto v2 = v(_1, 0, _2, 2);
    std::cout << v2 << std::endl;
    std::cout << v2(_2, _1) << std::endl;
}

int main()
{
    std::cout << "               naive       func naive   blocked   view blocked   func view blocked  func loop convert\n";
    invoke<64>();
    invoke<128>();
    invoke<256>();
    invoke<512>();
    invoke<1024>();
}