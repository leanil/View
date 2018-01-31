#pragma once

#include "SYCL/sycl.hpp"
#include <algorithm>
#include <chrono>

using Time = decltype(std::chrono::high_resolution_clock::now());

auto ms(Time const& t0, Time const& t1) { return std::chrono::duration_cast<std::chrono::microseconds>(std::max(t1, t0) - std::min(t1, t0)).count() / 1000.0; }

auto ms(cl::sycl::event const& e) {
    cl::sycl::cl_ulong t0, t1;
    auto status1 = clGetEventProfilingInfo(e.get(), CL_PROFILING_COMMAND_START, sizeof(cl::sycl::cl_ulong), &t0, nullptr);
    auto status2 = clGetEventProfilingInfo(e.get(), CL_PROFILING_COMMAND_END, sizeof(cl::sycl::cl_ulong), &t1, nullptr);
    return status1 == CL_SUCCESS && status2 == CL_SUCCESS ? (t1 - t0) / 1'000'000.0 : -1.0;
}

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
