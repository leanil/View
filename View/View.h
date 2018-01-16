#pragma once

#include "List.h"
#include <functional>
#include <iostream>
#include <type_traits>
#include <utility>

template<size_t D, size_t S>
struct P {
    static constexpr size_t dim = D;
    static constexpr size_t stride = S;
};

template<typename Ptr, typename T, typename Ds> class View;

template<typename Ptr, typename T, typename D, typename Ds>
class View<Ptr, T, List<D, Ds>> {
public:
    View(const Ptr& data, size_t base = 0) : data(data), base(base) {}

    View<Ptr, T, List<D, Ds>>& operator=(const View<Ptr, T, List<D, Ds>>& other) const {
        for (size_t i = 0; i < D::dim; ++i) {
            (*this)[i] = other[i];
        }
        return *this;
    }

    auto operator[](size_t idx) const {
        return View<Ptr, T, Ds>(data, base + idx*D::stride);
    }

    template<typename I, typename... Idx>
    auto operator()(I i, Idx... idx) const {
        return slice<EmptyList, List<D, Ds>>(0, i, idx...);
    }

private:
    template<typename SortedDims, typename Dims, typename I, typename... Idx, std::enable_if_t<std::is_placeholder_v<I> != 0, int> = 0>
    auto slice(size_t offset, I i, Idx... idx) const {
        return slice<typename Insert<SortedDims, KeyValuePair<Int<std::is_placeholder_v<I>>, typename Dims::head>>::type, typename Dims::tail>(offset, idx...);
    }

    template<typename SortedDims, typename Dims, typename I, typename... Idx, std::enable_if_t<std::is_placeholder_v<I> == 0, int> = 0>
    auto slice(size_t offset, I i, Idx... idx) const {
        return slice<SortedDims, typename Dims::tail>(offset + i * Dims::head::stride, idx...);
    }

    template<typename SortedDims, typename Dims>
    auto slice(size_t offset) const {
        return View<Ptr, T, typename Concat<typename Map<SortedDims, get_value>::type, Dims>::type>(data, base + offset);
    }

    Ptr data;
    size_t base;
};

template<typename Ptr, typename T, typename Dims>
std::ostream& operator<<(std::ostream& out, const View<Ptr, T, Dims>& v) {
    out << "{";
    for (size_t i = 0; i < Dims::head::dim; ++i) {
        out << (i ? "," : "") << v[i];
    }
    return out << "}";
}

template<typename Ptr, typename T>
class View<Ptr, T, EmptyList> {
public:
    View(Ptr data, size_t base = 0) : data(data), base(base) {}

    View<Ptr, T, EmptyList>& operator=(const View<Ptr, T, EmptyList>& other) const {
        data[base] = other.data[other.base];
        return *this;
    }

    void operator=(T x) const {
        data[base] = x;
    }

    operator T() const {
        return data[base];
    }
private:
    Ptr data;
    size_t base;
};

template<typename Ptr, typename T>
std::ostream& operator<<(std::ostream& out, const View<Ptr, T, EmptyList>& v) {
    return out << (T)v;
}
