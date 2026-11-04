#pragma once

// Minimal fork-join helper used by the extractor and the analysis pipelines.
//
// parallelFor(count, body) calls body(i) for every i in [0, count) using a
// short-lived pool of std::threads with dynamic (atomic counter) scheduling.
// Results must be written to per-index slots, so output order never depends
// on scheduling. Nested calls run serially on the calling worker, which keeps
// e.g. "one mesh per k in parallel" from oversubscribing the CPU when each
// mesh build would otherwise also parallelize its own grid sampling.

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace iso_parallel {

namespace detail {
inline thread_local int regionDepth = 0;

struct RegionGuard {
    RegionGuard() { ++regionDepth; }
    ~RegionGuard() { --regionDepth; }
    RegionGuard(const RegionGuard&) = delete;
    RegionGuard& operator=(const RegionGuard&) = delete;
};
} // namespace detail

inline unsigned hardwareThreads() {
    const unsigned n = std::thread::hardware_concurrency();
    return n == 0 ? 1u : n;
}

// Exceptions thrown by body are captured; the first one is rethrown on the
// calling thread after all workers have joined.
template <typename Body>
void parallelFor(std::size_t count, Body&& body, unsigned maxThreads = 0) {
    if (count == 0) return;
    unsigned threads = maxThreads == 0 ? hardwareThreads() : maxThreads;
    threads = static_cast<unsigned>(std::min<std::size_t>(threads, count));
    if (threads <= 1 || detail::regionDepth > 0) {
        detail::RegionGuard guard;
        for (std::size_t i = 0; i < count; ++i) body(i);
        return;
    }

    std::atomic<std::size_t> next{0};
    std::exception_ptr firstError;
    std::mutex errorMutex;
    auto worker = [&]() {
        detail::RegionGuard guard;
        for (;;) {
            const std::size_t i = next.fetch_add(1, std::memory_order_relaxed);
            if (i >= count) return;
            try {
                body(i);
            } catch (...) {
                std::lock_guard<std::mutex> lock(errorMutex);
                if (!firstError) firstError = std::current_exception();
                next.store(count, std::memory_order_relaxed); // stop handing out work
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(threads - 1);
    for (unsigned t = 1; t < threads; ++t) pool.emplace_back(worker);
    worker(); // the calling thread participates
    for (std::thread& thread : pool) thread.join();
    if (firstError) std::rethrow_exception(firstError);
}

} // namespace iso_parallel
