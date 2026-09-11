#include <sycl/sycl.hpp>
#include "ParallelHDKP1.h"
#include "seq/HDKP1.h"
#include <iostream>
#include <vector>
#include <future>
#include <thread>
#include <cmath>
#include <limits>
#include <algorithm>

void warmupSYCL(sycl::queue &q) {
    uint32_t dummy = 0;
    sycl::buffer<uint32_t, 1> dummy_buf(&dummy, sycl::range<1>(1));
    q.submit([&](sycl::handler &h) {
        auto acc = dummy_buf.get_access<sycl::access::mode::read_write>(h);
        h.single_task([=]() { acc[0] += 1; });
    }).wait();
}

/**
 * Naive Embarrassing Parallel HDKP1:
 * - Swaps the outer loop over points p in K_A for a parallel loop across threads.
 * - Each worker maintains a private pNN cache for Pruning Rule 1 (dist(p, pNN) <= cmax).
 * - Calls the original, untouched sequential nnMax() from seq/NNMAX.cpp as a black box.
 * - Shared lower bound cmax is monotonically updated via sycl::atomic_ref::fetch_max.
 */
double parallelHDKP1(const std::vector<Point> &pointsA, MREP2 *B, sycl::queue *pQueue, unsigned int num_threads) {
    size_t N = pointsA.size();
    if (N == 0) return 0.0;

    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 4;
    }

    sycl::queue local_q;
    if (!pQueue) {
        try {
            local_q = sycl::queue(sycl::cpu_selector_v);
        } catch (...) {
            local_q = sycl::queue(sycl::default_selector_v);
        }
        pQueue = &local_q;
    }
    sycl::queue &q = *pQueue;

    // Shared lower bound cmax in USM memory
    auto usm_deleter = [q](double *ptr) { sycl::free(ptr, q); };
    std::shared_ptr<double> cmax_ptr(sycl::malloc_shared<double>(1, q), usm_deleter);
    *cmax_ptr = 0.0;

    // Static range partition: worker t processes points in [start, end)
    std::vector<std::future<void>> workers;
    workers.reserve(num_threads);

    size_t chunk_size = (N + num_threads - 1) / num_threads;

    for (unsigned int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, N);
        if (start >= N) break;

        workers.push_back(std::async(std::launch::async, [&pointsA, B, start, end, raw_cmax = cmax_ptr.get()]() {
            auto ref_cmax = sycl::atomic_ref<double, sycl::memory_order::relaxed, sycl::memory_scope::system, sycl::access::address_space::global_space>(*raw_cmax);

            Point pNN;
            pNN.setX(std::numeric_limits<uint32_t>::max() / 2);
            pNN.setY(std::numeric_limits<uint32_t>::max() / 2);

            for (size_t i = start; i < end; ++i) {
                Point p = pointsA[i];
                double minNN = distSinSqrt(p, pNN);
                double current_cmax = ref_cmax.load();

                // Pruning Rule 1
                if (minNN > current_cmax) {
                    // Calls original untouched sequential nnMax from seq/NNMAX.cpp
                    double res = nnMax(B, p, current_cmax, minNN, pNN);
                    if (res > current_cmax) {
                        ref_cmax.fetch_max(res);
                    }
                }
            }
        }));
    }

    for (auto &w : workers) {
        w.get();
    }

    return std::sqrt(*cmax_ptr);
}

double parallelHDKP1(MREP2 *A, MREP2 *B, sycl::queue *q, unsigned int num_threads) {
    if (!A || !B) return 0.0;
    // Extract points from k2-tree A in Z-order
    std::vector<Point> pointsA = extractPointK2tree(A, 0);
    return parallelHDKP1(pointsA, B, q, num_threads);
}

double symmetricParallelHDKP1(MREP2 *A, MREP2 *B, sycl::queue *q, unsigned int num_threads) {
    double hAB = parallelHDKP1(A, B, q, num_threads);
    double hBA = parallelHDKP1(B, A, q, num_threads);
    return std::max(hAB, hBA);
}
