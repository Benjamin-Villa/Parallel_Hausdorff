#ifndef PARALLEL_HDKP1_H
#define PARALLEL_HDKP1_H

#include <sycl/sycl.hpp>
#include <vector>
#include "Util/Point.h"
#include "K2tree.h"

/**
 * Initializes and warms up the SYCL runtime outside of benchmark timers.
 */
void warmupSYCL(sycl::queue &q);

/**
 * Naive Embarrassing Parallel HDKP1:
 * Computes directed Hausdorff distance h(A, B) = max_{p in A} min_{q in B} dist(p, q)
 * by extracting points from K_A and parallelizing the outer point loop across workers.
 *
 * Each worker maintains a private pNN cache for Pruning Rule 1 and calls the original,
 * untouched sequential nnMax() from seq/NNMAX.cpp as a black box.
 * Shared lower bound cmax is monotonically updated via sycl::atomic_ref::fetch_max.
 *
 * @param A            k2-tree representing point set A (K_A)
 * @param B            k2-tree representing point set B (K_B)
 * @param q            Optional pointer to existing SYCL queue (nullptr creates a new CPU queue)
 * @param num_threads  Number of parallel workers (0 = auto-detect hardware concurrency)
 * @return             Directed Hausdorff distance h(A, B) (Euclidean distance)
 */
double parallelHDKP1(MREP2 *A, MREP2 *B, sycl::queue *q = nullptr, unsigned int num_threads = 0);

// Overload taking an explicit list of points for set A
double parallelHDKP1(const std::vector<Point> &pointsA, MREP2 *B, sycl::queue *q = nullptr, unsigned int num_threads = 0);

// Symmetrical Hausdorff distance: max(h(A, B), h(B, A))
double symmetricParallelHDKP1(MREP2 *A, MREP2 *B, sycl::queue *q = nullptr, unsigned int num_threads = 0);

#endif // PARALLEL_HDKP1_H
