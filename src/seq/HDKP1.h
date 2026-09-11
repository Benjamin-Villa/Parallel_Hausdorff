#ifndef HDKP1_H
#define HDKP1_H

#include "common_hdk.h"
#include "NNMAX.h"

/**
 * Algorithm 1: HDKP1(K_A, K_B)
 *
 * Computes directed Hausdorff distance h(A, B) = max_{p \in A} min_{q \in B} dist(p, q)
 * where points are stored in k2-trees K_A and K_B.
 *
 * Traverses points in K_A in depth-first (Z-order), maintains the global nearest neighbor
 * cache pNN, and invokes NNMAX only when Dist(p, pNN) > cmax (Pruning Rule 1).
 *
 * @param A   k2-tree representing point set A (K_A)
 * @param B   k2-tree representing point set B (K_B)
 * @return    Hausdorff distance h(A, B) (Euclidean distance)
 */
double hdkp1(MREP2 *A, MREP2 *B);

// Overload taking an explicit list of points for set A
double hdkp1(const std::vector<Point> &pointsA, MREP2 *B);

// Symmetrical Hausdorff distance: max(h(A, B), h(B, A))
double symmetricHDKP1(MREP2 *A, MREP2 *B);

#endif // HDKP1_H
