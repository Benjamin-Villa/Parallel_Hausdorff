#ifndef HDKP2_H
#define HDKP2_H

#include "common_hdk.h"
#include "NNMAX.h"
#include "IsCandidate.h"

/**
 * Algorithm 3: HDKP2(K_A, K_B)
 *
 * Computes directed Hausdorff distance h(A, B) using hierarchical branch-and-prune
 * over both point sets A and B stored in k2-trees.
 *
 * Uses a max-heap pQ to explore the most promising quadrants of K_A first.
 * Discards entire quadrants of K_A using IsCandidate (Pruning Rule 4).
 * For leaf points, queries K_B using NNMAX (Pruning Rules 1, 2, 3).
 * Terminates early when the maximum possible distance of the top quadrant is <= cmax.
 *
 * @param A   k2-tree representing point set A (K_A)
 * @param B   k2-tree representing point set B (K_B)
 * @return    Hausdorff distance h(A, B) (Euclidean distance)
 */
double hdkp2(MREP2 *A, MREP2 *B);

// Core internal worker tracking maximum heap statistics
void hdkMaxHeapv2(MREP2 *A, MREP2 *B, double &supermax);

// Symmetrical Hausdorff distance: max(h(A, B), h(B, A))
double symmetricHDKP2(MREP2 *A, MREP2 *B);

#endif // HDKP2_H
