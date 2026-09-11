#ifndef NNMAX_H
#define NNMAX_H

#include "common_hdk.h"

/**
 * Algorithm 2: NNMAX (Maximization of the Hausdorff distance)
 *
 * Computes the nearest neighbor of point p in k2-tree K_B, optimizing the Hausdorff distance candidate cmax.
 * Uses a min-heap prioritized by maxDist(p, quadrant).
 *
 * Applies Pruning Rules:
 *  - Rule 2: If an internal node has maxDist(p, h) <= cmax, returns cmax (early break).
 *  - Rule 3: If child node has minDist(p, h) > minNN, child is discarded (not inserted into heap).
 *  - Rule 1: If leaf node (point) has dist(p, pb) <= cmax, returns cmax.
 *
 * @param rep   k2-tree representing point set B (K_B)
 * @param q     Query point p from set A
 * @param cmax  Current Hausdorff candidate distance (squared Euclidean)
 * @param minNN Upper bound for distance between p and its nearest neighbor (squared Euclidean)
 * @param pNN   (Output) Coordinates of nearest neighbor in B found so far
 * @return      Updated candidate distance minNN or cmax
 */
double nnMax(MREP2 *rep, Point q, double cmax, double minNN, Point &pNN);

// In-place overload where max is both candidate cmax and updated nearest neighbor distance
void nnMax(MREP2 *rep, Point q, double &max);

#endif // NNMAX_H
