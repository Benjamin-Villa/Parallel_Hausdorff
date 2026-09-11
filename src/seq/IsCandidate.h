#ifndef ISCANDIDATE_H
#define ISCANDIDATE_H

#include "common_hdk.h"

/**
 * Algorithm 4: IsCandidate(nodeA, K_B, cmax)
 *
 * Determines if quadrant/region nodeA in set A contains candidate points that can
 * improve the Hausdorff distance candidate cmax.
 *
 * Uses a min-heap pR over nodes of K_B prioritized by MAXMAXDIST(nodeA, subregion).
 * Applies Pruning Rule 4:
 *  - If MAXMAXDIST(nodeA, child_B) <= cmax, no point in nodeA can have a nearest
 *    neighbor distance > cmax in set B, so the entire region nodeA is discarded (-1).
 *
 * @param nodeArea  Bounding box of region in set A (nodeA)
 * @param B         k2-tree representing point set B (K_B)
 * @param cmax      Current Hausdorff candidate distance (squared Euclidean)
 * @return          MAXMAXDIST(nodeA, S) as refined upper bound, or -1 if pruned
 */
double isCandidate(SpatialArea nodeArea, MREP2 *B, double cmax);

// Overload taking ElementQueue for direct integration with HDKP2 queue
double isCandidate(ElementQueue eq, MREP2 *B, double &supermax);

#endif // ISCANDIDATE_H
