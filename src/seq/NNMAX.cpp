#include "NNMAX.h"

/**
 * Algorithm 2: NNMAX(K_B, p, cmax, minNN)
 * Described in: "Efficient algorithms to calculate the Hausdorff distance on point sets represented by a k2-tree"
 * (Domínguez et al., GeoInformatica 2025)
 */
double nnMax(MREP2 *rep, Point q, double cmax, double minNN, Point &pNN) {
    // Line 1: pQ = CreateMin-Heap()
    std::priority_queue<ElementQueue, std::vector<ElementQueue>, minHeapComparator> pQueue;

    // Line 2-3: d = MaxDist(p, K_B); Insert(pQ, <K_B, d>)
    ElementQueue e;
    e.sTree = -1; // Root node
    e.area.x1 = 0;
    e.area.x2 = pow(K, rep->maxLevel + 1) - 1;
    e.area.y1 = e.area.x1;
    e.area.y2 = e.area.x2;
    e.distancia = maxDistSinSqrt(q, e.area);
    pQueue.push(e);

    // Line 4: while (not Empty(pQ)) do
    while (!pQueue.empty()) {
        // Line 5: e = Extract-Min(pQ)
        e = pQueue.top();
        pQueue.pop();

        // Line 6: if not IsLeaf(e.n) then
        if (!isLeaf(rep, e.sTree)) {
            // Line 7-9: if (e.d <= cmax) then return cmax [Pruning Rule 2]
            if (e.distancia <= cmax) {
                return cmax;
            }

            // Line 10: for all Node h child of e.n do
            long int subTreeIndex = firstChild(rep, e.sTree);
            for (int i = 0; i < (K * K); i++) {
                // Line 11: if (hasChildren(h)) then
                if (!esCero(rep, subTreeIndex + i)) {
                    ElementQueue child;
                    child.sTree = subTreeIndex + i;
                    child.area = getSubArea(e.area, i);

                    // Eager Rule 2 check: if child's max distance is <= cmax,
                    // any point inside child cannot improve cmax.
                    child.distancia = maxDistSinSqrt(q, child.area);
                    if (child.distancia <= cmax) {
                        return cmax;
                    }

                    // Line 12-13: hminDist = MinDist(p, h); if (hminDist < minNN) then [Pruning Rule 3]
                    double hminDist = minDistSinSqrt(q, child.area);
                    if (hminDist < minNN) {
                        // Line 14-15: hmaxDist = MaxDist(p, h); Insert(pQ, <h, hmaxDist>)
                        pQueue.push(child);
                    }
                }
            }
        } else {
            // Line 19-22: Leaf level (individual point in set B)
            // if (e.d <= cmax) return cmax [Pruning Rule 1]
            if (e.distancia <= cmax) {
                return cmax;
            }

            // Line 23-26: if e.d < minNN then minNN = e.d; pNN = e.n
            if (e.distancia < minNN) {
                minNN = e.distancia;
                pNN.setX(e.area.x1);
                pNN.setY(e.area.y1);
            }
        }
    }

    // Line 29: return minNN
    return minNN;
}

// In-place overload matching the signature used in HDKP2
void nnMax(MREP2 *rep, Point q, double &max) {
    Point dummyPNN;
    double currentMax = max;
    double res = nnMax(rep, q, currentMax, std::numeric_limits<double>::max(), dummyPNN);
    if (res > max) {
        max = res;
    }
}
