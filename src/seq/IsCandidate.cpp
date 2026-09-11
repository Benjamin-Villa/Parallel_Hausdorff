#include "IsCandidate.h"

/**
 * Algorithm 4: IsCandidate(nodeA, K_B, cmax)
 * Described in: "Efficient algorithms to calculate the Hausdorff distance on point sets represented by a k2-tree"
 * (Domínguez et al., GeoInformatica 2025)
 */
double isCandidate(SpatialArea nodeArea, MREP2 *B, double cmax) {
    // Line 1: pR = CreateMin-Heap()
    std::priority_queue<ElementQueue, std::vector<ElementQueue>, minHeapComparator> pR;

    // Line 2-3: d = MAXMAXDIST(nodeA, K_B); Insert(pR, <K_B, d>)
    ElementQueue e;
    e.sTree = -1; // Root of B
    e.area.x1 = 0;
    e.area.x2 = pow(K, B->maxLevel + 1) - 1;
    e.area.y1 = e.area.x1;
    e.area.y2 = e.area.x2;
    // Properly initialized MAXMAXDIST for the root bounding box (fixes uninitialized bug from legacy codebase)
    e.distancia = maxDistAreas(nodeArea, e.area);
    pR.push(e);

    // Line 4: while (not Empty(pR)) do
    while (!pR.empty()) {
        // Line 5: bb = Extract-Min(pR)
        ElementQueue bb = pR.top();
        pR.pop();

        // Line 6-7: if IsLeaf(bb.n) then return bb.d
        if (isLeaf(B, bb.sTree)) {
            return bb.distancia;
        }

        // Line 8-18: Non-leaf internal node
        long int subTreeIndex = firstChild(B, bb.sTree);
        for (int i = 0; i < (K * K); i++) {
            // Line 10: if (hasChildren(h)) then
            if (!esCero(B, subTreeIndex + i)) {
                ElementQueue h;
                h.sTree = subTreeIndex + i;
                h.area = getSubArea(bb.area, i);

                // Line 11: maxDist = MAXMAXDIST(nodeA, h)
                h.distancia = maxDistAreas(nodeArea, h.area);

                // Line 12-14: if (maxDist <= cmax) then return -1 [Pruning Rule 4]
                if (h.distancia <= cmax) {
                    return -1.0;
                }

                // Line 15: Insert(pR, <h, maxDist>)
                pR.push(h);
            }
        }
    }

    // Line 19: Heap exhausted without reaching a leaf point
    return -1.0;
}

// Overload taking ElementQueue matching HDKP2 usage
double isCandidate(ElementQueue eq, MREP2 *B, double &supermax) {
    return isCandidate(eq.area, B, supermax);
}
