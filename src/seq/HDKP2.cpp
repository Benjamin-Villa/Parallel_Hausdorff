#include "HDKP2.h"

namespace {

// Helper to find the first point in K_A to seed cmax (Algorithm 3, lines 1-2)
bool findFirstPoint(MREP2 *A, long int sTree, SpatialArea area, Point &p) {
    if (isLeaf(A, sTree)) {
        p.setX(area.x1);
        p.setY(area.y1);
        return true;
    }
    long int subTreeIndex = firstChild(A, sTree);
    for (int i = 0; i < (K * K); i++) {
        if (!esCero(A, subTreeIndex + i)) {
            if (findFirstPoint(A, subTreeIndex + i, getSubArea(area, i), p)) {
                return true;
            }
        }
    }
    return false;
}

} // anonymous namespace

/**
 * Algorithm 3: HDKP2(K_A, K_B)
 * Described in: "Efficient algorithms to calculate the Hausdorff distance on point sets represented by a k2-tree"
 * (Domínguez et al., GeoInformatica 2025)
 */
void hdkMaxHeapv2(MREP2 *A, MREP2 *B, double &supermax) {
    // Line 3: pQ = CreateMax-Heap()
    std::priority_queue<ElementQueue, std::vector<ElementQueue>, maxHeapComparatorHD> pQueue;

    // Root quadrant of K_A
    ElementQueue aa;
    aa.sTree = -1;
    aa.area.x1 = 0;
    uint32_t N = pow(K, A->maxLevel + 1);
    aa.area.x2 = N - 1;
    aa.area.y1 = 0;
    aa.area.y2 = N - 1;

    // Line 4-5: d = (N * sqrt(2))^2 = 2 * N^2 in squared Euclidean distance
    // Ensures consistent squared metric comparison with child distances from isCandidate
    aa.distancia = 2.0 * ((double)N * N);
    pQueue.push(aa);

    // Line 6: while (not Empty(pQ)) do
    while (!pQueue.empty()) {
        // Line 7: aa = Extract-Max(pQ)
        aa = pQueue.top();
        pQueue.pop();

        // Line 8-10: if aa.d <= cmax then return cmax
        // If the maximum possible distance in quadrant aa cannot exceed supermax,
        // no quadrant remaining in the max-heap can improve the solution.
        if (aa.distancia <= supermax) {
            return;
        }

        // Line 11: if IsLeaf(aa.n) then
        if (isLeaf(A, aa.sTree)) {
            Point q;
            q.setX(aa.area.x1);
            q.setY(aa.area.y1);

            // Line 12-15: nn = NNMAX(K_B, aa.n, cmax, inf); if nn > cmax then cmax = nn
            double aux = supermax;
            nnMax(B, q, supermax);
            if (aux >= supermax) {
                supermax = aux;
            }
        } else {
            // Line 16-24: Internal node: expand children and filter through IsCandidate
            long int subTreeIndex = firstChild(A, aa.sTree);
            for (int i = 0; i < (K * K); i++) {
                // Line 18: if (hasChildren(h)) then
                if (!esCero(A, subTreeIndex + i)) {
                    ElementQueue child;
                    child.sTree = subTreeIndex + i;
                    child.area = getSubArea(aa.area, i);

                    // Line 19: dh = IsCandidate(h, K_B, cmax)
                    // Line 20-21: if (dh != -1) then Insert(pQ, <h, dh>)
                    double dh = isCandidate(child.area, B, supermax);
                    if (dh != -1.0) {
                        child.distancia = dh;
                        pQueue.push(child);
                    }
                }
            }
        }
    }
}

double hdkp2(MREP2 *A, MREP2 *B) {
    double supermax = 0.0;

    // Line 1-2: Seed cmax with nearest neighbor of any point p in K_A
    SpatialArea rootArea;
    rootArea.x1 = 0;
    rootArea.x2 = pow(K, A->maxLevel + 1) - 1;
    rootArea.y1 = 0;
    rootArea.y2 = rootArea.x2;

    Point initialPoint;
    if (findFirstPoint(A, -1, rootArea, initialPoint)) {
        nnMax(B, initialPoint, supermax);
    }

    // Run HDKP2 max-heap search
    hdkMaxHeapv2(A, B, supermax);

    // Line 27: return cmax (Euclidean distance)
    return std::sqrt(supermax);
}

double symmetricHDKP2(MREP2 *A, MREP2 *B) {
    double hAB = hdkp2(A, B);
    double hBA = hdkp2(B, A);
    return std::max(hAB, hBA);
}
