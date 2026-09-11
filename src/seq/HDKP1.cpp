#include "HDKP1.h"

namespace {

/**
 * Depth-first traversal of K_A visiting points in canonical Z-order without
 * storing large point lists in memory.
 */
void hdkp1DFS(MREP2 *A, long int sTree, SpatialArea area, MREP2 *B, double &cmax, Point &pNN) {
    // If leaf level in K_A, we reached an individual point p
    if (isLeaf(A, sTree)) {
        Point p;
        p.setX(area.x1);
        p.setY(area.y1);

        // Line 4: minNN = Dist(p, pNN)
        double minNN = distSinSqrt(p, pNN);

        // Line 5: if minNN > cmax then [Pruning Rule 1]
        if (minNN > cmax) {
            // Line 6: cmax = NNMAX(K_B, p, cmax, minNN)
            cmax = nnMax(B, p, cmax, minNN, pNN);
        }
        return;
    }

    // Internal node: explore non-empty children in Z-order (quadrants 0, 1, 2, 3)
    long int subTreeIndex = firstChild(A, sTree);
    for (int i = 0; i < (K * K); i++) {
        if (!esCero(A, subTreeIndex + i)) {
            hdkp1DFS(A, subTreeIndex + i, getSubArea(area, i), B, cmax, pNN);
        }
    }
}

} // anonymous namespace

/**
 * Algorithm 1: HDKP1(K_A, K_B)
 * Described in: "Efficient algorithms to calculate the Hausdorff distance on point sets represented by a k2-tree"
 * (Domínguez et al., GeoInformatica 2025)
 */
double hdkp1(MREP2 *A, MREP2 *B) {
    // Line 1: pNN = (inf, inf)
    Point pNN;
    pNN.setX(std::numeric_limits<uint32_t>::max() / 2);
    pNN.setY(std::numeric_limits<uint32_t>::max() / 2);

    // Line 2: cmax = 0
    double cmax = 0.0;

    SpatialArea rootArea;
    rootArea.x1 = 0;
    rootArea.x2 = pow(K, A->maxLevel + 1) - 1;
    rootArea.y1 = 0;
    rootArea.y2 = rootArea.x2;

    // Line 3-8: for Each point p in K_A do
    hdkp1DFS(A, -1, rootArea, B, cmax, pNN);

    // Line 9: return cmax (Euclidean distance)
    return std::sqrt(cmax);
}

double hdkp1(const std::vector<Point> &pointsA, MREP2 *B) {
    Point pNN;
    pNN.setX(std::numeric_limits<uint32_t>::max() / 2);
    pNN.setY(std::numeric_limits<uint32_t>::max() / 2);
    double cmax = 0.0;

    for (size_t i = 0; i < pointsA.size(); i++) {
        Point p = pointsA[i];
        double minNN = distSinSqrt(p, pNN);
        if (minNN > cmax) {
            cmax = nnMax(B, p, cmax, minNN, pNN);
        }
    }
    return std::sqrt(cmax);
}

double symmetricHDKP1(MREP2 *A, MREP2 *B) {
    double hAB = hdkp1(A, B);
    double hBA = hdkp1(B, A);
    return std::max(hAB, hBA);
}
