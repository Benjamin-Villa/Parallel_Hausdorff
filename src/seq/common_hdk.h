#ifndef COMMON_HDK_H
#define COMMON_HDK_H

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <queue>
#include <vector>
#include <limits>
#include <Point.h>
#include <K2tree.h>

#ifndef K
#define K 2
#endif

// 2D bounding box representing a quadrant / submatrix in the k2-tree
typedef struct sArea {
    uint32_t x1, x2, y1, y2;
} SpatialArea;

// Element for priority queues (min-heap or max-heap)
typedef struct sElementQueue {
    long int sTree;      // Subtree index / position in bitmap (-1 for root)
    SpatialArea area;    // Spatial bounding region
    double distancia;    // Distance priority (squared Euclidean distance)
} ElementQueue;

// Min-heap comparator (smallest distance at top)
struct minHeapComparator {
    bool operator()(const ElementQueue &i, const ElementQueue &j) const {
        return i.distancia > j.distancia;
    }
};

// Max-heap comparator (largest distance at top)
struct maxHeapComparatorHD {
    bool operator()(const ElementQueue &i, const ElementQueue &j) const {
        return i.distancia < j.distancia;
    }
};

// Subdivides an area into one of K x K quadrants (0: top-left, 1: bottom-left, 2: top-right, 3: bottom-right)
inline SpatialArea getSubArea(SpatialArea a, int i) {
    SpatialArea resp;
    uint32_t midX = (a.x1 + a.x2) / 2;
    uint32_t midY = (a.y1 + a.y2) / 2;
    switch (i) {
        case 0:
            resp.x1 = a.x1;
            resp.x2 = midX;
            resp.y1 = a.y1;
            resp.y2 = midY;
            break;
        case 1:
            resp.x1 = a.x1;
            resp.x2 = midX;
            resp.y1 = midY + 1;
            resp.y2 = a.y2;
            break;
        case 2:
            resp.x1 = midX + 1;
            resp.x2 = a.x2;
            resp.y1 = a.y1;
            resp.y2 = midY;
            break;
        case 3:
            resp.x1 = midX + 1;
            resp.x2 = a.x2;
            resp.y1 = midY + 1;
            resp.y2 = a.y2;
            break;
        default:
            resp = a;
            break;
    }
    return resp;
}

// maxDist(p, S): Maximum squared Euclidean distance between point p and any point in subregion a
inline double maxDistSinSqrt(Point p, SpatialArea a) {
    long long diff_x1 = (long long)p.getX() - a.x1;
    long long diff_x2 = (long long)p.getX() - a.x2;
    uint64_t dx = std::max((uint64_t)std::abs(diff_x1), (uint64_t)std::abs(diff_x2));

    long long diff_y1 = (long long)p.getY() - a.y1;
    long long diff_y2 = (long long)p.getY() - a.y2;
    uint64_t dy = std::max((uint64_t)std::abs(diff_y1), (uint64_t)std::abs(diff_y2));

    return ((double)dx * dx) + ((double)dy * dy);
}

// minDist(p, S): Minimum squared Euclidean distance between point p and any point in subregion a
inline double minDistSinSqrt(Point p, SpatialArea a) {
    uint64_t x = p.getX();
    uint64_t y = p.getY();

    if (x >= a.x1 && x <= a.x2 && y >= a.y1 && y <= a.y2) {
        return 0.0;
    }

    uint64_t dx = 0;
    if (x < a.x1) dx = a.x1 - x;
    else if (x > a.x2) dx = x - a.x2;

    uint64_t dy = 0;
    if (y < a.y1) dy = a.y1 - y;
    else if (y > a.y2) dy = y - a.y2;

    return ((double)dx * dx) + ((double)dy * dy);
}

// MAXMAXDIST(R, S): Maximum squared Euclidean distance between any point in region s and any point in region b
inline double maxDistAreas(SpatialArea s, SpatialArea b) {
    Point p;
    double d = 0.0, dTest;

    // Corner 1: (s.x1, s.y1)
    p.setX(s.x1); p.setY(s.y1);
    d = maxDistSinSqrt(p, b);

    // Corner 2: (s.x1, s.y2)
    p.setX(s.x1); p.setY(s.y2);
    dTest = maxDistSinSqrt(p, b);
    if (dTest > d) d = dTest;

    // Corner 3: (s.x2, s.y1)
    p.setX(s.x2); p.setY(s.y1);
    dTest = maxDistSinSqrt(p, b);
    if (dTest > d) d = dTest;

    // Corner 4: (s.x2, s.y2)
    p.setX(s.x2); p.setY(s.y2);
    dTest = maxDistSinSqrt(p, b);
    if (dTest > d) d = dTest;

    return d;
}

// Squared Euclidean distance between two points
inline double distSinSqrt(Point p1, Point p2) {
    long long dx = (long long)p1.getX() - (long long)p2.getX();
    long long dy = (long long)p1.getY() - (long long)p2.getY();
    return ((double)dx * dx) + ((double)dy * dy);
}

#endif // COMMON_HDK_H
