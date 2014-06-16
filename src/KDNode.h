#ifndef KDNODE_H
#define KDNODE_H
// A kd-tree structure to store triangles, similar to the one used to
// store photons but with explicit pointers.

#include "AABB.h"
#include "Triangle.h"
#include "Intersection.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include <set>
using namespace std;

class KDNode {
public:
    AABB boundingBox;
    KDNode* left;
    KDNode* right;
    vector<Triangle*> triangles;
    
    KDNode() {}
    
    static KDNode* build(vector<Triangle*>&triangles, int depth);
    
    Intersection getFirstIntersection(Ray r, double planeDist);
    
    set<Triangle*> getItems() {
        if (left == NULL) return set<Triangle*>(triangles.begin(), triangles.end());
        
        set<Triangle*> leftI = left->getItems();
        set<Triangle*> rightI = right->getItems();
        leftI.insert(rightI.begin(), rightI.end());
        
        return leftI;
    }
};

#endif