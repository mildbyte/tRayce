#ifndef KDNODE_H
#define KDNODE_H
// A kd-tree structure to store triangles, similar to the one used to
// store photons but with explicit pointers.

#include "Box.h"
#include "Triangle.h"
#include <vector>
#include <algorithm>
using namespace std;

class KDNode {
public:
    Box boundingBox;
    KDNode* left;
    KDNode* right;
    vector<Triangle*> triangles;
    
    KDNode() {}
    
    KDNode* build(vector<Triangle*>&triangles, int depth) const;
};

#endif