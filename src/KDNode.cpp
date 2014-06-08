#include "KDNode.h"

KDNode* KDNode::build(vector<Triangle*>& triangles, int depth) const {
    KDNode* node = new KDNode();
    node->triangles = triangles;
    node->left = NULL;
    node->right = NULL;
    node->boundingBox = new Box();
    
    //Calculate the bounding box
    
    //Base case: return
    
    //Recursive case: split into 2 triangle vectors (using SAH), build subtrees, set the children
}