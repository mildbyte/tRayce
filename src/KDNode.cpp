#include "KDNode.h"

KDNode* KDNode::build(vector<Triangle*>& triangles, int depth) const {
    KDNode* node = new KDNode();
    node->triangles = triangles;
    node->left = NULL;
    node->right = NULL;
    node->boundingBox = Box();
    
    //TODO: only need to store the triangles in the leaves
    
    //Calculate the bounding box
    for (auto it : triangles) {
        node->boundingBox.addBox(it->getBoundingBox());
    }
    
    //Leaf node: fewer than 16 triangles
    if (triangles.size() < 16) return node;
    
    //Choose a greatest spread axis and sort the triangles by it
    int ax = node->boundingBox.getGreatestSpread();
    sort(triangles.begin(), triangles.end(), 
        [ax](Triangle* t1, Triangle* t2) {
            return t1->getMidpoint()[ax] < t2->getMidpoint()[ax]; 
        });
    
    //Split so that the surface areas in the left and the right child
    //are similar (surface area heuristic)
    vector<double> cumulSA(triangles.size());
    cumulSA[0] = triangles[0]->getSurfaceArea();
    for (unsigned int i = 1; i < triangles.size(); i++) 
        cumulSA[i] = cumulSA[i-1] + triangles[i]->getSurfaceArea();
        
    double split = cumulSA[cumulSA.size()-1] / 2.0;
    int splitPos = distance(cumulSA.begin(), upper_bound(cumulSA.begin(), cumulSA.end(), split));
    
    double splitCoordinate = triangles[splitPos]->getBoundingBox().getEndpoint()[ax];

    vector<Triangle*> left;
    vector<Triangle*> right;
    
    for(auto t : triangles) {
        //If the triangle's bounding box ends before the split coordinate,
        //add it to the left subtree.
        if (t->getBoundingBox().getEndpoint()[ax] <= splitCoordinate)
            left.push_back(t);
        //If the triangle's BB begins after the split coordinate, add it to
        //the right subtree.
        if (t->getBoundingBox().getPosition()[ax] <= splitCoordinate)
            right.push_back(t);
        //This way, a triangle whose BB straddles the coordinage ends up in
        //both subtrees.
    }
    
    node->left = KDNode::build(left, depth+1);
    node->right = KDNode::build(right, depth+1);
    
    return node;
}