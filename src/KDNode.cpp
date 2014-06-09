#include "KDNode.h"

KDNode* KDNode::build(vector<Triangle*>& triangles, int depth) {
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
    
    double splitCoordinate = triangles[splitPos]->getMidpoint()[ax];

    vector<Triangle*> left;
    vector<Triangle*> right;
    
    int straddling = 0;
    
    for(auto t : triangles) {
        //If the triangle's bounding box starts before the split coordinate,
        //add it to the left subtree.
        bool inLeft = false;
        if (t->getBoundingBox().getPosition()[ax] <= splitCoordinate) {
            inLeft = true;
            left.push_back(t);
        }
        //If the triangle's BB ends after the split coordinate, add it to
        //the right subtree.
        if (t->getBoundingBox().getEndpoint()[ax] >= splitCoordinate) {
            right.push_back(t);
            if (inLeft) straddling++;
        }
        //This way, a triangle whose BB straddles the coordinate ends up in
        //both subtrees.
    }
    
    //If more than 50% of the triangles end up in both subtrees, make this node a leaf.
    if (straddling * 2 > triangles.size()) return node;
    
    //If the L/R nodes are a subset of R/L nodes, we get infinite recursion and so
    //turn the node into a leaf.
    if (straddling == right.size() || straddling == left.size()) return node;
    
    printf("L%d R%d S%d\n", left.size(), right.size(), straddling);
    
    node->left = KDNode::build(left, depth+1);
    node->right = KDNode::build(right, depth+1);
    
    return node;
}

//ignore hits that happened less than planeDist away from the ray origin
Intersection KDNode::getFirstIntersection(Ray r, double planeDist) {
    //If the ray doesn't intersect the bounding box, return
    if (!boundingBox.intersects(r)) {
        Intersection result;
        result.happened = false;
        return result;
    }
    
    //Base case: do an O(n) search through the triangles
    if (left == NULL && right == NULL) {
        Intersection bestInter;
        Intersection currInter;

        double mindist = 0;
        bool found = false;

        for(auto t : triangles) {
            currInter = t->getIntersection(r);
            
            if (currInter.happened) { 
                //Ignore hits that happened before the image plane
                if (currInter.distance < planeDist) continue;
                if (!found || (currInter.distance < mindist)) {
                    bestInter = currInter;
                    mindist = currInter.distance;
                    found = true;
                }
            }
        }

        if (!found) bestInter.happened = false;

        return bestInter;
    }
    
    //Recursive case: try the first child, then the second one
    //(TODO: try the closest one first instead)
    Intersection interL = left->getFirstIntersection(r, planeDist);
    Intersection interR = right->getFirstIntersection(r, planeDist);
    if (!interL.happened || (interR.happened && interL.distance > interR.distance))
        return interR;
    else return interL;
}