#include "KDNode.h"

KDNode* KDNode::build(vector<Triangle*>& triangles, int depth) {
    KDNode* node = new KDNode();
    node->triangles = triangles;
    node->left = NULL;
    node->right = NULL;
    
    //Use the first triangle's bounding box as base since 0,0,0->0,0,0 still
    //keeps the end at 0,0,0 even if all triangles have negative coordinates.
    //Are we guaranteed to have > 0 triangles?
    node->boundingBox = triangles[0]->getBoundingBox();
    //TODO: only need to store the triangles in the leaves
    
    printf("level %d, %d triangles\n", depth, triangles.size());
    
    //Calculate the bounding box
    for (auto it = triangles.begin() + 1; it != triangles.end(); it++) {
        node->boundingBox.addBox((*it)->getBoundingBox());
    }
    
    //Leaf node: fewer than 16 triangles
    if (triangles.size() < 16) return node;
    
    //Choose a greatest spread axis and sort the triangles by it
    int ax = node->boundingBox.getGreatestSpread();
    sort(triangles.begin(), triangles.end(), 
        [ax](Triangle* t1, Triangle* t2) {
            return t1->getMidpoint()[ax] < t2->getMidpoint()[ax]; 
        });
        
    printf("splitting on %d\n", ax);
    printf("bounding box: from ");
    node->boundingBox.getPosition().print();
    printf(" to ");
    node->boundingBox.getEndpoint().print();
    printf("\n");
    
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
    
    // Try the closest-intersecting box first, if nothing, move on to the second one
    double lDist;
    double rDist;
    
    /*
    printf("Ray: ");
    r.origin.print();
    printf(" -> ");
    r.direction.print();
    printf("\n");
    
    printf("Left box: ");
    left->boundingBox.getPosition().print();
    printf(" -> ");
    left->boundingBox.getEndpoint().print();
    
    printf("\nRight box: ");
    right->boundingBox.getPosition().print();
    printf(" -> ");
    right->boundingBox.getEndpoint().print();
    */
    
    bool lInter = left->boundingBox.intersects(r, lDist);
    bool rInter = right->boundingBox.intersects(r, rDist);
    
    KDNode* first;
    KDNode* second;
    
    if (lInter && (lDist <= rDist || !rInter)) {
    //    printf("\nLeft intersected first\n");
        first = left;
        if (!rInter) second = NULL;
        else second = right;
    } else if (rInter && (rDist < lDist || !lInter)) {
    //    printf("\nRight intersected first\n");
        first = right;
        if (!lInter) second = NULL;
        else second = left;
    } else {
    //    printf("\nNo intersections\n");
        Intersection inter;
        inter.happened = false;
        return inter;
    }
    
    Intersection inter = first->getFirstIntersection(r, planeDist);
    if (!inter.happened && second != NULL) inter = second->getFirstIntersection(r, planeDist);
    
    return inter;
}