#include "KDNode.h"
#include <cassert>

KDNode* KDNode::build(vector<Triangle*>& triangles, int depth) {
    KDNode* node = new KDNode();
    node->left = NULL;
    node->right = NULL;

#ifdef _DEBUG
    printf("level %d, %d triangles\n", depth, triangles.size());
#endif
    
    //Calculate the bounding box
    bool firstTriangle = true;
    
    for (auto it : triangles) {
        if (firstTriangle) {
            firstTriangle = false;
            //Use the first triangle's bounding box as base since 0,0,0->0,0,0 still
            //keeps the end at 0,0,0 even if all triangles have negative coordinates.
            //Are we guaranteed to have > 0 triangles?
            node->boundingBox = it->getBoundingBox();
        }
        node->boundingBox.addBox(it->getBoundingBox());
    }
    
    //Leaf node: fewer than 16 triangles
	if (triangles.size() < 4) {
		node->triangles = triangles;
		return node;
	}
    
	
	int bestAx = -1;
	int bestStraddle = 0;

	//Find the axis with which as few triangles straddle the boundary as possible
	for (int ax = 0; ax < 3; ax++) {
		sort(triangles.begin(), triangles.end(),
			[ax](Triangle* t1, Triangle* t2) {
			return t1->getMidpoint()[ax] < t2->getMidpoint()[ax];
		});

		//Split so that the surface areas in the left and the right child
		//are similar (surface area heuristic)
		vector<double> cumulSA(triangles.size());
		cumulSA[0] = triangles[0]->getSurfaceArea();
		for (unsigned int i = 1; i < triangles.size(); i++)
			cumulSA[i] = cumulSA[i - 1] + triangles[i]->getSurfaceArea();

		double split = cumulSA[cumulSA.size() - 1] / 2.0;
		int splitPos = distance(cumulSA.begin(), upper_bound(cumulSA.begin(), cumulSA.end(), split));

		double splitCoordinate = triangles[splitPos]->getBoundingBox().getStartpoint()[ax];
		int straddling = 0;

		for (auto t : triangles) {
			if ((t->getBoundingBox().getStartpoint()[ax] <= splitCoordinate)
				&& (t->getBoundingBox().getEndpoint()[ax] >= splitCoordinate)) straddling++;
		}

		if (bestAx == -1 || straddling < bestStraddle) {
			bestAx = ax;
			bestStraddle = straddling;
		}	
	}

#ifdef _DEBUG
	printf("splitting on %d\n", bestAx);
	printf("bounding box: from ");
	node->boundingBox.getStartpoint().print();
	printf(" to ");
	node->boundingBox.getEndpoint().print();
	printf("\n");
#endif
	

	//Now do it for real
	//TODO: reuse the results of the best split

	sort(triangles.begin(), triangles.end(),
		[bestAx](Triangle* t1, Triangle* t2) {
		return t1->getMidpoint()[bestAx] < t2->getMidpoint()[bestAx];
	});

	//Split so that the surface areas in the left and the right child
	//are similar (surface area heuristic)
	vector<double> cumulSA(triangles.size());
	cumulSA[0] = triangles[0]->getSurfaceArea();
	for (unsigned int i = 1; i < triangles.size(); i++)
		cumulSA[i] = cumulSA[i - 1] + triangles[i]->getSurfaceArea();

	double split = cumulSA[cumulSA.size() - 1] / 2.0;
	int splitPos = distance(cumulSA.begin(), upper_bound(cumulSA.begin(), cumulSA.end(), split));
    
    double splitCoordinate = triangles[splitPos]->getBoundingBox().getStartpoint()[bestAx];

    vector<Triangle*> left;
    vector<Triangle*> right;
    
    for(auto t : triangles) {
        //If the triangle's bounding box starts before the split coordinate,
        //add it to the left subtree.
        bool inLeft = false;
        
        if (t->getBoundingBox().getStartpoint()[bestAx] <= splitCoordinate) {
            left.push_back(t);
        }
        //If the triangle's BB ends after the split coordinate, add it to
        //the right subtree.
        if (t->getBoundingBox().getEndpoint()[bestAx] >= splitCoordinate) {
            right.push_back(t);
        }
        //This way, a triangle whose BB straddles the coordinate ends up in
        //both subtrees.
    }
    
    //If more than 50% of the triangles end up in both subtrees, make this node a leaf.
	//If the L/R nodes are a subset of R/L nodes, we get infinite recursion and so
	//turn the node into a leaf.
	if (bestStraddle * 2 > triangles.size() || bestStraddle == right.size() || bestStraddle == left.size()) {
		node->triangles = triangles;
		return node;
	}
    
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
#ifdef _DEBUG
                AABB bb = t->getBoundingBox();
                double d;
                if (!bb.intersects(r, d)) {
                    printf("ASSERTION FAILED: triangle ");
                    t->print();
                    printf("; ray ");
                    r.origin.print();
                    printf(" -> ");
                    r.direction.print();
                    printf("; box ");
                    bb.getStartpoint().print();
                    printf(" -> ");
                    bb.getEndpoint().print();
                    printf(": Triangle intersected, BB wasn't\n");
                }
#endif
            }
        }

#ifdef _DEBUG
		bestInter.intersectedTriangles = triangles.size();
#endif

        if (!found) bestInter.happened = false;

        return bestInter;
    }
    
//#ifdef _DEBUG
//	printf("Ray: ");
//    r.origin.print();
//    printf(" -> ");
//    r.direction.print();
//    printf("\n");
//    
//    printf("Left box: ");
//    left->boundingBox.getStartpoint().print();
//    printf(" -> ");
//    left->boundingBox.getEndpoint().print();
//    
//    printf("\nRight box: ");
//    right->boundingBox.getStartpoint().print();
//    printf(" -> ");
//    right->boundingBox.getEndpoint().print();
//#endif
    
    // Try the closest-intersecting box first, if nothing, move on to the second one
	// TODO: doesn't completely work, sometimes returns later intersections if the ray goes
	// through where two boxes merge.

    double lDist;
    double rDist;
    
    bool lInter = left->boundingBox.intersects(r, lDist);
    bool rInter = right->boundingBox.intersects(r, rDist);
	/*
    KDNode* first;
    KDNode* second;

    Intersection inter;
    inter.happened = false;
    
    if (lInter && (lDist <= rDist || !rInter)) {
        inter = left->getFirstIntersection(r, planeDist);
        if (!inter.happened && rInter) inter = right->getFirstIntersection(r, planeDist);
    } else if (rInter && (rDist <= lDist || !lInter)) {
        inter = right->getFirstIntersection(r, planeDist);
        if (!inter.happened && lInter) inter = left->getFirstIntersection(r, planeDist);
    }
    
    return inter;
	*/
    
    Intersection i1, i2;
    i1.happened = false;
    i2.happened = false;

#ifdef _DEBUG
	i1.intersectedTriangles = 0;
	i2.intersectedTriangles = 0;
#endif
    
    if (lInter) i1 = left->getFirstIntersection(r, planeDist);
    if (rInter) i2 = right->getFirstIntersection(r, planeDist);

#ifdef _DEBUG
	i1.intersectedTriangles += i2.intersectedTriangles;
	i2.intersectedTriangles = i1.intersectedTriangles;
#endif
    
	if (!i1.happened) return i2;
	if (!i2.happened) return i1;
	
	if (i1.distance < i2.distance) return i1; else return i2;
}