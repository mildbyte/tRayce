#include "KDNode.h"
#include <cassert>

double surfaceArea(AABB box) {
	Vector size = box.getSize();

	return 2 * (size[0] * size[1] + size[1] * size[2] + size[2] * size[0]);
}

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
    
	//Best currently: don't split
	int bestAx = -1;
	int bestPos = -1;
	double bestSAH = triangles.size() * surfaceArea(node->boundingBox);
	int bestStraddling = 0;

	for (int ax = 0; ax < 3; ax++) {

		//TODO: improve to a plane-sweeping algo
		/*sort(triangles.begin(), triangles.end(),
			[ax](Triangle* t1, Triangle* t2) {
			return t1->getMidpoint()[ax] < t2->getMidpoint()[ax];
			});*/

		//Iterate through all possible split positions
		for (int pos = 0; pos < triangles.size(); pos++) {
			double splitCoordinate = triangles[pos]->getMidpoint()[ax];

			int leftCount = 0;
			int rightCount = 0;

			AABB leftAABB;
			AABB rightAABB;

			int straddling = 0;

			//Count the hypothetical number of triangles in left and right halves
			for (auto t : triangles) {
				AABB box = t->getBoundingBox();

				bool inLeft = false;

				if (box.getStartpoint()[ax] <= splitCoordinate) {
					if (leftCount++ == 0) leftAABB = box; else leftAABB.addBox(box);
					inLeft = true;
				}
				if (box.getEndpoint()[ax] >= splitCoordinate) {
					if (rightCount++ == 0) rightAABB = box; else rightAABB.addBox(box);
					if (inLeft) straddling++;
				}
			}

			//Calculate the surface area heuristic
			double sah = surfaceArea(leftAABB) * leftCount + surfaceArea(rightAABB) * rightCount;

			//If we can improve and it's not a degenerate case (one side is a subset of the other one), record
			if (sah < bestSAH && straddling != leftCount && straddling != rightCount) {
				bestAx = ax;
				bestSAH = sah;
				bestPos = pos;

				bestStraddling = straddling;
			}

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

	//If no better split found, make ourselves a leaf.
	if (bestAx == -1) {
		node->triangles = triangles;
		return node;
	}

	//Recreate the split
	vector<Triangle*> left;
	vector<Triangle*> right;

	double splitCoordinate = triangles[bestPos]->getMidpoint()[bestAx];

	for (auto t : triangles) {
		AABB box = t->getBoundingBox();

		if (box.getStartpoint()[bestAx] <= splitCoordinate) left.push_back(t);
		if (box.getEndpoint()[bestAx] >= splitCoordinate) right.push_back(t);
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
	
    KDNode* first;
    KDNode* second;

    Intersection inter;
    inter.happened = false;
    
    if (lInter && (lDist - rDist < -100*EPSILON || !rInter)) {
        inter = left->getFirstIntersection(r, planeDist);
        if (!inter.happened && rInter) inter = right->getFirstIntersection(r, planeDist);
	}
	else if (rInter && (rDist - lDist < -100 * EPSILON || !lInter)) {
		inter = right->getFirstIntersection(r, planeDist);
		if (!inter.happened && lInter) inter = left->getFirstIntersection(r, planeDist);
	}
	else if (!rInter && !lInter) return inter;
	else {
		Intersection i1, i2;
		i1.happened = false;
		i2.happened = false;
		if (lInter) i1 = left->getFirstIntersection(r, planeDist);
		if (rInter) i2 = right->getFirstIntersection(r, planeDist);
		if (!i1.happened) return i2;
		if (!i2.happened) return i1;

		if (i1.distance < i2.distance) return i1; else return i2;
	}

    
    return inter;
}