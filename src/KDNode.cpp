#include "KDNode.h"
#include <cassert>

double surfaceArea(AABB box) {
	Vector size = box.getSize();

	return 2 * (size[0] * size[1] + size[1] * size[2] + size[2] * size[0]);
}

double cost(double ratioLeft, double ratioRight, int left, int right) {
	return ratioLeft * left + ratioRight * right;
}

pair<double, SplitSide> SAH(SplitPlane plane, AABB box, int left, int right, int planar) {
	pair<AABB, AABB> boxes = box.split(plane);
	double area = surfaceArea(box);
	double ratioLeft = surfaceArea(boxes.first) / area;
	double ratioRight = surfaceArea(boxes.second) / area;

	double costLeft = cost(ratioLeft, ratioRight, left + planar, right);
	double costRight = cost(ratioLeft, ratioRight, left, planar + right);

	if (costLeft < costRight) return make_pair(costLeft, LEFT);
	else return make_pair(costRight, RIGHT);
}

//Performs a plane sweep across the triangles to find the best splitting plane using the SAH.
//From "On building fast kd-Trees for Ray Tracing, and on doing that in O(N log N)" by I. Wald and V. Havran
pair<pair<SplitPlane, SplitSide>, double> KDNode::findPlane(vector<Triangle*>& triangles) {
	double bestCost = INFINITY;
	SplitPlane bestPlane(0, 0);
	SplitSide bestSide;

	for (int dim = 0; dim < 3; dim++) {
		vector<SweepEvent> events;

		for (auto t : triangles) {
			AABB b = t->getBoundingBox().intersection(boundingBox);

			if (b.getSize()[dim] == 0.0) {
				events.push_back(SweepEvent(t, b.getStartpoint()[dim], PLANAR));
			}
			else {
				events.push_back(SweepEvent(t, b.getStartpoint()[dim], START));
				events.push_back(SweepEvent(t, b.getEndpoint()[dim], END));
			}
		}
		sort(events.begin(), events.end());

		int left = 0;
		int planar = 0;
		int right = triangles.size();

		for (int i = 0; i < events.size(); i++) {
			double coordinate = events[i].coordinate;
			SplitPlane plane(coordinate, dim);

			int pAdd = 0; int pPlan = 0; int pRem = 0;

			while (i < events.size() && events[i].coordinate == coordinate && events[i].type == END) {
				pRem++; i++;
			}
			while (i < events.size() && events[i].coordinate == coordinate && events[i].type == PLANAR) {
				pPlan++; i++;
			}
			while (i < events.size() && events[i].coordinate == coordinate && events[i].type == START) {
				pAdd++; i++;
			}

			planar = pPlan; right -= pPlan; right -= pRem;

			pair<double, SplitSide> result = SAH(plane, boundingBox, left, right, planar);
			
			if (result.first < bestCost) {
				bestCost = result.first;
				bestPlane = plane;
				bestSide = result.second;
			}
			left += pAdd; left += pPlan; planar = 0;
		}
	}
	return make_pair(make_pair(bestPlane, bestSide), bestCost);
}

KDNode* KDNode::limitedBuild(vector<Triangle*>& triangles, int depth, int depthLimit) {
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

	//Leaf node: fewer than 4 triangles or reached sufficient depth
	if (triangles.size() < 4 || depth >= depthLimit) {
		node->triangles = triangles;
		return node;
	}


	pair<pair<SplitPlane, SplitSide>, double> split = node->findPlane(triangles);

	//If can't get a better SAH, by splitting, make a leaf
	if (split.second >= surfaceArea(node->boundingBox) * triangles.size()) {
		node->triangles = triangles;
		return node;
	}

	double splitCoordinate = split.first.first.getCoordinate();
	vector<Triangle*> left;
	vector<Triangle*> right;

	for (auto t : triangles) {
		double triangleStart = t->getBoundingBox().getStartpoint()[split.first.first.getAxis()];
		double triangleEnd = t->getBoundingBox().getEndpoint()[split.first.first.getAxis()];

		if (triangleStart == triangleEnd && triangleStart == splitCoordinate)
			if (split.first.second == LEFT) left.push_back(t); else right.push_back(t);
		else {
			if (triangleStart < splitCoordinate) left.push_back(t);
			if (triangleEnd > splitCoordinate) right.push_back(t);
		}
	}

	//TODO deal with one side being a subset of another.
	node->left = limitedBuild(left, depth + 1, depthLimit);
	node->right = limitedBuild(right, depth + 1, depthLimit);


#ifdef _DEBUG
	printf("bounding box: from ");
	node->boundingBox.getStartpoint().print();
	printf(" to ");
	node->boundingBox.getEndpoint().print();
	printf("\n");
#endif

	return node;

}

KDNode* KDNode::build(vector<Triangle*>& triangles) {
	int limit = 1; int tmp = 1;
	while (tmp < triangles.size()) {
		tmp *= 2; limit++;
	}
	return limitedBuild(triangles, 0, limit);
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