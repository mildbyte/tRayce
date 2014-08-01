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

bool planeIntersectsRay(SplitPlane plane, Ray r, double &dist) {
	double rayO = r.origin[plane.getAxis()];
	double rayC = r.direction[plane.getAxis()];
	if (abs(rayC) < EPSILON) return false;

	dist = (plane.getCoordinate() - rayO) / rayC;
	return true;
}

//Performs a plane sweep across the triangles to find the best splitting plane using the SAH.
//From "On building fast kd-Trees for Ray Tracing, and on doing that in O(N log N)" by I. Wald and V. Havran
pair<pair<SplitPlane, SplitSide>, double> KDNode::findPlane(vector<Triangle*>& triangles, AABB boundingBox) {
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
	AABB boundingBox;

	for (auto it : triangles) {
		if (firstTriangle) {
			firstTriangle = false;
			//Use the first triangle's bounding box as base since 0,0,0->0,0,0 still
			//keeps the end at 0,0,0 even if all triangles have negative coordinates.
			//Are we guaranteed to have > 0 triangles?
			boundingBox = it->getBoundingBox();
		}
		boundingBox.addBox(it->getBoundingBox());
	}

	//Leaf node: fewer than 4 triangles or reached sufficient depth
	if (triangles.size() < 4 || depth >= depthLimit) {
		node->triangles = triangles;
		return node;
	}


	pair<pair<SplitPlane, SplitSide>, double> split = node->findPlane(triangles, boundingBox);

	//If can't get a better SAH, by splitting, make a leaf
	if (split.second >= surfaceArea(boundingBox) * triangles.size()) {
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
	node->plane = split.first.first;

	return node;

}

KDNode* KDNode::build(vector<Triangle*>& triangles) {
	int limit = 1; int tmp = 1;
	while (tmp < triangles.size()) {
		tmp *= 2; limit++;
	}
	return limitedBuild(triangles, 0, limit + 2);
}

//ignore hits that happened less than planeDist away from the ray origin
Intersection KDNode::getFirstIntersection(Ray r, double tMin, double tMax) {
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
                if (currInter.distance < tMin) continue;
                if (!found || (currInter.distance < mindist)) {
                    bestInter = currInter;
                    mindist = currInter.distance;
                    found = true;
                }
            }
        }

#ifdef _DEBUG
		bestInter.intersectedTriangles = triangles.size();
#endif

        if (!found) bestInter.happened = false;

        return bestInter;
    }

	double tSplit;
	planeIntersectsRay(plane, r, tSplit);

	KDNode* near;
	KDNode* far;
	if (r.origin[plane.getAxis()] < plane.getCoordinate()) near = left, far = right; else near = right, far = left;


	if (tSplit > tMax) {
		return near->getFirstIntersection(r, tMin, tMax);
	}
	else if (tSplit < tMin) {
		if (tSplit > 0) return far->getFirstIntersection(r, tMin, tMax);
		else if (tSplit < 0) return near->getFirstIntersection(r, tMin, tMax);
		else {
			if (r.direction[plane.getAxis()] < 0) return far->getFirstIntersection(r, tMin, tMax);
			else return near->getFirstIntersection(r, tMin, tMax);
		}
	}
	else {
		if (tSplit > 0) {
			Intersection inter = near->getFirstIntersection(r, tMin, tSplit);
			if (inter.happened) {
				if (inter.distance <= tSplit) return inter;
				Intersection inter2 = far->getFirstIntersection(r, tSplit, tMax);

#ifdef _DEBUG
				inter.intersectedTriangles += inter2.intersectedTriangles;
				inter2.intersectedTriangles = inter.intersectedTriangles;
#endif

				if (inter2.happened && inter2.distance < inter.distance) return inter2; else return inter;
			}
			return far->getFirstIntersection(r, tSplit, tMax);
		}
		else return near->getFirstIntersection(r, tSplit, tMax);
	}
}