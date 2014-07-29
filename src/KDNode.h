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

typedef enum {
	END = 0,
	PLANAR,
	START
} SweepEventType;

typedef struct SweepEvent {
	Triangle* t;
	double coordinate;
	SweepEventType type;

	SweepEvent(Triangle* t, double coordinate, SweepEventType type) : t(t), coordinate(coordinate), type(type) {}

	bool operator< (SweepEvent& event) {
		return (coordinate < event.coordinate) || (coordinate == event.coordinate && type < event.type);
	}
} SweepEvent;

typedef enum { LEFT, RIGHT } SplitSide;

class KDNode {
	pair<pair<SplitPlane, SplitSide>, double> findPlane(vector<Triangle*>& triangles, AABB boundingBox);
	static KDNode* limitedBuild(vector<Triangle*>&triangles, int depth, int limit);
public:
	SplitPlane plane;
    KDNode* left;
    KDNode* right;
    vector<Triangle*> triangles;
    
    KDNode() {}
    
    static KDNode* build(vector<Triangle*>&triangles);
    
    Intersection getFirstIntersection(Ray r, double planeDist, double tMin, double tMax);
    
    set<Triangle*> getItems() {
		if (left == NULL)
		if (triangles.size() == 0) {
			return set<Triangle*>();
		}
		else return set<Triangle*>(triangles.begin(), triangles.end());
        
        set<Triangle*> leftI = left->getItems();
        set<Triangle*> rightI = right->getItems();
        leftI.insert(rightI.begin(), rightI.end());
        
        return leftI;
    }
};

#endif