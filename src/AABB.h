#ifndef AABB_H
#define AABB_H

#include "Vector.h"
#include "Ray.h"
#include <algorithm>
using namespace std;

class SplitPlane {
private:
	double coordinate;
	int axis;
public:
	double getCoordinate() { return coordinate; }
	int getAxis() { return axis; }
	SplitPlane(double coordinate, int axis) : coordinate(coordinate), axis(axis) {}
	SplitPlane() {}
};

class AABB {
private:
    Vector start;
    Vector end;
	bool empty;
	bool planar;

	Vector size;
public:
	AABB() {}
	AABB(Vector startPoint, Vector endPoint) : start(startPoint), end(endPoint) {
		empty = false;
		planar = false;
		for (int i = 0; i < 3; i++) {
			if (startPoint[i] > endPoint[i]) empty = true;
			if (startPoint[i] == endPoint[i]) planar = true;
			size[i] = max(0.0, endPoint[i] - startPoint[i]);
		}
	}
	bool intersects(Ray ray, double &dist);
	void addBox(AABB b);
	int getGreatestSpread();

	bool isEmpty() {
		return empty;
	}

	bool isPlanar() {
		return planar;
	}

	Vector getStartpoint() { return start; }
	Vector getSize() { return size; }
	Vector getEndpoint() { return end; }

	AABB intersection(AABB box) {
		Vector start;
		Vector end;

		for (int i = 0; i < 3; i++) {
			start[i] = max(this->start[i], box.start[i]);
			end[i] = min(this->end[i], box.end[i]);
		}

		//No intersection if any start[i] > end[i]
		return AABB(start, end);
	}

	pair<AABB, AABB> split(SplitPlane plane) {
		Vector end1 = end;
		Vector start2 = start;
		end1[plane.getAxis()] = plane.getCoordinate();
		start2[plane.getAxis()] = plane.getCoordinate();

		return make_pair(AABB(start, end1), AABB(start2, end));
	}
};

#endif
