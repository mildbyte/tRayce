#ifndef AABB_H
#define AABB_H

#include "Vector.h"
#include "Ray.h"
#include <algorithm>
using namespace std;

class AABB {
private:
    Vector start;
    Vector end;
public:
	AABB() {}
	AABB(Vector startPoint, Vector endPoint) : start(startPoint), end(endPoint) {}
	bool intersects(Ray ray, double &dist);
	void addBox(AABB b);
	int getGreatestSpread();

	Vector getStartpoint() { return start; }
	Vector getSize() { return end - start; }
	Vector getEndpoint() { return end; }
};

#endif
