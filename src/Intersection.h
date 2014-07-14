//Intersection structure
//Stores data about an intersection.
#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "Vector.h"

class Renderable;

struct Intersection {
    bool happened;
    Vector coords;
    Vector normal;
    double distance;
    Renderable* object;
	int intersectedTriangles;
};

#endif
