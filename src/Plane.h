//An infinite plane, defined by a point on it and a normal
#ifndef PLANE_H
#define PLANE_H

#include "Vector.h"
#include "Intersection.h"
#include "Ray.h"
#include "Renderable.h"
#include <cmath>
#include <limits>

class Plane: public Renderable {
private:
    Vector point_;
    Vector normal_;
public:
    Plane(Vector point, Vector normal) {point_ = point; normal_ = normal;}
    Plane(){}

    virtual bool intersects(Ray ray);
    virtual Intersection getIntersection(Ray ray);
    virtual Vector sampleSurface();
    virtual Vector getNormalAt(Vector position);
    virtual double getSurfaceArea();
};

#endif
