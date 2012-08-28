//A sphere - the simplest object for raytracing
#ifndef SPHERE_H
#define SPHERE_H

#include "Renderable.h"
#include "Vector.h"
#include "Ray.h"

class Sphere: public Renderable {
private:
    double radius_;
    Vector position_;
public:
    Sphere() {}
    Sphere(double radius, Vector position) {radius_ = radius; position_ = position;}
    virtual bool intersects(Ray ray);
    virtual Intersection getIntersection (Ray ray);
};

#endif
