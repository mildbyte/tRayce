#ifndef BOX_H
#define BOX_H

#include "Renderable.h"
#include "Vector.h"
#include "Ray.h"
#include "AABB.h"

class Box: public Renderable {
private:
    AABB _aabb;
public:
    Box() {}
    Box(Vector position, Vector size) {_aabb.point = position; _aabb.size = size;}
    virtual bool intersects(Ray ray);
    virtual Intersection getIntersection (Ray ray);
};

#endif

