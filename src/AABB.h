//Axis-aligned bounding box, used to subdivide a scene into a grid

#ifndef AABB_H
#define AABB_H

struct AABB {
    Vector point;
    Vector size;
    AABB() {point.set(0, 0, 0); size.set(0, 0, 0);}
    AABB(Vector point, Vector size) {this->point = point; this->size = size;}
};

#endif
