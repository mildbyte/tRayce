//Triangles. It begins.
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Renderable.h"
#include "Vector.h"
#include "Ray.h"
#include "Random.h"
#include <cstdlib>

class Triangle: public Renderable {
private:
    Vector v1;
    Vector v2;
    Vector v3;
    Vector e1;
    Vector e2;
    Vector normal;
    double area;
public:
    Triangle() {}
    Triangle(Vector v1, Vector v2, Vector v3)
    {
        this->v1 = v1; this->v2 = v2; this->v3 = v3;
        e1 = v2 - v1;
        e2 = v3 - v1;
        normal = e1.cross(e2);
        area = normal.normalize() / 2;
    }
    
    virtual bool intersects(Ray ray);
    virtual Intersection getIntersection (Ray ray);
    virtual Vector sampleSurface();
    virtual Vector getNormalAt(Vector position);
    virtual double getSurfaceArea();
};

#endif
