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
    
    bool customNormals = false;
    Vector n1; //if !customNormals, holds the normal derived via cross product
    Vector n2;
    Vector n3;
    
    double area;
public:
    Triangle() {}
    Triangle(Vector v1, Vector v2, Vector v3, Vector n1, Vector n2, Vector n3) {
        this->v1 = v1; this->v2 = v2; this->v3 = v3;
        this->n1 = n1; this->n2 = n2; this->n3 = n3;
        customNormals = true;
        
        e1 = v2 - v1;
        e2 = v3 - v1;
        area = e1.cross(e2).normalize() / 2;
    }
    
    Triangle(Vector v1, Vector v2, Vector v3)
    {
        this->v1 = v1; this->v2 = v2; this->v3 = v3;
        e1 = v2 - v1;
        e2 = v3 - v1;
        n1 = e1.cross(e2);
        area = n2.normalize() / 2;
    }
    
    virtual bool intersects(Ray ray);
    virtual Intersection getIntersection (Ray ray);
    virtual Vector sampleSurface();
    virtual Vector getNormalAt(Vector position);
    virtual double getSurfaceArea();
};

#endif
