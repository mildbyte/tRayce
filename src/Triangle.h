//Triangles. It begins.
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Renderable.h"
#include "Vector.h"
#include "Ray.h"
#include "Random.h"
#include "Box.h"
#include <cstdlib>
#include <algorithm>

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
    Box boundingBox;
    
    void calcBB() {
        Vector minV = v1;
        
        minV.setX(min(minV.getX(), v2.getX()));
        minV.setY(min(minV.getY(), v2.getY()));
        minV.setZ(min(minV.getZ(), v2.getZ()));
        
        minV.setX(min(minV.getX(), v3.getX()));
        minV.setY(min(minV.getY(), v3.getY()));
        minV.setZ(min(minV.getZ(), v3.getZ()));
        
        Vector maxV = v2;
        
        maxV.setX(max(maxV.getX(), v2.getX()));
        maxV.setY(max(maxV.getY(), v2.getY()));
        maxV.setZ(max(maxV.getZ(), v2.getZ()));
        
        maxV.setX(max(maxV.getX(), v3.getX()));
        maxV.setY(max(maxV.getY(), v3.getY()));
        maxV.setZ(max(maxV.getZ(), v3.getZ()));
        
        boundingBox = Box(minV, maxV - minV);
    }
    
public:
    Triangle() {}
    Triangle(Vector v1, Vector v2, Vector v3, Vector n1, Vector n2, Vector n3) {
        this->v1 = v1; this->v2 = v2; this->v3 = v3;
        this->n1 = n1; this->n2 = n2; this->n3 = n3;
        customNormals = true;
        
        e1 = v2 - v1;
        e2 = v3 - v1;
        area = e1.cross(e2).normalize() / 2;
        calcBB();
    }
    
    Triangle(Vector v1, Vector v2, Vector v3)
    {
        this->v1 = v1; this->v2 = v2; this->v3 = v3;
        e1 = v2 - v1;
        e2 = v3 - v1;
        n1 = e1.cross(e2);
        area = n2.normalize() / 2;
        calcBB();
    }
    
    Vector getMidpoint() {
        return (v1 + v2 + v3) * (1.0/3.0);
    }
    
    Box getBoundingBox() {return boundingBox;}
    
    virtual bool intersects(Ray ray);
    virtual Intersection getIntersection (Ray ray);
    virtual Vector sampleSurface();
    virtual Vector getNormalAt(Vector position);
    virtual double getSurfaceArea();
};

#endif
