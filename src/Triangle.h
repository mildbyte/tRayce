//Triangles. It begins.
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Renderable.h"
#include "Vector.h"
#include "Ray.h"
#include "Random.h"
#include "AABB.h"
#include <cstdlib>
#include <algorithm>
using namespace std;

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

	bool hasTextures = false;
	Vector t1;
	Vector t2;
	Vector t3;
    
    double area;
    AABB boundingBox;
    
    void calcBB() {
		Vector minV, maxV;

		for (int i = 0; i < 3; i++) {
			minV[i] = min(v1[i], min(v2[i], v3[i]));
			maxV[i] = max(v1[i], max(v2[i], v3[i]));
		}
        
        boundingBox = AABB(minV, maxV);
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
        area = n1.normalize() / 2;
        calcBB();
    }

	void setTextureCoords(Vector t1, Vector t2, Vector t3) {
		hasTextures = true;
		this->t1 = t1;
		this->t2 = t2;
		this->t3 = t3;
	}
    
    Vector getMidpoint() {
        return (v1 + v2 + v3) * (1.0/3.0);
    }
    
    AABB getBoundingBox() {return boundingBox;}
    
    void print() { v1.print(); printf(" -> "); v2.print(); printf(" -> "); v3.print(); }
    
    virtual bool intersects(Ray ray);
    virtual Intersection getIntersection (Ray ray);
    virtual Vector sampleSurface();
    virtual Vector getNormalAt(Vector position);
    virtual double getSurfaceArea();
	virtual bool getUVAt(Vector position, double *u, double *v);
};

#endif
