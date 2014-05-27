#include "Triangle.h"
    
bool Triangle::intersects(Ray ray)
{
    Vector P = ray.direction.cross(e2);
    double det = P.dot(e1);
    
    if (det > -0.0001 && det < 0.0001) return false;
    
    double invDet = 1.0/det;
    
    Vector T = ray.origin - v1;
    double u = T.dot(P) * invDet;
    if (u < 0 || u > 1) return false;
    
    Vector Q = T.cross(e1);
    double v = ray.direction.dot(Q) * invDet;
    if (u < 0 || u + v > 1) return false;
    
    double t = e2.dot(Q) * invDet;
    if (t > 0.0001) return true;
    
    return false;
}

// Implements the Möller–Trumbore ray-triangle intersection algorithm.
Intersection Triangle::getIntersection (Ray ray)
{
    Intersection result;
    result.happened = false;
    
    Vector P = ray.direction.cross(e2);
    double det = P.dot(e1);
    
    if (det > -0.0001 && det < 0.0001) return result;
    
    double invDet = 1.0/det;
    
    Vector T = ray.origin - v1;
    double u = T.dot(P) * invDet;
    if (u < 0 || u > 1) return result;
    
    Vector Q = T.cross(e1);
    double v = ray.direction.dot(Q) * invDet;
    if (v < 0 || u + v > 1) return result;
    
    double t = e2.dot(Q) * invDet;
    if (t > 0.0001) {
        result.happened = true;
        result.object = this;
        result.normal = normal;
        result.distance = t;
        result.coords = ray.origin + ray.direction * t;
    }
    
    return result;
}

Vector Triangle::sampleSurface()
{
    double u;
    double v;
    
    do {
        u = drand();
        v = drand();
    } while (u + v > 1.0);
    
    return v1 * (1 - u - v) + v2 * u + v3 * v;
}

Vector Triangle::getNormalAt(Vector position)
{
    return normal;
}

double Triangle::getSurfaceArea()
{
    return area;
}