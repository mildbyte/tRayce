#include "Triangle.h"
    
bool Triangle::intersects(Ray ray)
{
    Vector P = ray.direction.cross(e2);
    double det = P.dot(e1);
    
    if (det >= -EPSILON && det <= EPSILON) return false;
    
    double invDet = 1.0/det;
    
    Vector T = ray.origin - v1;
    double u = T.dot(P) * invDet;
    if (u < EPSILON || u > 1 + EPSILON) return false;
    
    Vector Q = T.cross(e1);
    double v = ray.direction.dot(Q) * invDet;
    if (u < EPSILON || u + v > 1 + EPSILON) return false;
    
    double t = e2.dot(Q) * invDet;
    if (t > EPSILON) return true;
    
    return false;
}



// Implements the Möller–Trumbore ray-triangle intersection algorithm.
Intersection Triangle::getIntersection (Ray ray)
{
    Intersection result;
    result.happened = false;
    
    Vector P = ray.direction.cross(e2);
    double det = P.dot(e1);
    
    if (det >= -EPSILON && det <= EPSILON) return result;
    
    double invDet = 1.0/det;
    
    Vector T = ray.origin - v1;
    double u = T.dot(P) * invDet;
    if (u < EPSILON || u > 1 + EPSILON) return result;
    
    Vector Q = T.cross(e1);
    double v = ray.direction.dot(Q) * invDet;
    if (v < EPSILON || u + v > 1 + EPSILON) return result;
    
    double t = e2.dot(Q) * invDet;
    if (t > EPSILON) {
        result.happened = true;
        result.object = this;
        
		if (customNormals) result.normal = n1 * (1 - u - v) + n2 * u + n3 * v;
        else result.normal = n1;
            
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
    if (customNormals) {
        //Interpolate the normal
        
        //Get the barycentric coordinates of the point
        //(from https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates)
        
        Vector relPos = position - v1;
        
        double d00 = e1.dot(e1);
        double d01 = e1.dot(e2);
        double d11 = e2.dot(e2);
        double d20 = relPos.dot(e1);
        double d21 = relPos.dot(e2);
        double denom = d00 * d11 - d01 * d01;
        double v = (d11 * d20 - d01 * d21) / denom;
        double w = (d00 * d21 - d01 * d20) / denom;
        double u = 1.0 - v - w;

        Vector n = n1 * w + n2 * u + n3 * v;
		n.normalize();
		return n;
    } else return n1;
}

double Triangle::getSurfaceArea()
{
    return area;
}