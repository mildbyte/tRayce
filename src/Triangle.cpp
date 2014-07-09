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
        //(from http://answers.unity3d.com/questions/383804/calculate-uv-coordinates-of-3d-point-on-plane-of-m.html)
        
		Vector p1 = position - v1;
		Vector p2 = position - v2;
		Vector p3 = position - v3;

		double a1 = p2.cross(p3).modulus() / area;
		double a2 = p3.cross(p1).modulus() / area;
		double a3 = p1.cross(p2).modulus() / area;

        Vector n = n1 * a1 + n2 * a2 + n3 * a3;
		n.normalize();
		return n;
    } else return n1;
}


bool Triangle::getUVAt(Vector position, double *uTex, double *vTex) {
	if (!hasTextures) return false;

	//Use the same code as for getNormalAt to get the barycentric points
	Vector p1 = position - v1;
	Vector p2 = position - v2;
	Vector p3 = position - v3;

	double a1 = p2.cross(p3).modulus() / area;
	double a2 = p3.cross(p1).modulus() / area;
	double a3 = p1.cross(p2).modulus() / area;

	*uTex = t1[0] * a1 + t2[0] * a2 + t3[0] * a3;
	*vTex = t1[1] * a1 + t2[1] * a2 + t3[1] * a3;

	return true;
}

double Triangle::getSurfaceArea()
{
    return area;
}