#include "Plane.h"

bool Plane::intersects(Ray ray) {
    double denominator = ray.direction.dot(normal_);
    double numerator = normal_.dot(point_ - ray.origin);

    //Numerator being zero could also mean that the ray lies completely in the
    //plane
    if (std::abs(denominator) < 0.0001) return false;

    double distance = numerator / denominator;
    if (distance < 0) return false;
    return true;
}

Intersection Plane::getIntersection(Ray ray) {
    Intersection result;

    double denominator = ray.direction.dot(normal_);
    double numerator = normal_.dot(point_ - ray.origin);

    if (std::abs(denominator) < 0.0001) {
        result.happened = false;
        return result;
    }
    
    //Get the intersection distance
    result.distance = numerator / denominator;
    if (result.distance < 0) {
        result.happened = false;
        return result;
    }
    
    //Fill the intersection information
    result.happened = true;
    result.coords = ray.origin + ray.direction * result.distance;
    result.normal = normal_;
    result.object = this;

    return result;
}

Vector Plane::sampleSurface() {
    return point_;
    // sampling the plane as a light emitter for the purposes of photon mapping is problematic
}

Vector Plane::getNormalAt(Vector position) {
    return normal_;
}
    
double Plane::getSurfaceArea() {
    return 1.0; // Treat the plane as a point light
}
