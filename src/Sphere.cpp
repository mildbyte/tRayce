#include "Sphere.h"

bool Sphere::intersects(Ray ray) {
    //Assumes the ray is normalized
    double b = 2 * (ray.direction.getX() * (ray.origin.getX() - position_.getX())
        + ray.direction.getY() * (ray.origin.getY() - position_.getY())
        + ray.direction.getZ() * (ray.origin.getZ() - position_.getZ())
    );

    double c = (ray.origin.getX() - position_.getX()) * (ray.origin.getX() - position_.getX()) +
        (ray.origin.getY() - position_.getY()) * (ray.origin.getY() - position_.getY()) +
        (ray.origin.getZ() - position_.getZ()) * (ray.origin.getZ() - position_.getZ()) -
        radius_ * radius_;

    double d = b * b - 4 * c;

    //Negative delta means no intersections
    if (d < 0.0) return false;

    double t = (-b - sqrt(d)) * 0.5;
    if (t > 0.0001) return true;
    t = (-b + sqrt(d)) * 0.5;
    if (t > 0.0001) return true;
    return false;
}

Intersection Sphere::getIntersection(Ray ray) {
    //Gets data about an intersection with a ray
    Intersection result;
    result.happened = false;

    double b = 2 * (ray.direction.getX() * (ray.origin.getX() - position_.getX())
        + ray.direction.getY() * (ray.origin.getY() - position_.getY())
        + ray.direction.getZ() * (ray.origin.getZ() - position_.getZ())
    );

    double c = (ray.origin.getX() - position_.getX()) * (ray.origin.getX() - position_.getX()) +
        (ray.origin.getY() - position_.getY()) * (ray.origin.getY() - position_.getY()) +
        (ray.origin.getZ() - position_.getZ()) * (ray.origin.getZ() - position_.getZ()) -
        radius_ * radius_;

    double d = b * b - 4 * c;

    //Negative delta means no intersections
    if (d < 0) return result;

    double t = (-b - sqrt(d)) * 0.5;
    if (t > 0.001) result.happened = true;
    else {
        t = (-b + sqrt(d)) * 0.5;
        if (t > 0.001)
            result.happened = true;
        else
            return result;
    }

    result.coords = Vector(ray.origin + ray.direction * t);
    result.normal = Vector(result.coords - position_);

    //The normal should not coincide with the view direction, if it does, the
    //hit happened from the inside
    if (result.normal.dot(ray.direction) > 0) {
        result.fromTheInside = true;
    } else result.fromTheInside = false;
    result.object = this;
    result.distance = t;

    return result;
}
