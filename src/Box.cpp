#include "Box.h"

bool Box::intersects(Ray ray, double &dist) {
    Vector boxEnd = _aabb.point + _aabb.size;

    double t1 = (_aabb.point.getX() - ray.origin.getX()) / ray.direction.getX();
    double t2 = (boxEnd.getX() - ray.origin.getX()) / ray.direction.getX();

    double t3 = (_aabb.point.getY() - ray.origin.getY()) / ray.direction.getY();
    double t4 = (boxEnd.getY() - ray.origin.getY()) / ray.direction.getY();

    double t5 = (_aabb.point.getZ() - ray.origin.getZ()) / ray.direction.getZ();
    double t6 = (boxEnd.getZ() - ray.origin.getZ()) / ray.direction.getZ();

    double tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
    double tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));
    
    if (tmax < 0) return false;
    if (tmin > tmax) return false;
    
    dist = tmin;
    return true;
}
/*
Intersection Box::getIntersection(Ray ray) {
    double distance[6] = {-1, -1, -1, -1, -1, -1};
    Vector boxEnd = _aabb.point + _aabb.size;

    distance[0] = (_aabb.point.getX() - ray.origin.getX()) / ray.direction.getX();
    distance[3] = (boxEnd.getX() - ray.origin.getX()) / ray.direction.getX();

    distance[1] = (_aabb.point.getY() - ray.origin.getY()) / ray.direction.getY();
    distance[4] = (boxEnd.getY() - ray.origin.getY()) / ray.direction.getY();

    distance[2] = (_aabb.point.getZ() - ray.origin.getZ()) / ray.direction.getZ();
    distance[5] = (boxEnd.getZ() - ray.origin.getZ()) / ray.direction.getZ();

    int closestFace = -1;
    double minDist = 0;

    for (int i = 0; i < 6; i++) {
        Vector rayToBox = ray.origin + ray.direction * distance[i];
        if ((rayToBox.getX() > _aabb.point.getX() - 0.001) &&
            (rayToBox.getX() < boxEnd.getX() + 0.001) &&
            (rayToBox.getY() > _aabb.point.getY() - 0.001) &&
            (rayToBox.getY() < boxEnd.getY() + 0.001) &&
            (rayToBox.getZ() > _aabb.point.getZ() - 0.001) &&
            (rayToBox.getZ() < boxEnd.getZ() + 0.001)
        ) {
            if ((closestFace == -1) || distance[i] < minDist) {
                closestFace = i;
                minDist = distance[i];
            }
        }
    }

    Intersection result;
    result.happened = (closestFace != -1);
    if (!result.happened) return result;

    result.distance = minDist;
    result.coords = ray.origin + ray.direction * (minDist + 0.0001);
    result.object = this;

    switch (closestFace) {
        case 0: result.normal = Vector(-1, 0, 0); break;
        case 1: result.normal = Vector(1, 0, 0); break;
        case 2: result.normal = Vector(0, -1, 0); break;
        case 3: result.normal = Vector(0, 1, 0); break;
        case 4: result.normal = Vector(0, 0, -1); break;
        case 5: result.normal = Vector(0, 0, 1); break;
    }

    return result;
}*/

void Box::addBox(Box b) {
    /*printf("Source box: ");
    getPosition().print();
    printf(" -> ");
    getEndpoint().print();
    printf("; added ");
    b.getPosition().print();
    printf(" -> ");
    b.getEndpoint().print();
    printf(" get ");*/
    
    Vector end = _aabb.point + _aabb.size;
    
    for (int i = 0; i < 3; i++) _aabb.point[i] = min(_aabb.point[i], b._aabb.point[i]);
    
    Vector bEnd = b._aabb.point + b._aabb.size;
    
    for (int i = 0; i < 3; i++) end[i] = max(end[i], bEnd[i]);
    
    _aabb.size = end - _aabb.point;
    
    /*
    getPosition().print();
    printf(" -> ");
    getEndpoint().print();
    printf("\n");*/
}