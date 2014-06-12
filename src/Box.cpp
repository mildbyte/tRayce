#include "Box.h"

bool Box::intersects(Ray ray, double &dist) {
    double distance[6] = {-1, -1, -1, -1, -1, -1};
    Vector boxEnd = _aabb.point + _aabb.size;

    distance[0] = (_aabb.point.getX() - ray.origin.getX()) / ray.direction.getX();
    distance[3] = (boxEnd.getX() - ray.origin.getX()) / ray.direction.getX();

    distance[1] = (_aabb.point.getY() - ray.origin.getY()) / ray.direction.getY();
    distance[4] = (boxEnd.getY() - ray.origin.getY()) / ray.direction.getY();

    distance[2] = (_aabb.point.getZ() - ray.origin.getZ()) / ray.direction.getZ();
    distance[5] = (boxEnd.getZ() - ray.origin.getZ()) / ray.direction.getZ();

    int closestFace = -1;
    dist = 0;

    for (int i = 0; i < 6; i++) {
        Vector rayToBox = ray.origin + ray.direction * distance[i];
        if (
            (rayToBox.getX() > _aabb.point.getX() - EPSILON) &&
            (rayToBox.getX() < boxEnd.getX() + EPSILON) &&
            (rayToBox.getY() > _aabb.point.getY() - EPSILON) &&
            (rayToBox.getY() < boxEnd.getY() + EPSILON) &&
            (rayToBox.getZ() > _aabb.point.getZ() - EPSILON) &&
            (rayToBox.getZ() < boxEnd.getZ() + EPSILON)
        ) {
            if (((closestFace == -1) || distance[i] < dist) && dist > EPSILON) {
                closestFace = i;
                dist = distance[i];
            }
        }
    }
    return (closestFace != -1);
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