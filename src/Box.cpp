#include "Box.h"

bool Box::intersects(Ray ray) {
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
        if (
            (rayToBox.getX() > _aabb.point.getX() - 0.001) &&
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
    Vector difference = b._aabb.point - _aabb.point;
    difference.setX(min(difference.getX(), 0.0));
    difference.setY(min(difference.getY(), 0.0));
    difference.setZ(min(difference.getZ(), 0.0));
    
    _aabb.point += difference;
    
    difference = b._aabb.point + b._aabb.size - _aabb.point - _aabb.size;
    difference.setX(max(difference.getX(), 0.0));
    difference.setY(max(difference.getY(), 0.0));
    difference.setZ(max(difference.getZ(), 0.0));
    
    _aabb.size += difference;
}