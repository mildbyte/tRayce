#include "Box.h"

bool Box::intersects(Ray ray) {
    Vector rayOrigin = ray.origin;
    Vector rayDirection = ray.direction;
    float distance[6] = {-1, -1, -1, -1, -1, -1};
    Vector boxEnd = _aabb.point + _aabb.size;



    distance[0] = (_aabb.point.getX() - rayOrigin.getX()) / ray.direction.getX();
    distance[3] = (boxEnd.getX() - rayOrigin.getX()) / ray.direction.getX();

    distance[1] = (_aabb.point.getY() - rayOrigin.getY()) / ray.direction.getY();
    distance[4] = (boxEnd.getY() - rayOrigin.getY()) / ray.direction.getY();

    distance[2] = (_aabb.point.getZ() - rayOrigin.getZ()) / ray.direction.getZ();
    distance[5] = (boxEnd.getZ() - rayOrigin.getZ()) / ray.direction.getZ();

    int closestFace = -1;
    float minDist = 0;

    for (int i = 0; i < 6; i++) {
        Vector rayToBox = rayOrigin + ray.direction * distance[i];
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

Intersection Box::getIntersection(Ray ray) {
    float distance[6] = {-1, -1, -1, -1, -1, -1};
    Vector boxEnd = _aabb.point + _aabb.size;

    distance[0] = (_aabb.point.getX() - ray.origin.getX()) / ray.direction.getX();
    distance[3] = (boxEnd.getX() - ray.origin.getX()) / ray.direction.getX();

    distance[1] = (_aabb.point.getY() - ray.origin.getY()) / ray.direction.getY();
    distance[4] = (boxEnd.getY() - ray.origin.getY()) / ray.direction.getY();

    distance[2] = (_aabb.point.getZ() - ray.origin.getZ()) / ray.direction.getZ();
    distance[5] = (boxEnd.getZ() - ray.origin.getZ()) / ray.direction.getZ();

    int closestFace = -1;
    float minDist = 0;

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

    result.fromTheInside = false;

    return result;
}
