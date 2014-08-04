#include "AABB.h"

int AABB::getGreatestSpread() {
	Vector boxSize = end - start;
	double size = 0;
	int axis = -1;
	for (int i = 0; i < 3; i++) {
		if (axis == -1 || boxSize[i] > size) {
			axis = i; size = boxSize[i];
		}
	}

	return axis;
}

bool AABB::intersects(Ray ray, double &dist) {
    double distance[6] = {-1, -1, -1, -1, -1, -1};

    distance[0] = (start.getX() - ray.origin.getX()) / ray.direction.getX();
    distance[3] = (end.getX() - ray.origin.getX()) / ray.direction.getX();

    distance[1] = (start.getY() - ray.origin.getY()) / ray.direction.getY();
    distance[4] = (end.getY() - ray.origin.getY()) / ray.direction.getY();

    distance[2] = (start.getZ() - ray.origin.getZ()) / ray.direction.getZ();
    distance[5] = (end.getZ() - ray.origin.getZ()) / ray.direction.getZ();

    int closestFace = -1;
    dist = 0;

    for (int i = 0; i < 6; i++) {
        Vector rayToBox = ray.origin + ray.direction * distance[i];
        if (
            (rayToBox.getX() > start.getX() - EPSILON) &&
            (rayToBox.getX() < end.getX() + EPSILON) &&
            (rayToBox.getY() > start.getY() - EPSILON) &&
            (rayToBox.getY() < end.getY() + EPSILON) &&
            (rayToBox.getZ() > start.getZ() - EPSILON) &&
            (rayToBox.getZ() < end.getZ() + EPSILON)
        ) {
            if (((closestFace == -1) || distance[i] < dist) && distance[i] > EPSILON) {
                closestFace = i;
                dist = distance[i];
            }
        }
    }
    return (closestFace != -1);
}

void AABB::addBox(AABB b) {
	for (int i = 0; i < 3; i++) {
		start[i] = min(start[i], b.start[i]);
		end[i] = max(end[i], b.end[i]);
		size[i] = max(0.0, end[i] - start[i]);
	}

}