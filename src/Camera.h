//Data about the image plane
#ifndef CAMERA_H
#define CAMERA_H

#include "Vector.h"

struct Camera {
    Vector position;
    Vector direction;
    float width, height;
    float planeDistance;
};

#endif
