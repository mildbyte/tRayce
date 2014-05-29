//Data about the image plane
#ifndef CAMERA_H
#define CAMERA_H

#include "Vector.h"

struct Camera {
    Vector position;
    Vector direction;
    double width, height;
    double planeDistance;
    double lensRadius; //Ray source randomly perturbed
};

#endif
