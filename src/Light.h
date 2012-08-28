//Stores data about a light (currently only point lights)
#ifndef LIGHT_H
#define LIGHT_H

#include "Vector.h"
#include "Renderable.h"

enum LightType {
    LT_NONE,
    LT_POINT,
    LT_AREA
};

struct Light {
    LightType type;
    Vector position;
    Light(){type = LT_NONE; brightness = 1; position.set(0, 0, 0);}
    double brightness;
};

struct PointLight: public Light {
    PointLight(){type = LT_POINT;}
};

struct AreaLight: public Light {
    //Area light (array of several point lights, can cast soft shadows)
    //Inspired by PovRay (http://www.povray.org/documentation/view/3.6.0/313/)
    //Two direction vectors and their lengths
    Vector dir1, dir2;
    double size1, size2;

    AreaLight(){type = LT_AREA; dir1.set(1, 0, 0); dir2.set(0, 0, -1);
        size1 = 1; size2 = 1;}
};

#endif
