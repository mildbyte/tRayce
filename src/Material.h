//A primitive's color, diffuse, specular, transparency etc.
#ifndef MATERIAL_H
#define MATERIAL_H

#include "Vector.h"

struct Material {
    Vector color;
    double diffuse;
    double specular;
    double specularExp;
    double transparency;
    double refrIndex;
    double reflectivity;
    Vector emittance;

    Material() {
        color.set(0, 0, 0);
        emittance.set(0, 0, 0);
        diffuse = 1;
        specular = 1;
        specularExp = 20;
        transparency = 0.0;
        reflectivity = 0.0;
        refrIndex = 1.0;
    }
};

#endif
