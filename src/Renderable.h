//An abstract class for all items than can be visible on the screen (spheres etc)
#ifndef RENDERABLE_H
#define RENDERABLE_H

#include "Intersection.h"
#include "Ray.h"
#include "Material.h"

class Renderable {
private:
    //Material material_;
public:
    Material material;
    //Material& getMaterial() {return material_;}

    virtual bool intersects(Ray ray) = 0;
    virtual Intersection getIntersection(Ray ray) = 0;
    virtual Vector sampleSurface() = 0;
    virtual Vector getNormalAt(Vector position) = 0;
    
    //Surface area for scaling in pathtracing
    virtual double getSurfaceArea() = 0;
};

#endif
