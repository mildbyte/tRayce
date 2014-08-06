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
    virtual Vector getNormalAt(Vector position) = 0;

	//Texture mapping
	virtual bool getUVAt(Vector, double*, double*) { return false; } //does not support by default
};

#endif
