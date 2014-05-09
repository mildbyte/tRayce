//A list of all objects on scene
#ifndef RENDERABLESLIST_H
#define RENDERABLESLIST_H

#include <list>
#include "Renderable.h"
#include "Intersection.h"
#include <cmath>

class RenderablesList {
private:
    std::list<Renderable*> rendlist_;
public:
    RenderablesList();
    ~RenderablesList();

    //Gets the closest intersection to the origin
    Intersection getFirstIntersection(Ray ray, double planeDistance);

    void addRenderable(Renderable* renderable);

    //Finds out if the ray intersects any object on the scene
    bool intersects(Ray ray);

    //Does the ray intersect an object closer than a certain value?
    bool intersectsCloser(Ray ray, double distance);
};

#endif
