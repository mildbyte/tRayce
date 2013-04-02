#include "RenderablesList.h"

RenderablesList::RenderablesList() {
    rendlist_.clear();
}

RenderablesList::~RenderablesList() {
    std::list<Renderable*>::iterator it = rendlist_.begin();
    while(it != rendlist_.end()) {
        delete (*it);
        it++;
    }
    rendlist_.clear();
}

bool RenderablesList::intersects(Ray ray) {
    std::list<Renderable*>::iterator it = rendlist_.begin();
    while(it != rendlist_.end()) {
        if ((*it)->intersects(ray)) return true;
        it++;
    }
    return false;
}

Intersection RenderablesList::getFirstIntersection(Ray ray, double planeDistance) {
    //Main bottleneck. Works with O(n), can be modified to reach O(logn)
    //Linear search to find the closest intersection
    Intersection bestInter;
    Intersection currInter;

    double mindist = 0;
    bool found = false;

    std::list<Renderable*>::iterator it = rendlist_.begin();
    while(it != rendlist_.end()) {
        currInter = (*it)->getIntersection(ray);
        it++;
        if (currInter.happened) { 
            //Ignore hits that happened before the image plane
            if (currInter.distance < planeDistance) continue;
            if (!found || (currInter.distance < mindist)) {
                bestInter = currInter;
                mindist = currInter.distance;
                found = true;
            }
        }
    }

    if (!found) bestInter.happened = false;

    return bestInter;
}

bool RenderablesList::intersectsCloser(Ray ray, double distance) {
    //Checks if the ray intersects an object closer that at a certain distance.
    //Used to check if a point on an object is shaded (faster than checking
    //every object on the scene.
    Intersection currInter;

    std::list<Renderable*>::iterator it = rendlist_.begin();
    while(it != rendlist_.end()) {
        currInter = (*it)->getIntersection(ray);
        if (currInter.happened) {
           if (currInter.distance <= distance && currInter.distance > 0) return true;
        }
        it++;
    }

    return false;
}

void RenderablesList::addRenderable(Renderable* renderable) {
    rendlist_.push_back(renderable);
}
