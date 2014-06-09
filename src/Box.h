#ifndef BOX_H
#define BOX_H

#include "Renderable.h"
#include "Vector.h"
#include "Ray.h"
#include "AABB.h"
#include <algorithm>
using namespace std;

class Box {
private:
    AABB _aabb;
public:
    Box() {}
    Box(Vector position, Vector size) {_aabb.point = position; _aabb.size = size;}
    bool intersects(Ray ray);
    //Intersection getIntersection (Ray ray);
    
    void addBox(Box b);
    
    int getGreatestSpread() {
        if (_aabb.size.getX() >= _aabb.size.getY() && _aabb.size.getX() >= _aabb.size.getZ()) return 0;
        else if (_aabb.size.getY() >= _aabb.size.getX() && _aabb.size.getY() >= _aabb.size.getZ()) return 1;
        else return 2;
    }
    
    Vector getPosition() {return _aabb.point;}
    Vector getSize() {return _aabb.size;}
    Vector getEndpoint() {return _aabb.point + _aabb.size;}
};

#endif

