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
    bool intersects(Ray ray, double &dist);
    //Intersection getIntersection (Ray ray);
    
    void addBox(Box b);
    
    int getGreatestSpread() {
        double size = 0;
        int axis = -1;
        for (int i = 0; i < 3; i++) {
            if (axis == -1 || _aabb.size[i] > size) {
                axis = i; size = _aabb.size[i];
            }
        }
        
        return axis;
    }
    
    Vector getPosition() {return _aabb.point;}
    Vector getSize() {return _aabb.size;}
    Vector getEndpoint() {return _aabb.point + _aabb.size;}
};

#endif

