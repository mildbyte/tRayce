#include "PhotonMap.h"
#include <limits>
#define DINFINITY std::numeric_limits<double>::infinity()

inline double sqr(double a) {return a*a;}

//Determines if two points are within a certain distance squared within each other
//returns the squared distance if they are
//returns the infinity if they aren't
double isWithin(Vector a, Vector b, double distsq) {
    double curr = sqr(a.getX() - b.getX());
    if (curr > distsq) return DINFINITY;

    curr += sqr(a.getY() - b.getY());
    if (curr > distsq) return DINFINITY;
    
    curr += sqr(a.getZ() - b.getZ());
    if (curr > distsq) return DINFINITY;

    return curr;
}

PhotonMap::PhotonMap(int size) {
    photons_ = new Photon[size];
    currPtr_ = 0;
}

PhotonMap::~PhotonMap() {
    delete(photons_);
}

void PhotonMap::addPhoton(Vector position, Vector direction, Vector energy) {
    Photon p;
    p.position = position; p.direction = direction; p.energy = energy;
    photons_[currPtr_++] = p;
}

//Gathers the photons in a given radius to determine the illumination of an entity at a certain point
//and a certain normal.
Vector PhotonMap::gatherPhotons(Vector point, Vector normal, double exposure, double sqRadius) {
    Vector result(0, 0, 0);

    for (int i = 0; i < currPtr_; i++) {
        double distancesq = isWithin(point, photons_[i].position, sqRadius);
        if (distancesq == DINFINITY) continue;

        double weight = -normal.dot(photons_[i].direction);
        if (weight < 0) continue;
        
        weight *= (1.0 - sqrt(distancesq)) / exposure;
        result += photons_[i].energy * weight;
    }

    return result;
}
