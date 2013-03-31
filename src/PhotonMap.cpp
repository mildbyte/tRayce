#include "PhotonMap.h"
#include <limits>
#define DINFINITY std::numeric_limits<double>::infinity()

//for printf debugging
#include <cstdio>

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
//    printf("OH HAI photon at %f, %f, %f from %f, %f, %f with %f, %f, %f!\n",
//        position.getX(), position.getY(), position.getZ(),
//        direction.getX(), direction.getY(), direction.getZ(),
//        energy.getX(), energy.getY(), energy.getZ());
    Photon p;
    p.position = position; p.direction = direction; p.energy = energy;
    photons_[currPtr_++] = p;
}

int PhotonMap::countPhotonsAt(Vector point, double sqDist) {
    int result = 0;
    for (int i = 0; i < currPtr_; i++) {
           if (isWithin(point, photons_[i].position, sqDist) != DINFINITY) {
               result++;
           }
    }
    return result;
}

int abs (int a) {return a < 0 ? -a : a;}

double PhotonMap::getDistance(Vector point, int noPhotons) {
    double left = 0;
    double right = 10000;

    while (true) {
        double pivot = (left + right) / 2;
//        printf("left=%f, right=%f\n", left, right);
        int currNoPhotons = countPhotonsAt(point, sqr(pivot));

//        printf("curr=%d, needed=%d\n", currNoPhotons, noPhotons);

        if (abs(noPhotons - currNoPhotons) < 2) {
            return pivot;
        } else if (currNoPhotons < noPhotons) {
            left = pivot;
        } else right = pivot;
    }
}

//Gathers the photons in a given radius to determine the illumination of an entity at a certain point
//and a certain normal.
Vector PhotonMap::gatherPhotons(Vector point, Vector normal, double exposure, double sqRadius, int noPhotons) {
    Vector result(0, 0, 0);
    double radius = getDistance(point, noPhotons);
    sqRadius = sqr(radius);
    double factor = 1/3.1415926 / sqRadius;

    for (int i = 0; i < currPtr_; i++) {
        double distancesq = isWithin(point, photons_[i].position, sqRadius);

        if (distancesq == DINFINITY) continue;

//        double weight = -normal.dot(photons_[i].direction);
//        if (weight < 0) continue;
        
//        weight *= (radius - sqrt(distancesq)) / exposure / radius;
        result += photons_[i].energy * factor;
    }

    return result;
}
