#include "PhotonMap.h"
#include <limits>
#define DINFINITY std::numeric_limits<double>::infinity()

//for printf debugging
#include <cstdio>

#include <cstdlib>

//Determines which axis we are currently using to compare photons
char currAxis;
//Has to be global because it is accessed from a comparator routine.

inline double sqr(double a) {return a*a;}

double getVectorComponent(Vector a, char axis) {
    switch(axis) {
        case 0: return a.getX(); break;
        case 1: return a.getY(); break;
        case 2: return a.getZ(); break;
    }
}

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
    neighbours_ = NULL;
    neighbourDists_ = NULL;
}

PhotonMap::~PhotonMap() {
    delete(photons_);
    delete(kdTree_);
}

int photonComparator (const void* p1, const void* p2) {
    return(getVectorComponent(((Photon*) p1)->position, currAxis) >
           getVectorComponent(((Photon*) p2)->position, currAxis));
}

void PhotonMap::balanceTree(int left, int right, int depth, int position) {
    printf("left: %d, right: %d, depth: %d\n", left, right, depth);

    if (right - left <= 1) {
        kdTree_[position] = photons_[left];
        return;
    }

    currAxis = depth % 3; 
    
    printf("entering qsort...\n");
    qsort(&(photons_[left]),
          right - left, sizeof(Photon), &photonComparator);
    printf("leaving qsort...\n");
    
    Photon median = photons_[(left + right)/2];

    kdTree_[position] = median;
   
    balanceTree(left+1, (left + right)/2 + 1, depth + 1, position * 2);
    balanceTree((left + right)/2 + 1, right, depth + 1, position * 2 + 1);
}

void PhotonMap::makeTree() {
    printf("Making the kd-tree...\n");
    kdTree_ = new Photon[currPtr_ + 1];
    balanceTree(0, currPtr_, 0, 1);
}

void PhotonMap::replaceMaxDist(Photon newPh, double newPhDist) {
    if (foundNeighbours_ < neighboursNeeded_) {
        neighbours_[foundNeighbours_] = newPh;
        neighbourDists_[foundNeighbours_] = newPhDist;
        foundNeighbours_++;
        return;
    }

    int maxPos = 0;
    double maxDist = neighbourDists_[0];

    for (int i = 1; i < foundNeighbours_; i++) {
        if (neighbourDists_[i] > maxDist) {
            maxDist = neighbourDists_[i];
            maxPos = i;
        }
    }

    if (neighbourDists_[maxPos] > newPhDist) {
        neighbours_[maxPos] = newPh;
        neighbourDists_[maxPos] = newPhDist;
    }
}

inline double sqDist(Vector a, Vector b) {
    Vector diff = b - a;
    return diff.dot(diff);
}

void PhotonMap::findNearestNeighbours(Vector point, int treePos) {
    if (treePos > currPtr_ / 2) {//No children
        replaceMaxDist(kdTree_[treePos], sqDist(point, kdTree_[treePos].position));
        return;
    }

    double subdivideLocation = getVectorComponent(kdTree_[treePos].position,
                                                  kdTree_[treePos].axis);

    double distToMedian = 
        getVectorComponent(point, kdTree_[treePos].axis) - subdivideLocation;
    
    if (distToMedian < 0) {
        findNearestNeighbours(point, 2*treePos);
    } else {
        findNearestNeighbours(point, 2*treePos+1);
    }
    
    double sqDistToMedian = sqr(distToMedian);
    
    for (int i = 0; i < foundNeighbours_; i++) {
        if (sqDistToMedian < neighbourDists_[i]) {
            if (distToMedian < 0) {
                findNearestNeighbours(point, 2*treePos+1);
            } else {
                findNearestNeighbours(point, 2*treePos);
            }
            break;
        }
    }
}

void PhotonMap::nearestNeighboursWrapper(Vector point, int amount) {
    if (neighbours_ == NULL) {
        neighbours_ = new Photon[amount];
    }
    if (neighbourDists_ == NULL) {
        neighbourDists_ = new double[amount];
    }
    
    neighboursNeeded_ = amount;
    foundNeighbours_ = 0;

    findNearestNeighbours(point, 1);
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


//Gathers the photons in a given radius to determine the illumination of an entity at a certain point
//and a certain normal.
Vector PhotonMap::gatherPhotons(Vector point, Vector normal, double exposure, double sqRadius, int noPhotons) {
    Vector result(0, 0, 0);
//    double radius = getDistance(point, noPhotons);
//    sqRadius = sqr(radius);
    nearestNeighboursWrapper(point, 10);

    sqRadius = 0;
    
//    printf("found %d neighbours\n", foundNeighbours_);

    for (int i = 0; i < foundNeighbours_; i++) {
        if (neighbourDists_[i] > sqRadius) sqRadius = neighbourDists_[i];
    }

    double factor = 1.0 / 9.424777959 / sqRadius;
    factor = 10;
    double radius = sqrt(sqRadius);

    for (int i = 0; i < foundNeighbours_; i++) {
        double weight = -normal.dot(neighbours_[i].direction);
        if (weight < 0) continue;
        
        weight *= (1 - sqrt(neighbourDists_[i]) / radius);

        result += photons_[i].energy * weight;
    }

    return result * factor;
}
