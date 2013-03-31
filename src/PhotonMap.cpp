#include "PhotonMap.h"
#include <limits>
#define DINFINITY std::numeric_limits<double>::infinity()

//for printf debugging
#include <cstdio>

#include <cstdlib>

//Determines which axis we are currently using to compare photons
char currAxis;
//The photons used in the comparator
Photon* globalPhotons;
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
//    delete(photons_);
//    delete(kdTree_);
}

int photonComparator (const void* p1, const void* p2) {
    Photon *photon1 = (Photon*) p1;
    Photon *photon2 = (Photon*) p2;
    return(getVectorComponent(photon1->position, currAxis) >
           getVectorComponent(photon2->position, currAxis));
}

void PhotonMap::balanceTree(int left, int right, int depth, int position) {
//    printf("left: %d, right: %d, depth: %d, position: %d\n", left, right, depth, position);
    if (right == left) return;

    if (right - left == 1) {
        kdTree_[position] = left;
//        printf("trying %d\n", left);
//        if (indices_[left]) printf("Error! %d already used!\n", left);
//        indices_[left] = 1;
//        printf("reached the base case\n");
        return;
    }

    currAxis = depth % 3; 
    
//    printf("entering qsort...\n");
    qsort(&(photons_[left]),
          right - left, sizeof(Photon), &photonComparator);
//    printf("leaving qsort...\n");
   
    int median = (left + right) / 2;
    photons_[median].axis = currAxis;

    kdTree_[position] = median;
//    printf("trying %d\n", median);
//    if (indices_[median]) printf("Error! %d already used!\n", median);

//    indices_[median] = 1;
   
    balanceTree(left, median, depth + 1, position * 2);
    balanceTree(median + 1, right, depth + 1, position * 2 + 1);
}

void dumpPhoton(Photon p) {
 printf("photon at %f, %f, %f from %f, %f, %f with %f, %f, %f, axis ",
        p.position.getX(), p.position.getY(), p.position.getZ(),
        p.direction.getX(), p.direction.getY(), p.direction.getZ(),
        p.energy.getX(), p.energy.getY(), p.energy.getZ());
    if (p.axis == 0) printf("X\n");
    if (p.axis == 1) printf("Y\n");
    if (p.axis == 2) printf("Z\n");

}

void PhotonMap::dumpList() {
    for (int i = 0; i < currPtr_; i++) {
        printf("%d ", i);
        dumpPhoton(photons_[i]);
    }
}

void PhotonMap::dumpTree() {
    for (int i = 1; i < kdTreeSize_; i++) {
        printf("%d ", i);
        if (kdTree_[i] == -1) printf("NULL\n");
        else dumpPhoton(photons_[kdTree_[i]]);
    }
}

void PhotonMap::dumpNeighbours() {
    for (int i = 0; i < foundNeighbours_; i++) {
        printf("sqd = %f; ", neighbourDists_[i]);
        dumpPhoton(photons_[neighbours_[i]]);
    }
}

void PhotonMap::makeTree() {
    printf("Making the kd-tree...\n");

//    dumpList();
    kdTreeSize_ = 1;
    while (kdTreeSize_ < currPtr_ + 1) kdTreeSize_ = kdTreeSize_ << 2;
    printf("Size %d\n", kdTreeSize_);

//    indices_ = new int[currPtr_];
//    for (int i = 0; i < currPtr_; i++) indices_[i] = 0;

    kdTree_ = new int[kdTreeSize_];
    for (int i = 0; i < kdTreeSize_; i++) kdTree_[i] = -1;
    balanceTree(0, currPtr_, 0, 1);

//    for (int i = 0; i < currPtr_; i++) if (!indices_[i])
//            printf("Error! %d wasn't touched!\n", i);

//    dumpList();
//    printf("\n\n\n");
//    dumpTree();
}

void PhotonMap::replaceMaxDist(int newPh, double newPhDist) {
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
//    printf("Looking for neighbours for %f, %f, %f, position %d\n", point.getX(), point.getY(), point.getZ(), treePos);
    if (treePos >= kdTreeSize_) return;
    if (kdTree_[treePos] == -1) return;
    
    replaceMaxDist(kdTree_[treePos], sqDist(point, photons_[kdTree_[treePos]].position));

    double subdivideLocation = getVectorComponent(photons_[kdTree_[treePos]].position,
                                                  photons_[kdTree_[treePos]].axis);

    double distToMedian = 
        getVectorComponent(point, photons_[kdTree_[treePos]].axis) - subdivideLocation;
    
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
        neighbours_ = new int[amount];
    }
    if (neighbourDists_ == NULL) {
        neighbourDists_ = new double[amount];
    }
    
    neighboursNeeded_ = amount;
    foundNeighbours_ = 0;

    findNearestNeighbours(point, 1);
}

void PhotonMap::addPhoton(Vector position, Vector direction, Vector energy) {
    Photon p;
    p.position = position; p.direction = direction; p.energy = energy;
    photons_[currPtr_++] = p;
}


//Gathers the photons in a given radius to determine the illumination of an entity at a certain point
//and a certain normal.
Vector PhotonMap::gatherPhotons(Vector point, Vector normal, int noPhotons) {
    Vector result(0, 0, 0);
//    double radius = getDistance(point, noPhotons);
//    sqRadius = sqr(radius);
    nearestNeighboursWrapper(point, noPhotons);

    double sqRadius = 0;
    
//printf("Looking for neighbours for %f, %f, %f\n\n\n", point.getX(), point.getY(), point.getZ());
//dumpNeighbours();
//    printf("found %d neighbours\n", foundNeighbours_);
/*
    double minDist = DINFINITY;
    int besti = 0;
    for (int i = 0; i < currPtr_; i++) {
        double currDist = sqDist(point, photons_[i].position);
        if (currDist < minDist) {minDist = currDist; besti = i;}
    }

    if (minDist < neighbourDists_[0]) {
        printf("ERROR!\n");
        printf("Brute force result: %f; ", minDist);
        dumpPhoton(photons_[besti]);

        for (int i = 1; i < currPtr_+1; i++) {
            if (kdTree_[i] == besti) {
                printf("Found in the tree at %d (in the photons array at %d", i, besti);
                break;
            }
        }
    }
*/
    for (int i = 0; i < foundNeighbours_; i++) {
        if (neighbourDists_[i] > sqRadius) sqRadius = neighbourDists_[i];
    }

//    double factor = 1.0 / 9.424777959 / sqRadius;
//    double factor = 1.0 / sqRadius;
    double factor = 1;
    double radius = sqrt(sqRadius);

    for (int i = 0; i < foundNeighbours_; i++) {
        double weight = -normal.dot(photons_[neighbours_[i]].direction);
        if (weight < 0) continue;
        
        weight *= (1 - sqrt(neighbourDists_[i]) / radius);

        result += photons_[neighbours_[i]].energy * weight;
    }

    return result * factor;
}
