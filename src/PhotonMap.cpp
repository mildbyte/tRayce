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
    printf("Allocating the space for the map (%d bytes per photon)...\n", sizeof(Photon));
    photons_ = new Photon[size];
    currPtr_ = 0;

    kdTreeVisited_ = 0;
}

PhotonMap::~PhotonMap() {
    delete(photons_);
    delete(kdTree_);
}

int photonComparator (const void* p1, const void* p2) {
    Photon *photon1 = (Photon*) p1;
    Photon *photon2 = (Photon*) p2;
    return(getVectorComponent(photon1->position, currAxis) >
           getVectorComponent(photon2->position, currAxis));
}

char PhotonMap::getGreatestSpreadAxis(int left, int right) {
    double greatestSpread = 0;
    char greatestSpreadAxis = 0;

    double min = DINFINITY;
    double max = 0;

    for (int i = 0; i < 3; i++) {
        min = DINFINITY;
        max = 0;

        for (int j = left; j < right; j++) {
            double position = getVectorComponent(photons_[j].position, i);
            if (position < min) min = position;
            if (position > max) max = position;
        }

        if (max - min > greatestSpread) {
            greatestSpread = max - min;
            greatestSpreadAxis = i;
        }
    }
//    printf("greatest spread for %d..%d : axis %d: %f to %f\n", left, right, greatestSpreadAxis, min, max);
    return greatestSpreadAxis;
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

    currAxis = getGreatestSpreadAxis(left, right);
    
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

/*void PhotonMap::dumpNeighbours() {
    for (int i = 0; i < foundNeighbours_; i++) {
        printf("sqd = %f; ", neighbourDists_[i]);
        dumpPhoton(photons_[neighbours_[i]]);
    }
}*/

void PhotonMap::makeTree() {
    printf("Making the kd-tree for %d photons...\n", currPtr_);

//    dumpList();
    kdTreeSize_ = 1;
    while (kdTreeSize_ < currPtr_ + 1) kdTreeSize_ = kdTreeSize_ << 1;

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

void PhotonMap::scalePhotonPower(double factor) {
    for (int i = 0; i < currPtr_; i++) {
        photons_[i].energy *= factor;
    }
}

void PhotonMap::addNearestNeighbour(int newPh, double newPhDist) {
//    printf("Adding %f to the list\n", newPhDist);
    if (neighbours_.size() < neighboursNeeded_) {
//        printf("List incomplete, added\n");
        Neighbour newNeighbour;
        newNeighbour.id = newPh;
        newNeighbour.distance = newPhDist;

        neighbours_.push(newNeighbour);
        return;
    }
    
    double maxDist = neighbours_.top().distance;

    if (maxDist > newPhDist) {
        Neighbour newNeighbour;
        newNeighbour.id = newPh;
        newNeighbour.distance = newPhDist;

        neighbours_.pop();
        neighbours_.push(newNeighbour);
//        printf("Greatest distance (%f) replaced\n", maxDist);
    }
}

inline double sqDist(Vector a, Vector b) {
    Vector diff = b - a;
    return diff.dot(diff);
}

void PhotonMap::findNearestNeighbours(Vector point, int treePos) {
//   printf("Looking for neighbours for %f, %f, %f, position %d\n", point.getX(), point.getY(), point.getZ(), treePos);
    kdTreeVisited_++;
    if (treePos >= kdTreeSize_) return;
    if (kdTree_[treePos] == -1) return;
    
    addNearestNeighbour(kdTree_[treePos], sqDist(point, photons_[kdTree_[treePos]].position));

    double subdivideLocation = getVectorComponent(photons_[kdTree_[treePos]].position,
                                                  photons_[kdTree_[treePos]].axis);

//    printf("Median ");
//    dumpPhoton(photons_[kdTree_[treePos]]);
    double distToMedian = 
        getVectorComponent(point, photons_[kdTree_[treePos]].axis) - subdivideLocation;
//    printf("Signed distance: %f\n", distToMedian);
    
//    printf("Going down the first subtree...\n");

    if (distToMedian < 0) {
        findNearestNeighbours(point, 2*treePos);
    } else {
        findNearestNeighbours(point, 2*treePos+1);
    }
//    printf("Came back\n");
    
    double sqDistToMedian = sqr(distToMedian);
    
    if (sqDistToMedian < neighbours_.top().distance) {
//            printf("Have to try the second tree: %f, %f to median\n", neighbours_.top().distance, sqDistToMedian);
        if (distToMedian < 0) {
            findNearestNeighbours(point, 2*treePos+1);
        } else {
            findNearestNeighbours(point, 2*treePos);
        }
    }
}

void PhotonMap::nearestNeighboursWrapper(Vector point, int amount) {
    neighboursNeeded_ = amount;

    findNearestNeighbours(point, 1);
}

void PhotonMap::addPhoton(Vector position, Vector direction, Vector energy) {
    Photon p;
    p.position = position; p.direction = direction; p.energy = energy;
    photons_[currPtr_++] = p;
}

double simpsonKernel(double sqx) {
    return 3.0 / PI * sqr(1 - sqx);
}


//Gathers the photons in a given radius to determine the illumination of an entity at a certain point
//and a certain normal.
Vector PhotonMap::gatherPhotons(Vector point, Vector normal, int noPhotons) {
    Vector result(0, 0, 0);
    nearestNeighboursWrapper(point, noPhotons);

    double sqRadius = neighbours_.top().distance;
    
//    double factor = 1.0 / 9.4247778 / sqRadius;
//    double factor = 3000 / 3.1415926 / sqRadius;
    double factor = 3.0 / (PI /* neighbours_.size()*/ * sqRadius);
//    factor *= 100000;
//    factor = 10;

    while (!neighbours_.empty()) { 
        Neighbour neighbour = neighbours_.top();
        neighbours_.pop();

        double weight = -normal.dot(photons_[neighbour.id].direction);
        if (weight < 0) continue;
//        weight = 1;
        
//        weight = (1 - sqrt(neighbour.distance / sqRadius));
        weight = simpsonKernel(neighbour.distance / sqRadius);

        result += photons_[neighbour.id].energy * weight;

    }

    return result * factor;
}
