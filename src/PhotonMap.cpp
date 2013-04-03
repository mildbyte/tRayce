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
    kdTreeVisited_++;
    if (treePos >= kdTreeSize_) return;
    if (kdTree_[treePos] == -1) return;
    
    addNearestNeighbour(kdTree_[treePos], sqDist(point, photons_[kdTree_[treePos]].position));

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
    
    if (sqDistToMedian < neighbours_.top().distance) {
        if (distToMedian < 0) {
            findNearestNeighbours(point, 2*treePos+1);
        } else {
            findNearestNeighbours(point, 2*treePos);
        }
    }
}

void PhotonMap::nearestNeighboursWrapper(Vector point, int amount) {
    neighboursNeeded_ = amount;
    neighbours_ = std::priority_queue<Neighbour>();

    findNearestNeighbours(point, 1);
}

void PhotonMap::addPhoton(Vector position, Vector direction, Vector energy, Vector normal) {
    Photon p;
    p.position = position; p.direction = direction; p.energy = energy; p.normal = normal;
    photons_[currPtr_++] = p;
}

double simpsonKernel(double sqx) {
    return 3.0 / PI * sqr(1 - sqx);
}

void PhotonMap::precalculateIrradiance(int frequency, int noPhotons) {
    irradiancePhotonFrequency_ = frequency;

    //Precompute the radiance for every frequencyth photon
    for (int i = 0; i < currPtr_; i+= frequency) {
        if (i % 1000 == 0) {
            printf ("\r%d out of %d photons processed (%f\%)\n", 
                i, currPtr_, (double)i/currPtr_*100);
        }
        photons_[i].irradiance = irradianceEstimate(
            photons_[i].position, photons_[i].normal, noPhotons);
    }
}

void PhotonMap::findIrradiancePhoton(Vector point, Vector normal,
                                     double threshold, int treePos) {
     
    if (treePos >= kdTreeSize_) return;
    if (kdTree_[treePos] == -1) return;

    //Set the photon as best if it's indeed an irradiance photon, its normal is close to
    //our needed normal and it's closer to the sought point than the previous best
    if (kdTree_[treePos] % irradiancePhotonFrequency_ == 0) {
        if (normal.dot(photons_[kdTree_[treePos]].normal) >= threshold) {
            double sqdistance = sqDist(point, photons_[kdTree_[treePos]].position);
            if (sqdistance < irradiancePhotonDist_) {
                irradiancePhotonDist_ = sqdistance;
                irradiancePhotonId_ = kdTree_[treePos];
            }
        }
    }
    
    double subdivideLocation = getVectorComponent(photons_[kdTree_[treePos]].position,
                                                  photons_[kdTree_[treePos]].axis);

    double distToMedian = 
        getVectorComponent(point, photons_[kdTree_[treePos]].axis) - subdivideLocation;
    
    if (distToMedian < 0) {
        findIrradiancePhoton(point, normal, threshold, 2*treePos);
    } else {
        findIrradiancePhoton(point, normal, threshold, 2*treePos + 1);
    }
    
    double sqDistToMedian = sqr(distToMedian);
    
    if (sqDistToMedian < irradiancePhotonDist_) {
        if (distToMedian < 0) {
            findIrradiancePhoton(point, normal, threshold, 2*treePos + 1);
        } else {
            findIrradiancePhoton(point, normal, threshold, 2*treePos);
        }
    }

//            findIrradiancePhoton(point, normal, threshold, 2*treePos);
//            findIrradiancePhoton(point, normal, threshold, 2*treePos + 1);
}

Vector PhotonMap::acceleratedIrradiance(Vector point, Vector normal, double threshold) {
    irradiancePhotonDist_ = DINFINITY;
    irradiancePhotonId_ = -1;

    findIrradiancePhoton(point, normal, threshold, 1);

    if (irradiancePhotonId_ == -1) return Vector(0, 0, 0);

    return photons_[irradiancePhotonId_].irradiance;
}

//Gathers the photons in a given radius to determine the illumination of an entity at a certain point
//and a certain normal.
Vector PhotonMap::irradianceEstimate(Vector point, Vector normal, int noPhotons) {
    Vector result(0, 0, 0);
    nearestNeighboursWrapper(point, noPhotons);

    double sqRadius = neighbours_.top().distance;
    
    double factor = 3.0 / (PI * sqRadius);
    while (!neighbours_.empty()) { 
        Neighbour neighbour = neighbours_.top();
        neighbours_.pop();

//        double weight = -normal.dot(photons_[neighbour.id].direction);
//        if (weight < 0) continue;
        
//        weight = (1 - sqrt(neighbour.distance / sqRadius));
        double weight = simpsonKernel(neighbour.distance / sqRadius);

        result += photons_[neighbour.id].energy * weight;

    }
    return result * factor;
}

//Returns a pixel on the gradient from black to white depending on the distance of the closest
//photon. The smaller the weight, the faster the falloff.
Vector PhotonMap::visualizePhoton(Vector point, double weight) {
    nearestNeighboursWrapper(point, 1);
    double distance = neighbours_.top().distance;
    double factor = weight / (weight + distance);
    
    return Vector(factor, factor, factor);
}
